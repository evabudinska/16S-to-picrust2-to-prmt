# 11_kegg_info.R — KEGG utilities (KO↔Rxn map + robust REST fetch s fallbackem)
suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -- KO↔Reaction map --------------------------------------------------------

ensure_kos_rxns <- function(path, overwrite = FALSE) {
  if (!overwrite && file.exists(path)) return(invisible(normalizePath(path)))
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)

  url <- "https://rest.kegg.jp/link/reaction/ko"
  con <- url(url); on.exit(try(close(con), silent = TRUE), add = TRUE)
  x <- readLines(con, warn = FALSE)
  if (!length(x)) stop("Empty response from KEGG link endpoint: ", url)

  dt <- fread(text = paste(x, collapse = "\n"), header = FALSE, sep = "\t",
              colClasses = "character")
  if (ncol(dt) != 2L) stop("Unexpected format from KEGG link endpoint")
  setnames(dt, c("KO","Rxn"))
  dt[, KO  := sub("^ko:", "", KO)]
  dt[, Rxn := sub("^rn:", "", Rxn)]
  dt <- unique(dt)

  fwrite(dt, path, sep = "\t")
  if (nrow(dt) < 1000) warning("Unusually small KO↔Reaction mapping: ", nrow(dt))
  invisible(normalizePath(path))
}

# -- Flatfile parser (ENTRY / EQUATION / ORTHOLOGY) -------------------------

.parse_kegg_flat <- function(lines) {
  if (!length(lines)) return(NULL)
  key <- trimws(substr(lines, 1, 12), which = "right")
  val <- substr(lines, 13, 1e6)

  end_idx <- which(trimws(lines) == "///")
  if (length(end_idx) == 0L) end_idx <- length(lines)
  entry_idx <- which(key == "ENTRY")
  if (!length(entry_idx)) return(NULL)

  out <- vector("list", length(entry_idx))
  for (i in seq_along(entry_idx)) {
    start <- entry_idx[i]
    stopi <- if (i < length(entry_idx)) entry_idx[i+1] - 1L else end_idx[length(end_idx)]
    k <- key[start:stopi]
    v <- val[start:stopi]

    acc <- list(); lastk <- NULL
    for (j in seq_along(k)) {
      kj <- k[j]; vj <- v[j]
      if (nzchar(kj)) { lastk <- kj; acc[[lastk]] <- c(acc[[lastk]], vj) }
      else if (!is.null(lastk)) acc[[lastk]] <- c(acc[[lastk]], vj)
    }

    entry_raw <- paste(acc$ENTRY %||% "", collapse = " ")
    rn <- sub(".*\\b(R[0-9]{5})\\b.*", "\\1", entry_raw)

    equation <- trimws(paste(acc$EQUATION %||% "", collapse = " "))

    orth_lines <- acc$ORTHOLOGY %||% character()
    kos <- unique(unlist(regmatches(orth_lines, gregexpr("K[0-9]{5}", orth_lines))))
    orth <- if (length(kos)) stats::setNames(kos, kos) else character()

    out[[i]] <- list(ENTRY = paste0(rn, "  Reaction"), EQUATION = equation, ORTHOLOGY = orth)
    names(out)[i] <- rn
  }
  out
}

# -- HTTP get helper (base R; bez dalších závislostí) -----------------------

.http_get_lines <- function(url, timeout_sec = 120) {
  o <- getOption("timeout"); on.exit(options(timeout = o), add = TRUE)
  options(timeout = max(getOption("timeout"), timeout_sec))
  con <- url(url); on.exit(try(close(con), silent = TRUE), add = TRUE)
  readLines(con, warn = FALSE)
}

# -- Adaptivní fetch: batch -> půlení -> single; checkpoint + fail log ------

.fetch_reactions <- function(rxn_ids,
                             chunk_size = 3,
                             sleep_sec = 1,
                             long_sleep_sec = 60,
                             max_retries = 6,
                             checkpoint_file = NULL,
                             existing = NULL) {

  rxn_ids <- unique(rxn_ids)
  if (is.null(existing)) existing <- list()
  have <- names(existing)
  todo <- setdiff(rxn_ids, have)
  failed <- character()

  # resume failed list z checkpointu (pokud existuje)
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    cp <- try(readRDS(checkpoint_file), silent = TRUE)
    if (!inherits(cp, "try-error") && is.list(cp)) {
      if (!is.null(cp$Reaction_info) && is.list(cp$Reaction_info)) {
        existing <- cp$Reaction_info
        have <- names(existing)
      }
      if (!is.null(cp$failed)) failed <- unique(c(failed, cp$failed))
    }
  }

  # helper: ulož checkpoint
  .save_cp <- function() {
    if (!is.null(checkpoint_file)) {
      saveRDS(list(Reaction_info = existing, failed = failed), checkpoint_file)
    }
  }

  # helper: vlož placeholder pro chybějící reakce (kvůli downstreamu)
  .ensure_placeholders <- function(ids) {
    for (id in ids) {
      if (is.null(existing[[id]])) {
        existing[[id]] <- list(ENTRY = paste0(id, "  Reaction"),
                               EQUATION = NULL, ORTHOLOGY = character())
      }
    }
  }

  # rekurzivní fetch bloku ids
  fetch_block <- function(ids, cs) {
    ids <- unique(ids)
    if (!length(ids)) return(invisible())

    # zkus batch fetch s několika pokusy
    for (attempt in seq_len(max_retries)) {
      url <- paste0("https://rest.kegg.jp/get/", paste(ids, collapse = "+"))
      res <- try(.http_get_lines(url), silent = TRUE)

      good <- (!inherits(res, "try-error")) && length(res) && any(grepl("^ENTRY", res))
      if (good) {
        parsed <- .parse_kegg_flat(res)
        if (!is.null(parsed) && length(parsed)) {
          existing[names(parsed)] <- parsed
          # doplň placeholdery pro ID, která server ve stejné odpovědi vynechal
          missing <- setdiff(ids, names(parsed))
          if (length(missing)) {
            failed <<- unique(c(failed, missing))
            .ensure_placeholders(missing)
          }
          .save_cp()
          Sys.sleep(sleep_sec)
          return(invisible())
        }
      }

      # backoff
      if (attempt < max_retries) {
        if (!inherits(res, "try-error") && grepl("ENTRY", paste(res, collapse = "\n"))) {
          Sys.sleep(max(2, sleep_sec * attempt))
        } else {
          message(sprintf("KEGG fetch error for [%s] (attempt %d/%d); backing off %ds",
                          paste(head(ids,3), collapse = ","),
                          attempt, max_retries, long_sleep_sec))
          Sys.sleep(long_sleep_sec)
        }
      }
    } # retries

    # batch selhal -> pokud jde rozdělit, rozděl; jinak jednostlivě
    if (length(ids) > 1) {
      mid <- floor(length(ids)/2)
      fetch_block(ids[seq_len(mid)], max(1, cs %/% 2))
      fetch_block(ids[(mid+1):length(ids)], max(1, cs %/% 2))
    } else {
      # single ID s per-ID retriem
      id <- ids[1]
      ok <- FALSE
      for (attempt in seq_len(max_retries)) {
        url <- paste0("https://rest.kegg.jp/get/", id)
        res <- try(.http_get_lines(url), silent = TRUE)
        good <- (!inherits(res, "try-error")) && length(res) && any(grepl("^ENTRY", res))
        if (good) {
          parsed <- .parse_kegg_flat(res)
          if (!is.null(parsed) && length(parsed)) {
            existing[id] <- parsed[id]
            ok <- TRUE; break
          }
        }
        if (attempt < max_retries) {
          Sys.sleep(if (grepl("Forbidden|403", paste(res, collapse = " "), ignore.case = TRUE))
                      long_sleep_sec else max(2, sleep_sec * attempt))
        }
      }
      if (!ok) {
        failed <<- unique(c(failed, id))
        .ensure_placeholders(id)
        .save_cp()
      }
    }
  }

  # hlavní smyčka přes počáteční chunking
  idx <- ceiling(seq_along(todo) / max(1, chunk_size))
  chunks <- split(todo, idx)
  for (i in seq_along(chunks)) {
    fetch_block(chunks[[i]], chunk_size)
  }

  # finální checkpoint + fail log (vedle checkpointu)
  .save_cp()
  if (length(failed)) {
    flog <- if (is.null(checkpoint_file)) "failed_rxns.txt"
            else file.path(dirname(checkpoint_file), "failed_rxns.txt")
    writeLines(sort(unique(failed)), flog)
    message("KEGG fetch finished with ", length(failed), " failed IDs. Logged: ", flog)
  } else {
    message("KEGG fetch finished with 0 failed IDs.")
  }

  existing
}

# -- Public builder ----------------------------------------------------------

get_kegg_reaction_info <- function(kos_to_rxns_file,
                                   kolist = NULL,
                                   save_out = TRUE,
                                   out_prefix = "",
                                   checkpoint_file = NULL,
                                   chunk_size = 3,
                                   sleep_sec = 1,
                                   long_sleep_sec = 60,
                                   max_retries = 6) {
  stopifnot(data.table::is.data.table(kos_to_rxns_file))
  if (!identical(names(kos_to_rxns_file), c("KO","Rxn"))) {
    setnames(kos_to_rxns_file, c("KO","Rxn"))
  }

  all_kegg <- list(
    KOs       = sort(unique(kos_to_rxns_file$KO)),
    Reactions = sort(unique(kos_to_rxns_file$Rxn))
  )

  if (!is.null(kolist) && length(kolist)) {
    kos_to_rxns_file <- kos_to_rxns_file[KO %in% kolist]
    all_kegg$KOs       <- intersect(all_kegg$KOs, kolist)
    all_kegg$Reactions <- intersect(all_kegg$Reactions, unique(kos_to_rxns_file$Rxn))
  }

  # načti případný checkpoint (kvůli resume)
  existing <- list()
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    cp <- try(readRDS(checkpoint_file), silent = TRUE)
    if (!inherits(cp, "try-error") && is.list(cp) && !is.null(cp$Reaction_info)) {
      existing <- cp$Reaction_info
      if (!is.list(existing)) existing <- list()
      message("Resuming KEGG fetch from checkpoint: ", checkpoint_file,
              " (", length(existing), " entries cached)")
    }
  }

  fetched <- .fetch_reactions(
    rxn_ids         = all_kegg$Reactions,
    chunk_size      = chunk_size,
    sleep_sec       = sleep_sec,
    long_sleep_sec  = long_sleep_sec,
    max_retries     = max_retries,
    checkpoint_file = checkpoint_file,
    existing        = existing
  )

  # srovnej pořadí a doplň placeholdery na jistotu
  missing_all <- setdiff(all_kegg$Reactions, names(fetched))
  if (length(missing_all)) {
    for (id in missing_all) {
      fetched[[id]] <- list(ENTRY = paste0(id, "  Reaction"),
                            EQUATION = NULL, ORTHOLOGY = character())
    }
  }
  all_kegg$Reaction_info <- fetched[all_kegg$Reactions]
  message("Done downloading reaction info: ", length(all_kegg$Reaction_info))

  all_kegg$kos_to_rxns <- kos_to_rxns_file

  if (isTRUE(save_out)) {
    fn <- paste0(out_prefix, "_KeggReactions.rda")
    dir.create(dirname(fn), showWarnings = FALSE, recursive = TRUE)
    save(all_kegg, file = fn)
    message("Saved all_kegg -> ", fn)
  }
  all_kegg
}

