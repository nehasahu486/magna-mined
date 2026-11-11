# Drop-in helper for MagnaGenes v2.0 (data-raw/prepare_data.R)
# Automatically selects the latest dated `fin_data_YYYY-MM-DD.tsv` in
# `data-raw/source/annotation/`, moves older files into `archive/`, and
# deduplicates by MD5 hash so `archive/` never accumulates duplicate files.

resolve_and_archive_fin_data <- function(ann_dir) {
  files <- list.files(ann_dir, pattern = "^fin_data(_[0-9]{4}-[0-9]{2}-[0-9]{2})?(_[0-9a-f]{8})?\\.tsv$", full.names = TRUE)
  if (!length(files)) stop("No fin_data*.tsv found in ", ann_dir)

  # Extract YYYY-MM-DD date from filename (undated fin_data.tsv sorts as 0000-00-00)
  extract_date <- function(p) {
    m <- regmatches(basename(p), regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", basename(p)))
    if (length(m) == 1 && nzchar(m)) m else "0000-00-00"
  }
  dates <- vapply(files, extract_date, character(1))
  mtimes <- file.info(files)$mtime
  ord <- order(dates, mtimes, decreasing = TRUE)
  files <- files[ord]

  latest <- files[1]
  older <- files[-1]

  if (length(older)) {
    arch_dir <- file.path(ann_dir, "archive")
    dir.create(arch_dir, recursive = TRUE, showWarnings = FALSE)
    existing_arch <- list.files(arch_dir, pattern = "\\.tsv$", full.names = TRUE)
    known_hashes <- if (length(existing_arch)) unname(tools::md5sum(existing_arch)) else character(0)

    for (f in older) {
      h <- unname(tools::md5sum(f))
      if (h %in% known_hashes) {
        unlink(f)  # Skip duplicate already in archive/
        next
      }
      dest <- file.path(arch_dir, basename(f))
      if (file.exists(dest)) {
        stem <- sub("\\.tsv$", "", basename(f))
        dest <- file.path(arch_dir, sprintf("%s_%s.tsv", stem, substr(h, 1, 8)))
      }
      file.rename(f, dest)
      known_hashes <- c(known_hashes, h)
    }
  }

  latest
}

# Usage inside data-raw/prepare_data.R (replaces hardcoded fin_data.tsv line):
# fin_data_path <- resolve_and_archive_fin_data(file.path(raw_dir, "annotation"))
# message("Using annotation file: ", basename(fin_data_path))
# ann <- fread(fin_data_path, sep = "\t", quote = "", encoding = "UTF-8")
