#!/usr/bin/env Rscript

# Save pipeline label objects for barrier follow-up scripts.
# Run this in the same R session, or after loading an .rdata that contains
# agree_dat and/or sd3_new_af2.

args <- commandArgs(trailingOnly = TRUE)
outfile <- if (length(args) >= 1 && nzchar(args[[1]])) args[[1]] else "barrier_label_objects.rdata"

objs_to_save <- intersect(c("agree_dat", "sd3_new_af2"), ls(envir = .GlobalEnv, all.names = TRUE))
if (length(objs_to_save) == 0L) {
  stop(
    "No label objects found in the global environment.\n",
    "Expected at least one of: agree_dat, sd3_new_af2",
    call. = FALSE
  )
}

save(list = objs_to_save, file = outfile, envir = .GlobalEnv)
message("Saved label objects to: ", normalizePath(outfile))
message("Objects saved: ", paste(objs_to_save, collapse = ", "))
