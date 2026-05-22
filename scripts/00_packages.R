# =============================================================================
# Project:  PROTECTS-Microbiome
# Script:   00_packages.R
# Analyst:  M. Mangalea
# Date:     2026-02-25
# Purpose:  Install and load all packages required for PROTECTS-Microbiome analyses.
#           Run this script once after opening the project on a new R session
#           to confirm all packages are loaded.
# =============================================================================

# -- 0. Confirm R version -----------------------------------------------------
stopifnot(
  "R >= 4.5.0 required" = getRversion() >= "4.5.0"
)
cat ("R version:", R.version$version.string, "\n")

# Use a personal library first to avoid permission issues on shared systems.
user_lib <- Sys.getenv(
  "R_LIBS_USER",
  unset = "~/R/x86_64-pc-linux-gnu-library/4.5"
)

dir.create(path.expand(user_lib), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(path.expand(user_lib), .libPaths()))

cat("Active library paths:\n")
print(.libPaths())

# --1. Helper Functions -------------------------------------------------------

# Remove stale package lock directories left by interrupted installs.
clean_stale_locks <- function(lib_path = .libPaths()[1]) {
  lock_dirs <- list.files(
    path = lib_path,
    pattern = "^00LOCK",
    full.names = TRUE
  )
  
  if (length(lock_dirs) > 0) {
    message(
      "Removing stale lock directories: ",
      paste(basename(lock_dirs), collapse = ", ")
    )
    unlink(lock_dirs, recursive = TRUE, force = TRUE)
  }
}

# Install packages from CRAN.
install_cran_packages <- function(pkgs,
                                  repos = "https://cran.rstudio.com/",
                                  lib = .libPaths()[1]) {
  if (length(pkgs) == 0) {
    return(invisible(NULL))
  }
  
  clean_stale_locks(lib)
  
  install.packages(
    pkgs,
    repos = repos,
    lib = lib,
    dependencies = TRUE
  )
}

# Install optional CRAN packages; continue if they fail.
install_optional_cran <- function(pkgs,
                                  repos = "https://cran.rstudio.com/",
                                  lib = .libPaths()[1]) {
  if (length(pkgs) == 0) {
    return(invisible(NULL))
  }
  
  clean_stale_locks(lib)
  
  for (pkg in pkgs) {
    if (pkg %in% rownames(installed.packages(lib.loc = .libPaths()))) {
      next
    }
    
    ok <- tryCatch({
      install.packages(
        pkg,
        repos = repos,
        lib = lib,
        dependencies = TRUE
      )
      pkg %in% rownames(installed.packages(lib.loc = .libPaths()))
    }, error = function(e) {
      message("Optional package '", pkg, "' failed: ", conditionMessage(e))
      FALSE
    })
    
    if (!ok) {
      warning(
        "Optional package '", pkg,
        "' could not be installed. Continuing without it."
      )
    }
  }
  
  invisible(NULL)
}

# Load packages with a clear error if any are missing.
load_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(
        "Package '", pkg, "' is not installed or not available on .libPaths().",
        call. = FALSE
      )
    }
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
  invisible(NULL)
}


# --2. Load CRAN packages -----------------------------------------------------

cran_pkgs <- c(
  "tidyverse",
  "ggrepel",
  "ggsignif",
  "RColorBrewer",
  "readxl",
  "vegan",
  "Polychrome",
  "ape",
  "ggalluvial",
  "circlize",
  "viridis",
  "fmsb",
  "forcats",
  "viridisLite",
  "corrplot",
  "Nonpareil",
  "knitr"
)

# Optional package with a compiled dependency chain that may fail on systems
# with `cannot find -lgfortran` errors. Analyses should still run without it.
optional_cran_pkgs <- c("rbiom")

installed_now <- rownames(installed.packages(lib.loc = .libPaths()))
to_install_cran <- setdiff(cran_pkgs, installed_now)


if (length(to_install_cran) > 0) {
  message(
    "Installing missing CRAN packages: ",
    paste(to_install_cran, collapse = ", ")
  )
  install_cran_packages(to_install_cran, lib = .libPaths()[1])
} else {
  message("All required CRAN packages already installed.")
}

install_optional_cran(optional_cran_pkgs, lib = .libPaths()[1])

load_packages(cran_pkgs)

# -- 3. Bioconductor setup -----------------------------------------------------

bioc_pkgs <- c(
  "ComplexHeatmap",
  "phyloseq",
  "microbiome",
  "Maaslin2",
  "EnhancedVolcano"
)

installed_now <- rownames(installed.packages(lib.loc = .libPaths()))
to_install_bioc <- setdiff(bioc_pkgs, installed_now)

if (length(to_install_bioc) > 0) {
  message(
    "Installing missing Bioconductor packages: ",
    paste(to_install_bioc, collapse = ", ")
  )
  
  clean_stale_locks(.libPaths()[1])
  
  BiocManager::install(
    pkgs = to_install_bioc,
    ask = FALSE,
    update = FALSE,
    lib = .libPaths()[1]
  )
} else {
  message("All required Bioconductor packages already installed.")
}

load_packages(bioc_pkgs)

# -- 5. Session info -----------------------------------------------------------

si <- sessionInfo()
print(si)

# Optional reproducibility log:
# dir.create("results", showWarnings = FALSE)
# writeLines(capture.output(si), con = file.path("results", "session_info.txt"))

message("\n00_packages.R complete — required packages installed and loaded.")


