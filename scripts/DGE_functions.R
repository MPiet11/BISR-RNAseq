# Automate package installation -------------------------------------------

if ("BiocManager" %in% rownames(installed.packages()) == FALSE) {
  install.packages("BiocManager")
}

ProjectLibraries <- function(pkgs) {
  library(BiocManager)

  check <-
    sapply(pkgs,
           require,
           warn.conflicts = TRUE,
           character.only = TRUE)
  if (any(!check)) {
    message("Loading missing packages for the project")
    pkgs.missing <- pkgs[!check]
    BiocManager::install(pkgs.missing, update = TRUE, ask = FALSE)
    check <-
      sapply(pkgs.missing,
             require,
             warn.conflicts = TRUE,
             character.only = TRUE)
  } else {
    print(check)
    message("All project packages are loaded")
  }

}
