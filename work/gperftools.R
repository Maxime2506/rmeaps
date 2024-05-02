# Set compiler variables, so we'll find the conda-based gperftools install
Sys.setenv(
  CPATH=sprintf("%s:%s/include", Sys.getenv("CPATH"), Sys.getenv("CONDA_DIR")),
  LIBRARY_PATH=sprintf("%s:%s/lib", Sys.getenv("CPATH"), Sys.getenv("CONDA_DIR"))
)
devtools::install_github("bnprks/Rgperftools")
