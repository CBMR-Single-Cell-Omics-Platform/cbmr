.onUnload <- function (libpath) {
  library.dynam.unload("cbmr", libpath)
}