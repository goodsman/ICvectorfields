# currently only contains functions to ensure proper clean up of compiled code
# when package is unloaded.

.onUnload <- function (libpath) {
  library.dynam.unload("ICvectorfields", libpath)
}
