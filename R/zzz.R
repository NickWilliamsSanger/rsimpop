.onLoad = function(libname, pkgname) {
  initSimPop(-1,bForce = TRUE)
}

.onUnload =function(libpath) {
  library.dynam.unload("rsimpop", libpath)
}
