
#' internal function for debugging purposes.
#'
#' @details This function saves the input to a function in a directory, so that the exact input triggering the bug
#'          can be examined.
dumpInput = function(Rdirectory, inputList) {
  catLog('Dumping input for debugging purposes..')
  debugDir = paste0(Rdirectory, '/debug')
  if ( !file.exists(debugDir) ) dir.create(debugDir)
  dumpDir = paste0(debugDir, '/inputDump-', make.names(Sys.time()), '-', abs(rnorm(1, 0, 0.1)))
  dir.create(dumpDir)
  for ( input in names(inputList) ) {
    catLog(input, '..', sep='')
    assign(input, inputList[[input]])
    save(list=input, file=paste0(dumpDir, '/', input, '.Rdata'))
  }
  catLog('done.\n')
}
