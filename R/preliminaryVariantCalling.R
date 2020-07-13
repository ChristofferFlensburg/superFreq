
performPreliminaryVariantCallingOnMissingVCFs = function(metaData, reference, cpus=1) {
  #check which VCF files are missing, place them somewhere smart, in plots or R?
  missing = !file.exists(metaData$VCF)
  if ( !any(missing) ) return()

  #check if varscan exists
  catLog('Check if the varscan command exists.\n')
  a = system('varscan')
  if ( a != 0 ) {
    catLog('Cant find varscan command, so cant perform preliminary variant calling of missing VCF.\n')
    warning('Cant find varscan command, so cant perform preliminary variant calling of missing VCF.')
    return()
  }
  catLog('Seems like varscan command exists.\n')
  catLog('Samtools shouldve been check earlier, so assuming that exists.\n')
  
  #use the one-liner from the superFreq README
  systemCalls = paste0('samtools mpileup -d 1000 -q 15 -Q 15 -A -f ', reference, ' ',  metaData$BAM[missing], ' | varscan mpileup2cns - --variants --strand-filter 0 --p-value 0.5 --min-var-freq 0.02 --output-vcf 1 > ', metaData$VCF[missing])
  ret = mclapply(systemCalls, function(call) {
    catLog(call, '\n\n')
    ret = system(call)
    return(ret)
  }, mc.cores=cpus)
  if ( all(ret == 0) ) catLog('Runs returned with exit status 0, seems good.\n')
  else {
    catLog('Not all runs returned with exit status 0. Might be a problem, but continuing with a warning.\n')
    warning('Not all varscan runs exited with status 0. Preliminary variant calling might not have been done properly. Exit statuses: ', ret)
  }

  return()
}
