#' Build chemical standards based MS library.
#' 
#' @param stdInfo A string, the path plus file name directing to a pool information file.
#' @param stdFiles A string, the path to a directory containing the peak list files corresponding with sample pools.
#' @param polarity A string, the MS ionization mode.
#' @return Returns the MS library object of chemical standards as specified in the input file.
#' @examples\donttest{
#' stdInfo <- system.file("extdata/samplePool_neg.tsv",package="matchWeiz")
#' stdFiles <- system.file("extdata/pool_data",package="matchWeiz")
#' testlib.neg <- buildMSlib (stdInfo,stdFiles,"negative")
#' stdInfo <- system.file("extdata/samplePool_pos.tsv",package="matchWeiz")
#' testlib.pos <- buildMSlib (stdInfo,stdFiles,"positive") 
#' }
#' @export
buildMSlib = function(stdInfo,
                      stdFiles,    
                      polarity                  
  ) {
  
  require(Rdisop)
  
  tryCatch({
    info = read.delim(stdInfo,as.is=TRUE)
  }, error=function(err) {
    stop ("Failed to read the stdInfo. file:", conditionMessage(err))
  }) 
  
  pools = unique(info$stdFile)  
  tryCatch ({
    stdData = lapply (file.path(stdFiles,pools),read.delim,as.is=TRUE) 
  }, error=function(err) {
    stop ("Failed to read peak lists of chemical standards:", conditionMessage(err))
  })

  sampleNames = unlist(lapply(pools,function(x) {strsplit(x,split=".(csv|tsv|txt)")[[1]]}))
  names(stdData) = sampleNames

  mainOut = data.frame()
  matched = 0
  total = 0

  for (i in 1:length(pools)) {
    poolData = info[info$stdFile==pools[i],]
    pl = stdData[[i]]
    print (paste("Processing pool:",i))
    results = apply(poolData,1,function(x) { 
      print( paste("Cehcking pool compound:",x["ID"]) )
      findCompound(x,pl,polarity) 
    } )
    matched = matched + sum(!unlist(sapply(results,is.null)))
    total = total + nrow(poolData)    

    poolResults = mapply ( function(x,i,file)
    {	
      if (!is.null(x)) {      
        peaks.out = convertResults(x, file, polarity, poolData[i,])
        return (peaks.out)
      } 
    }, results, c(1:nrow(poolData)), SIMPLIFY=FALSE )
  } # end: main loop
  
  mainOut = rbind (mainOut, do.call ("rbind", c(poolResults,deparse.level=1) ))
  # write output results without the name of the row names  
  row.names (mainOut) = NULL
  #write.table (mainOut,file=paste("msLib",polarity,".tsv",sep="_"),sep="\t",row.names=F)
  MSlib = buildLibobject (mainOut, polarity)
  #save (MSlib,file=paste("MSlib",polarity,"RData",sep="."))
  return (MSlib)
} # end of function 'buildMBlib'


