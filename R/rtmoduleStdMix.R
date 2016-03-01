

findFragments = function(fragList, mzrt, refRt, refMz, rtWin=7, mzDiff=0.01) {
	# Detect list of fragments in the given peak list.
	# Return the fragment indices, or zero if not found.
	frags = unlist( strsplit(fragList, split="[,;]") )
#	# round m/z values to 2 decimal digits
	if (sum (floor(as.numeric(frags))==floor(refMz))) {
		# remove parent ion from the list of fragments
		frags = frags[-(floor(as.numeric(frags))==floor(refMz))]
	}
	allMz = floor(mzrt$mz*100)/100
	frags.found = vector()
	for (i in 1:length (frags)) {
		if( length(grep("\\(", frags[i])) ) {
			fragment = as.numeric(strsplit(frags[i],"[()]")[[1]][1])
		} else {
			fragment = as.numeric(frags[i])
		}		
		mzMatch = which( fragment>allMz-mzDiff & fragment<allMz+mzDiff )
		if ( length(mzMatch) ) {
			mzrtMatch = which((refRt >= mzrt$rt[mzMatch]-rtWin) & (refRt <= mzrt$rt[mzMatch]+rtWin))	
			if ( length(mzrtMatch) ) {
				frags.found = c(frags.found,mzMatch[mzrtMatch])
			} else {
				frags.found = c(frags.found,0)	
			}
		}
	}
	return( frags.found )
}

# findParentPeak(pl,refTable[i,],refMz,ms2.parent,intCol1,intCol2,massTol)

findParentPeak = function (
  pl,
  refData,
  refMz,
  intCol1,
  intCol2,
  massTol,
  rtTol=120) { 
  
	# detect RT tags in the given peak list. 
  refRT = refData["rt.lib"] * 60 
  rt.min = refRT - rtTol 
  rt.max = refRT + rtTol 
   
  # get sequence of low/high energy columns  
  lowE_idx = seq (intCol1,(intCol2-1),by=2)
  highE_idx = seq ((intCol1+1),intCol2,by=2)
  # case of one sample
  if (length(lowE_idx)==1 & length(highE_idx)==1) {
    mzrt = data.frame ( 
      mz = pl$mz, 
      rt = pl$rt,
      lowE = pl[,lowE_idx],
      highE = pl[,highE_idx]
    )
  } else {
  # case of PL containing multiple samples - use mean intensity values  
    mzrt = data.frame ( 
      mz = pl$mz, 
      rt = pl$rt,
      lowE = rowMeans(pl[,lowE_idx],na.rm=TRUE),
      highE = rowMeans(pl[,highE_idx],na.rm=TRUE)
    )
  }
  mzrt$massTol = round(massTol*mzrt$mz/1e6,4)
	output = list()
	mzMatch = which ( (refMz > (mzrt$mz - mzrt$massTol)) & (refMz < (mzrt$mz + mzrt$massTol)) )
#browser()  
# check mass match
	if ( length(mzMatch) ) {
	# check mass & RT match for all mass-matched peaks
	mzrtMatch = which( (mzrt[mzMatch,"rt"] > rt.min) & (mzrt[mzMatch,"rt"] < rt.max) )
	matchedPeaks = mzrt[mzMatch[mzrtMatch],]
	if ( length(mzrtMatch) ) {		
		# check if lowE channel intensity is not lower than the highE channel intensity
		chChk = which (matchedPeaks$lowE > matchedPeaks$highE) 
		if( !length(chChk) ) {
				output = NA						
		} else {
			matchedPeaks = matchedPeaks[chChk,]
		# if there are several peaks which correspond to all conditions -         
			for( i in 1:nrow(matchedPeaks) ) {
			  # if there are MS/MS - find matching fragments 
			  if( !is.na(refData["msms.frgs"]) ) {
			    fragment.ions = findFragments(as.character(refData["msms.frgs"]),mzrt,matchedPeaks$rt[i],refMz)
			  } else { fragment.ions=0 }
			  if (any (fragment.ions != 0)) {			    
					output[[i]] = c(
						# store index of main peak
						as.numeric(row.names(matchedPeaks[i,])),	 
						# store indices of fragment peaks 
						fragment.ions[which(fragment.ions != 0)] 
					)
			 } else {
          # no fragments = no annotation
					  output[[i]] = NULL
			}
		}
	}   
	} # end: if ( length(mzrtMatch) )
	}
#  print (output)
	out.len = sapply (output, length)
	if (!(length(output))) {
    # remove selections without detected fragment(s)!
		output = NA
	} else 
  if (length (output)>1) {
		# if all selection have the same number of detected fragment -
		# then select the one with higher intensity values.
		if (all(out.len)==out.len[which.max(out.len)]) {
			int = sapply (output, function(i) mzrt$lowE[i])
			output = output[which.max (int)]
		} else {
		# otherwise return the selection with most detected fragments
			output = output[which.max (out.len)]
		}
	}
  # return just ID of primary ion 
	return(output[[1]][1])	
}
  
#########################################################################
# Main RT module 							                                          #
# Input: reads a table of reference compounds 				                  #		
# (containing chemical formula,retention time & msms mass fragments) 	  #
# detects them in the input peak list. It then builds a robust linear 	#
# model of detected vs. reference RT, removes outliers and builds an 	  #
# improved new model. It then shifts RT of all peak list features and	  #
# calculates RT deviation limits on the 95% confidence limit of the 	  #
# model residuals. 							                                        #
# Output: shifted RT model data with CIs for the mean 	                #
#########################################################################

rtModule = function (plList,ionMode,refTable,project,massTol,confint=0.95,massErrorSurface=NULL) {

	minRT = 55
  
###########################################
# Get input and set variables
###########################################
	matched.log = data.frame ()
	matched = 0
	peaks.out = list()	# output
###########################################
# Find RT tags:
###########################################
	print ("###############")
	print ("Getting RT tags")
	print ("###############")
	for ( i in 1:nrow(refTable) ) { # 
	# get theoretical mass values using the input formula!
	if (is.na (refTable[i,"formula"])) { next() }
	tryCatch ({
		mol = Rdisop::getMolecule(refTable[i,"formula"])
	}, error = function(err) {
		print (paste("Formula",err))
		stop ("Reference table error - check chemical formulaS")
	} )
	M = round (mol$exactmass,4)
	if(is.na(M)) { next() } # skip compound - no theoretical mass!
  if ("msms.parent" %in% names(refTable)) {
	  ms2.parent = refTable$msms.parent[i]
  } else { ms2.parent=NA }
	# Set the reference mass value, depanding on ion mode and adducts
	if ( ionMode == "positive" ) {
	# protonation
	  refMz = M + protonMass
	  mainAdduct = M + sodiumMass 	      
	} else if ( ionMode == "negative" ) {
	# deprotonation
	  refMz = M - protonMass
	  mainAdduct = M + formicAcidMass 
	} else {
		stop  (paste("\nError:",ionMode))
	}
#	refRT = floor (refTable[i,"lib.rt"]*60)		
#	refRT = floor (refTable[i,"rt"]*60)
	compound = as.character (refTable[i,"compound"])
	print (paste("Searching:",compound))
	print (paste("Referece M/Z:",refMz))
  # read peak lists:
  matchedList = lapply (plList, function(pl) {
    
    intCol1 = min (grep(".*0[1|2]$", names(pl)))
    intCol2 = max (grep(".*0[1|2]$", names(pl)))  

    peakData = findParentPeak(pl,refTable[i,],refMz,intCol1,intCol2,massTol)
#browser()     
	  if (is.na(peakData)) {	
      # Try annotating using main adduct
	    peakData = findParentPeak(pl,refTable[i,],mainAdduct,intCol1,intCol2,massTol)	
	  }
  	if (!is.na(peakData)) {    
 		  cbind (
			"comp.idx"=i,
      "name"=compound,
			"rt.lib"=round(refTable[i,"rt.lib"]*60),
			"rt.samp"=round(pl[peakData[[1]],"rt"]),
			"rt.diff"=round(pl[peakData[[1]],"rt"]-refTable[i,"rt.lib"]*60)
 		  )}   
  })

	matchedList = do.call('rbind',matchedList) 
 
	if (!is.null(matchedList)) {
		print (paste("Found:",nrow(matchedList),"RT tags for compound:",compound))	
		matched = matched + nrow(matchedList)	
		matched.log = data.frame(rbind(matched.log,matchedList))
	}	
  matched.log$rt.lib = as.numeric(matched.log$rt.lib)
  matched.log$rt.samp = as.numeric(matched.log$rt.samp)
}	# end ( i in 1:nrow(refTable.superpools) )
#	names (matched.log) = c("comp.idx","rt.lib","rt.samp","peak.idx")

###########################################
# Generate model, short RT values
###########################################
	print ("#########################")
  print (paste("Found total:",matched,"RT-markers!"))
  print ("Generating RT shift model")
	print ("#########################")	

	maxRt = max (matched.log$rt.lib,matched.log$rt.samp,na.rm=TRUE)

	if (nrow(matched.log)) {
  	# Build a robust LM model:
  	rt.shift.model = lm (rt.samp~rt.lib, data=matched.log)
  	# Define a confidence interval for RT search
  	CI = mad (residuals(rt.shift.model))*1.46*qt(confint,df=(nrow(matched.log)-2)) 
    # Shift RT tags RT (use as a control for plotting - later apply to whole peak list)
  	matched.log$rt.shifted = matched.log$rt.samp*1/rt.shift.model$coefficients[2]-rt.shift.model$coefficients[1]
  	matched.log$rt.shifted[matched.log$rt.shifted<minRT] = minRT
  	# Write a report of the RT tags
  	write.table (matched.log,file=paste("matched.log",project,ionMode,Sys.Date(),"tsv",sep="."),row.names=FALSE,sep="\t")
  	inSpan = length (which (
  		matched.log$rt.shifted<=(matched.log$rt.lib+round(CI)) & matched.log$rt.shifted>=(matched.log$rt.lib-round(CI))
  		))
  	png (width=800,height=800,file=paste("rt.shift",project,ionMode,Sys.Date(),"png",sep="."))
  	plot (matched.log$rt.lib, matched.log$rt.lib,xlim=c(0,maxRt),ylim=c(0,maxRt),cex=0.7,xlab="RT.lib",ylab="RT.samp")
  	abline (0,1)
  	# Plot corrected data points for RT tags 		
  	points (matched.log$rt.lib, matched.log$rt.shifted,cex=0.7,col=2)
  	points (matched.log$rt.lib, matched.log$rt.samp,pch='+')
  	# Add CIs
  	lines (matched.log$rt.lib, matched.log$rt.lib+CI,col=3,lty=3)
  	lines (matched.log$rt.lib, matched.log$rt.lib-CI,col=3,lty=3)
  	text (0,maxRt,labels=paste("Model offset:",round(rt.shift.model$coefficients[1])," sec"),col=2,pos=4)
  	text (0,maxRt-40,labels=paste("CI =",round(CI),"sec"),col=2,pos=4)
  	text (0,maxRt-80,labels=paste("Observations:",nrow(matched.log)),col=2,pos=4)		
  	text (0,maxRt-120,labels=paste("Observations in CI span:",inSpan),col=2,pos=4)		
  	dev.off()	
  	model = list("model"=rt.shift.model,"CI"=CI)
  	print ("###############")
  	print("Done: RT Module")
  	print ("###############")
	}
 	return (model)
}
