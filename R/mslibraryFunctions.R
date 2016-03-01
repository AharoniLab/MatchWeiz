#####################
# common adducts		#
#####################

# protonMass = 1.007278	# [M-H]- or [M+H]+
# sodiumMass = 22.9892	# [M+Na]+
# potassiumMass = 38.9632	# [M+K]+
# formicAcidMass = 44.9982 	# [M-H+HCOOH]-
# oxygenMass = 16.000748 	# [M-2H+H2O]

#######################
# other global vars		#
#######################

# MAX_RANK = 15 

getAtoms = function(formula) {
	# parse chemical formula and return the atomic composition
	output = vector()
	numOfAtoms = length(strsplit(formula, split="([A-Z][a-z]*\\d*)")[[1]])
	for ( i in 1:(numOfAtoms) ) {
		elem = gsub(pattern=".*([A-Z][a-z]*).*", replacement="\\1", formula)
		formula = strsplit(formula, split=elem)[[1]][1]
		output = c(output,elem)
	}
	return(output)
}


# parseCameraAdducts = function ( rule, mainIon ) {
# 	# parse the Camera adduct formula into relevant symbols:
# 	# 1. Number of main ion specie
# 	# 2a. The adduct symbol
# 	# 2b. The adduct type (adduct or fragment)
# 	# 3. The number of charges
# 	# The formula prototype:
# 	#	[M-2H+NH4]-
# 	# 	"["[Mcount]"M"([adductType][adductSimbol])*"][chargeCount]*"-"
# 
# 	# Here we go:
# 	output = character()
# #	charge = rule$charge
# 	Mcount = rule$nmol	
# 	for (i in 1:Mcount) { output = paste(output,mainIon,sep="") }
# 
# #print (output)
# #
# #	# Get the charge number
# 	formula = rule$name
# 	# Get the adduct type and symbol
# 	adduct = sub (x=formula, pattern="\\[.*M(.*)\\].*$", replacement="\\1")
# #print (adduct)
# 	for (i in 1:(length (strsplit (adduct,split="[A-Z][a-z]?\\d*")[[1]]))) {
# print (i)
# 		# Get adduct type
# 		adductType = sub (x=adduct, pattern="^([-|+]).*", replacement="\\1")
# print (adductType)
# 		# Remove adduct type chars from string
# #		adduct = sub (x=adduct, pattern="^[-|+]", replacement="")
# print (adduct)
# 		# Get the first element
# 		adductSymb = sub (pattern="^([-|+]\\d?[A-Z][a-z]?\\d*).*$", replacement="\\1",x=adduct)
# print (adductSymb)
# #print (adduct)		
# 		# Get the rest
# 		adduct = sub (adductSymb,"",adduct, fixed=TRUE) # strsplit(adduct, split=adductSymb)[[1]][2]
# 		# Add the sign
# 		if (!length(grep("^[-|+]",adduct))) { adduct = paste (adductType,adduct,sep="") }
# #		if (!length(grep("^[-|+]",adduct)) & adductType=='-') { adduct = paste (adductType,adduct,sep="") }
# #		else if (length(grep("^[+]",adduct))) { adduct = sub (x=adduct, pattern="^[-|+]", replacement="") }
# print (adduct)		
# #stop()
# 		if (adductType=='-') { 
# print (adductSymb)
# 			if (adductSymb == "-H") { adductNum=1 } 
# 			else if (adductSymb == "-2H") { adductNum=2 } 				
# 			else if (adductSymb == "-3H") { adductNum=3 } 
# 			else {
# 				# look for "[notDigit][digit]" pattern
# 				adductNum = sub (x=adductSymb,pattern="^-\\D+(\\d*)$",replacement="\\1")
# 			}
# 			if (adductNum=="") {
# #				switch (adductSymb, "-H" = 1, "-2H" = 2, "-3H" = 3)
# 				adductNum = 1
# #				 
# 			}
# print (adductNum)
# #
# 			adductSymb = sub (x=adductSymb,pattern="-\\d?(\\D+)\\d*$",replacement="\\1")
# 			adductSymb = paste ("(",adductSymb,"-",adductNum,")",sep="")
# print (adductSymb)
# 		}
# 		adductSymb = sub (x=adductSymb, pattern="^[+]", replacement="")
# 		output = paste(output,adductSymb,sep="")
# 	}
# 
# 	print (paste("Converted to:",output))
# 	#output = as.chemical.formula(makeup (output))
# 	output = CHNOSZ::as.chemical.formula(CHNOSZ::makeup(output))
# 	print (output) 
# 	return (output)	
# }


# peaks.out = convertResults(x, file, ionMode, poolData[i,])
convertResults = function(peakData,mol,polarity,stdInfo) {
  
	lowE = 1
	highE = 2
	intVec = grep("0[1|2]$",names(peakData[[1]]))
	peaks.out = data.frame()	# final output structure
	mainPeak = peakData[[1]]
	mainPeakTag = names(peakData)[1]
	peaks.out = rbind(peaks.out, cbind(
				"ID" = stdInfo["ID"],
				"name" = stdInfo["compound"],
				"peak.tag" = mainPeakTag,
				"formula" = stdInfo["formula"],
				"M" = peakData[["MM"]],
				"rt.win" = stdInfo["rtWin"], 
				"XCMS" = "$$$", 	# data separator
				"row.idx" = as.numeric(rownames(mainPeak)),
				"mz" = mainPeak$mz, 
				"rt.med" = round(mainPeak$rt/60,2),
				"rt.min" = round(mainPeak$rtmin/60,2),
				"rt.max" = round(mainPeak$rtmax/60,2),
				"lowE" = floor(as.numeric(mainPeak[intVec[lowE]])),
				"highE" = floor(as.numeric(mainPeak[intVec[highE]])),
				"isotopes" = mainPeak$isotopes,
				mainPeak["pcgroup"],
				"peak.table" = stdInfo$stdFile)
			)
	
	isotopes = peakData$Isotopes
	if ( !is.null(isotopes) ) {				
		for( isoNum in 1:length(isotopes) ) { 
				peakTag = names(isotopes)[isoNum]
				data = data.frame( rep(as.data.frame(""), length.out=ncol(peaks.out)) )
				names( data ) = names(peaks.out)
				data$ID = stdInfo$ID
				data$peak.tag = peakTag
				data$XCMS = "$$$"
				data$rt.med = round(isotopes[[isoNum]]$rt/60, digits=2)
				data$rt.min = round(isotopes[[isoNum]]$rtmin/60,2)
				data$rt.max = round(isotopes[[isoNum]]$rtmax/60,2)
				data$mz = isotopes[[isoNum]]$mz	
				data$row.idx = as.numeric(rownames(peakData$Isotopes[[isoNum]]))
				data$pcgroup = as.numeric(peakData$Isotopes[[isoNum]]$pcgroup)
				data$peak.table = stdInfo$stdFile
				data$lowE = floor(as.numeric(isotopes[[isoNum]][intVec[lowE]]))
				data$highE = floor(as.numeric(isotopes[[isoNum]][intVec[highE]]))
				data$isotopes = isotopes[[isoNum]]$isotopes
				# bind this data to the output table
				peaks.out = rbind(peaks.out, data)
			}
		} 

	adducts = peakData$Adducts
	if( length(adducts) ) {
		for( adductIdx in 1:length(adducts) ) 
		{
			adductName = names(peakData[["Adducts"]])[adductIdx]
			# get initial data headers from the table of standards
			data = data.frame( rep(as.data.frame(""), length.out=ncol(peaks.out)) )
			names( data ) = names(peaks.out)
			data$ID = stdInfo$ID
			data$peak.tag = adductName	# peak class tag 
			data$XCMS = "$$$"
			data$mz = as.numeric(adducts[[adductIdx]]$mz)
			data$rt.med = round(adducts[[adductIdx]]$rt/60, 2) # xcms RT 
			data$rt.min = round(adducts[[adductIdx]]$rtmin/60, 2) # xcms RT 			
			data$rt.max = round(adducts[[adductIdx]]$rtmax/60, 2) # xcms RT 
			data$row.idx = as.numeric(rownames(adducts[[adductIdx]]))
			data$pcgroup = as.numeric(adducts[[adductIdx]]$pcgroup)
			data$peak.table = stdInfo$stdFile
			data$lowE = floor(as.numeric(adducts[[adductIdx]][intVec[lowE]]))
			data$highE = floor(as.numeric(adducts[[adductIdx]][intVec[highE]]))
			data$isotopes = adducts[[adductIdx]]$isotopes
			peaks.out = rbind(peaks.out, cbind(data))
		}
	}
 
	ionCluster = peakData$pcgrp
	if (length(ionCluster)) {
		for( ionIdx in 1:nrow(ionCluster) ) 
		{
			# add ion cluster group data
			data = data.frame( rep(as.data.frame(""), length.out=ncol(peaks.out)) )
			names( data ) = names(peaks.out)
			data$ID = stdInfo$ID
			data$peak.tag = "pcgrp"	
			data$XCMS = "$$$"
			data$rt.med = round(ionCluster$rt[ionIdx]/60,2)
			data$rt.min = round(ionCluster$rtmin[ionIdx]/60,2)
			data$rt.max = round(ionCluster$rtmax[ionIdx]/60,2)
			data$mz = ionCluster$mz[ionIdx]
			data$row.idx = as.numeric(rownames(ionCluster[ionIdx,]))
			data$pcgroup = as.numeric(ionCluster$pcgroup[ionIdx])
			data$peak.table = stdInfo$stdFile
			data$lowE = floor(as.numeric(ionCluster[ionIdx,intVec[lowE]]))
			data$highE = floor(as.numeric(ionCluster[ionIdx,intVec[highE]]))
			data$isotopes = ionCluster[ionIdx,"isotopes"]
			# bind this data to the output table
			peaks.out = rbind(peaks.out, data)
		}
	}
#browser()
	rownames(peaks.out) = NULL
	return(peaks.out)
}	

findCompound = function(poolData,pl,polarity) 
{
##############################################################
# An intermediate function to extract data of pooled samples #
# Pass data to a generic peak-picking function 'peak detect' #
##############################################################

##
## 	Estimate RT range by RT window number
##
if (is.na(poolData["RT"])) {
  	rtWin = as.integer(poolData["rtWin"]) 
  	rtBins = seq(1,39,by=2)
  	estimatedRtRange = c(rtBins[max(1,rtWin-2)], rtBins[min(20,rtWin+3)]) * 60
  	print ("Estimating RT range by RT bins")
  	refRT = (rtWin*2-1)*60
} else {
    refRT = poolData["RT"]*60 # manually observed RT in minutes
    estimatedRtRange = c(refRT-30,refRT+30)
    print ("Estimating RT range by observed RT")
}  
##	
## 	Calculate exact mass using the 'Rdisop' package:
##		check for the cases where the compound is a natural cation in solution
##		in these cases for positive mode: mass = [M], and for negative mode: M = [M-2H]
##

  chemFormula = poolData["formula"]

	if (length(grep("\\+$",chemFormula))) { 
#		print(poolData["formula"])
#browser()
		# remove trailing plus sign 
	#	poolData["chemFormula"] = substring(poolData["chemFormula"], 1,nchar(poolData["chemFormula"])-1)		
		# set mode
		naturalCation = 1 
	} else { 
		naturalCation = 0 
	}
	elem = c(Rdisop::initializeCHNOPS(),Rdisop::initializeElements("Cl"),Rdisop::initializeCharges())	
	mol = Rdisop::getMolecule(chemFormula, elements=elem)
	Mneutral = mol$exactmass # get exact neutral mass
	##
	## NOTE: generalize these function
	if (polarity=="negative") {
		if (naturalCation) {
			# get number of Hydrogens in original formula, then substract two (double deprotonation)
			Hnum = as.numeric( sub( ".*H(\\d*).*","\\1",chemFormula) ) - 2			
		} else {
			# get number of Hydrogens in original formula, then substract one (deprotonation)
			Hnum = as.numeric( sub( ".*H(\\d*).*","\\1",chemFormula) ) - 1
		}
	} else if (polarity=="positive") {
		if (naturalCation) {
			# no need to protonate
			Hnum = 0
		} else {
			# get number of Hydrogens in original formula, then add one (protonation)
			Hnum = as.numeric( sub( ".*H(\\d*).*","\\1",chemFormula) ) + 1
			print(Hnum)			
		}
	}
  print(paste("Chem Formula:",chemFormula))
  refFormula = sub("H(\\d*)",paste("H",Hnum,sep=""),chemFormula) # adjust chemical formula
	print(paste("Ref Formula:",refFormula))

	output = peakDetect_2ch(
		pl=pl,
		Mneutral=Mneutral,
		refRT=refRT,
		lowBound=estimatedRtRange[1], 
		highBound=estimatedRtRange[2],
		polarity=polarity,
		chemFormula=chemFormula,
		refFormula = refFormula,
		naturalCation=naturalCation
	)
  output = c(output,"MM"=Mneutral)
	return (output)
}

peakDetect_2ch = function(pl,
                          Mneutral,
                          refRT,
                          lowBound,
                          highBound,
                          polarity,
                          chemFormula,
                          refFormula,
                          naturalCation,
                          refErr=20,
                          spec=NULL) 
{
	#########################################################################
	# Detect peak row in xcms peak table, based on following procedures: 	#
	# 1) Mass accuracy range (M/z Da.)					#
	# 2) Estimated RT and RT boundaries (time set in seconds)		#
	# 3) Peak signal to noise ratio (XCMS calculated S/N) 			#
	# 4) Camera peaks group - minimum is two members			#
	# 5) Main peak isotopes - minimum is one extra (M+1)			#
	# If none or multiple peaks are detected - return NULL			#
	# * Add a module which returns wild Mz matches across all chromatogram	#
	#########################################################################

  print(paste(refRT,lowBound,highBound))
#browser()
	Mp = NA 
	intCol = grep("01$",names(pl))	# get the low energy intensity column idx
##
## Set ion mode dependant parameters, including reference mass
##
  if (polarity=="positive") { 
		if (naturalCation) { 
        Mtag="[M]+" 
        refMz = Mneutral
		} else {
        Mtag = "[M+H]+"
        refMz = Mneutral + protonMass
		} 
		ions = ruleset.pos
	} else if (polarity=="negative") { 
		if (naturalCation) { 
        Mtag="[M-2H]-" 
        refMz = Mneutral - 2*protonMass
		} else { 
        Mtag="[M-H]-" 
        refMz = Mneutral - protonMass
		} 
		ions = ruleset.neg
	}
  massIons = Mneutral + ions$massdiff 
  names(massIons) = ions$name
##
## 	Mass2mass matching
##
  massDelta = round(refMz*refErr/1e+06,digits=4)	# convert PPM window to mass units (Da)
  mzMatch = which( (refMz > pl[,"mz"]-massDelta) & (refMz < pl[,"mz"]+massDelta) )  # Mz match by mass accuracy
	
	if (!length(mzMatch)) {
    # look for another 
	 massMatches = lapply(massIons, function(refMz) { 
		mzMatch = which((refMz > pl[,"mz"]-massDelta) & (refMz < pl[,"mz"]+massDelta) ) 
		if (length(mzMatch)) {
		  mzrtMatch = which ((pl[mzMatch,"rt"] < highBound) & (pl[mzMatch,"rt"] > lowBound))	
		  if (length(mzrtMatch)) {       
  
		    intValid = pl[mzMatch[mzrtMatch],intCol] > pl[mzMatch[mzrtMatch],(intCol+1)]
		    mzrtMatch = mzrtMatch[intValid]		
		    mainIon = mzMatch[mzrtMatch[which.max(pl[mzMatch[mzrtMatch],intCol])]]
		    if (length(mainIon)) {
			pc = pl$pcgroup[mainIon]
			pcgrp = pl[pl$pcgroup==pc,]
			if (nrow(pcgrp)>1) {
			  if (pl[mainIon,intCol] < max(pcgrp[,intCol])) {
			 	return (NULL)
			  } else {			  
			     	return (mainIon)
			  }
			} else {
				return (NULL)
			}
		      } else {
		  	  return (NULL)	
		      }
		    }  
		  }
	    })	   	 		

		  massMatches = do.call ("rbind",massMatches)

		  if (is.null(massMatches)) { return (NULL) }
    
		  if (nrow(massMatches)>1) { massMatches = as.matrix(massMatches[which.max(pl[massMatches,intCol]),]) }  
    
    mainIon = pl[massMatches,]
		Mtag = rownames(massMatches)
		# get the ion rule
		rule = ions[ions$name==Mtag,]
		# adjust formula
		refFormula = parseCameraAdducts (rule, chemFormula)	
  #browser() 
	} else {
	  ## main parent MZ match:
	  mzrtMatch = which ((pl[mzMatch,"rt"] < highBound) & (pl[mzMatch,"rt"] > lowBound))	
	  if (length(mzrtMatch)) { 
	    intValid = pl[mzMatch[mzrtMatch],intCol] > pl[mzMatch[mzrtMatch],(intCol+1)]
	    mzrtMatch = mzrtMatch[intValid]		
      if (length(mzrtMatch)) {
	      mainIon = pl[mzMatch[mzrtMatch[which.max(pl[mzMatch[mzrtMatch],intCol])]],]
      } else { return(NULL) }
	  } else {
      return (NULL)
	  }
	}

	pc = mainIon$pcgroup
	pcgrp = pl[pl$pcgroup==pc,]
  maxPeak = max(pcgrp[,intCol])
	print(mainIon)
	isoNum = gsub("^\\[(\\d+).*$","\\1",mainIon$isotopes)
	isotopes = grep(paste("\\[",isoNum,"\\]",sep=""), pl[,"isotopes"])
	##
	## 	Check Isotopes
	##
  if (length(isotopes)) {
	  iso = chkIsotopes(mainIon, pl[isotopes,], polarity, refFormula, intCol, Mneutral)
  } else {
    iso = NULL
  }

		##
		## 	Get Camera adducts and grouping		
		##
			# if pcgroup is small and there are no isotopes - return NULL
			if ((nrow(pcgrp)==1) & (is.null(iso))) { 
				print("Warning: pcgroup too small. No Isotopes Found... Returning NULL")	
        
  #>>> might add here an accuracy + Mtag filter to allow single peak enteries - as we saw some FNs 
  # (but add only when ppm<=5 and Mtag=[parent ion])
        
 				return(NULL)
			} else if ( nrow(pcgrp) > 1 ) {
				if (!is.null(iso)) {
					# remove found isotopes from peak groups
					exclude = which(row.names(pcgrp) %in% as.numeric(lapply(iso,rownames)))
					if (length(exclude)) {
						pcgrp = pcgrp[-exclude,] 
					}
				}
				# try to detect some main ion adducts and attach the pcgrp peaks
				adducts = attachAdducts(mainIon, pcgrp, Mneutral)
				# redefine pcgrp, excluding annotated adducts
				pcgrp = adducts[["pcgrp"]]
				# define separatly the annotated adducts, otherwise as NULL
				if ( length(which(names(adducts)!="pcgrp")) ) {
					adducts = adducts[which(names(adducts)!="pcgrp")]
				} else {
					adducts = NULL
				}
			} else {
				# there are isotopes, but no adducts 
				adducts = NULL
			}
##
## 	Remove main ion peak from pcgrp 
##
#print(pcgrp) length(pcgrp)
			if (!is.null(pcgrp)) {
			# remove the main ion from the pcgrp
				if ( length(which(row.names(pcgrp) == mainIon)) ) {
					pcgrp = pcgrp[-which(row.names(pcgrp) == mainIon),]					
				}	

    # apply relative intensity filter
				intVec = pmax(pcgrp[,intCol],pcgrp[,(intCol+1)])
				lowAdducts =  which (intVec < 0.05*maxPeak & pcgrp$mz > Mneutral)
				lowFrags =  which (intVec < 0.01*maxPeak & pcgrp$mz < (Mneutral-8))
				toRemove = c(lowFrags,lowAdducts)
				if (length(toRemove)) {
				   if (length(toRemove)) { pcgrp = pcgrp[-toRemove,] }
				}
				if (!nrow(pcgrp)) { pcgrp=NULL }
			}	                 					
			output = list(mainIon,				
				            iso,
				            adducts,
				            pcgrp
                )
			names (output) = c(Mtag,"Isotopes","Adducts","pcgrp")
			return (output)		

	# no proper peak matching found
	print("No xcms peak match or multiple matches or main ion not found - returning NULL!")
	return (NULL)					 
} 


chkIsotopes = function(mainIonPeak,peaks,polarity,chemFormula=NA,intCol=NA,M,refErr=30) {
	# use the Rdisop package 
	# get the calculates theoretical isotopes
	# find the corresponding peaks in xcms peak list
	# decompose found peaks
	# compare highest scoring (most probable) chemical furmula to the given formula.
	# return a list of annotated isotope peaks.	
	if (is.na(chemFormula)) { return(NULL) }
	print("==========================")
	print("Cheking possible isotopes:")
	print("==========================")
	isotopes = list()
	minIso = 2
	maxIso = 4
	rtBounds = 4
	## use Rdisop to extract isotopes
		# remove protonation sign, if exists
	if (length(grep("\\+$",chemFormula))) {
	  chemFormula = substring(chemFormula, 1,nchar(chemFormula)-1)
	}
	# define the elements for isotope decomposition
  t = getAtoms(as.character(chemFormula))
	essentialElements = Rdisop::initializeElements (rev(t))		
#	print( sapply(essentialElements, function(x) { x$name }) )	
	# adjust molecular formula
	if (polarity == "positive") { z = 1 } else if (polarity == "negative") { z = -1 }
	print(paste("Using:",chemFormula))
	# create Rdisop object
	mol = Rdisop::getMolecule(chemFormula, elements=essentialElements, z=z)
  # theoretical mass values calculated by formula
	isotopicMasses = mol$isotopes[[1]][1,1:maxIso]
	if (isotopicMasses[1] > 900) { decompPPM = 5 } else { decompPPM = 30 }

	# set the first iso peak as main ion
	isotopeMatch = data.frame ( 
			"mz" = peaks$mz,
			"int" = as.numeric(peaks[,intCol]),
			"peakIdx" = as.numeric(row.names(peaks))
	)

	# prepare vectors of detected masses, intensities and corresponding peak indices
	masses = isotopeMatch[,1]
	intensities = isotopeMatch[,2]
	peakIds = isotopeMatch[,3]

	# decompose detected signals and check that chemical formula is possible
	iso = Rdisop::decomposeIsotopes(masses, intensities, ppm=decompPPM, maxisotopes=4, elements=essentialElements)

  if( length(grep(chemFormula,getFormula(iso))) ) {
			for (i in 2:length(masses)) {
				M_Annot = paste("M+",i-1,sep="")
				isotopes[[M_Annot]] = peaks[i,]
			}
			print("**** Isotopes Check - OK ****")
			return(isotopes)
	} 
	return(NULL)
}

attachAdducts = function(mainIon, peaks, Mneutral, adduct=NA) {
	print("==========================")
	print("Checking possible adducts:")
	print("==========================")
	adducts = list()		
	# exclude peak of the main ion from the ion cluster
	peaks = peaks[-which(row.names(peaks)==row.names(mainIon)),]
	# find the manually detected adduct or principle ion, if exists
	if (!is.na(adduct) ) {
		mzMatch = which( round(peaks$mz,1) == round(detectedAdduct,1) )
		if (length(mzMatch)) {
		# add this adduct as manually detected feature
			adducts[["DA"]] = peaks[mzMatch,]
		# remove the adduct from the ion cluster
			peaks = peaks[-mzMatch,]
		}			
	}
	refMass = floor(Mneutral)
	mzMatch = grep(refMass, peaks$adduct)
	# get the Camera adducts that match the ref mass
	if (length(mzMatch) ) {	
	  print("**** Added adducts ****")
		adductNames = unlist( lapply( peaks$adduct[mzMatch], function(x) { 
			adductName = strsplit(as.character(x), split=" ")[[1]][1]
		}) )
		for (i in 1:length(mzMatch)) {
			adducts[[adductNames[i]]] = peaks[mzMatch[i],]
		}		
		peaks = peaks[-mzMatch,]
	} else {
	  print("**** Not found ****")
	}
	if (nrow(peaks)>0) { 
		# annotated all remaining peaks as "pcgrp" members
		adducts[["pcgrp"]] = peaks	
	}
	return (adducts)
}


buildLibobject = function(DBpeaks, polarity, massTol=0.01) {

  
  print ("Building library object from table of peaks")
  Libobj = list(peaks=DBpeaks, polarity=polarity)
  names(Libobj$peaks)[names(Libobj$peaks)=="rt.med"] = "rt.db"
  names(Libobj$peaks)[names(Libobj$peaks)=="mz"] = "mz.db"
  maxInt.vec = unlist( apply(Libobj$peaks[,c("lowE","highE")], 1, function(x) { max(x) } ) )
  maxInt.vec[maxInt.vec==0] = NA
  ##	get mass error tolerance by automatic calculation using mass and intensity - TO ADD
  ### Libobj$peaks$massTol = calcMassErr(Libobj$peaks$mz, maxInt.vec, massErrorSurface)
  ## set mass tolerance by the given threshold
  Libobj$peaks$massTol = massTol
  ##	Calculate mass error tolerance for the main ions
  print ("Calculating mass error tolerances...")
  annotatedPeaks = which( !is.na( as.numeric(Libobj$peaks$M) ) )
  Libobj$peaks$massTol[annotatedPeaks] = apply( Libobj$peaks[annotatedPeaks,], 1, function(peak) {
    M = as.numeric(peak["M"])
    mz = as.numeric(peak["mz.db"])	
    # calc mass diff for the neutral specie
    if (floor(M) == floor(mz)) { abs(M - mz) }
    # calc mass diff for the M-H specie
    else if (floor(M-protonMass) == floor(mz)) { abs(M-protonMass-mz) }
    # calc mass diff for the M+H specie
    else if (floor(M+protonMass) == floor(mz)) { abs(M+protonMass-mz) }
    else { massTol }    
  }) 
  ##	rank peaks high energy ramp channel (DB putative fragments)
  print ("Ranking fragments...")
  dbIds = unique(Libobj$peaks$ID)
  Libobj$ranked = lapply(dbIds, function(id) {
    lib = Libobj$peaks[Libobj$peaks$ID==id,]
    frgs = lib[lib$mz < (as.numeric(lib$M[1])-4*protonMass),]
    # remove any CAMERA annotated isotopes
    whichIso = grep("M\\+",frgs$isotopes[1:nrow(lib)])
    frgs = frgs[-whichIso,c("mz.db","lowE","highE")]
    names(frgs) = c("mz","lowE","highE")
    # add maxInt column (for int. threshold filtering and order peaks by the highE col.
    if (nrow(frgs)) {
      frgs$max = pmax(frgs[,2],frgs[,3])
      frgs[order(frgs$highE, decreasing=T),][1:min(MAX_RANK,nrow(frgs)),]
    }
  })
#  browser()  
  names(Libobj$ranked) = dbIds
  rownames(Libobj$peaks) = NULL
  print ("Done: library object")
  return ( Libobj )	
}