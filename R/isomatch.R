#################
# IsoMatch module
#################

getIsotopes = function (peaks, intCol) {
	# Expect a sub set of the xcms peak list with Camera annotations
	output = list()
	i = 0
	# Get all peaks with isotope IDs
	isoIds = grep ("^\\[\\d+\\]",peaks[,"isotopes"])
	if (length(isoIds)<2) { return (NULL) }
	peaks = peaks[isoIds,]
	# Recursively extract all relevent groups
	while (nrow(peaks)>1) {
		# highest intensity ion peak
		theMain = which.max (peaks[,intCol])
		# Take all isotope peaks within the same camera group 
 		theGroup = which (peaks[,"pcgroup"]==peaks[theMain,"pcgroup"])
		if (length(theGroup)>1) {
			intDiff = diff (peaks[theMain:max(theGroup),intCol])
			mzDiff = diff (peaks[theMain:max(theGroup),"mz"])
			realIsos = theMain
			# test differences in intensity and mass
			for (j in seq(along.with=intDiff)) {
			# add ions with decreasing intensity and mass
				if (intDiff[j]>0) { break() }
				if (round (mzDiff[j],0) != 1) { break() }				
				realIsos = c(realIsos,(theMain+j))
			}
      # add identified isotopes to output  
			if (length(realIsos)>1) {	
				i = i+1;				
				output[[i]] = peaks[realIsos,]	}
			# remove the identified isotopes from peak list
			peaks = peaks[-realIsos,]
		} else {
			peaks = peaks[-theGroup,]
		}
	} 
	if (!length (output)) { output = NULL }
	return (output)
}

getPossibleIon = function (rules,peak,neutralIonMass) {

	# Get peaks possibly related to the main ion formula	
	#  by mass differences specified by the 'rules' table.
	# Return the valid peaks, with associated formula annotation.
	
	# all valid mass diffrences from the neutral mass ion
	massDiffs = rules$nmol*neutralIonMass + rules$massdiff
	# all valid peaks
	ion = character()
  
	massMatch = which(abs(massDiffs-peak$mz) < peak$massTol)
	if (length(massMatch)) {
			ion = as.character(rules$name[massMatch])
	}	  
	return (rules[rules$name==ion,])
}

parseCameraAdducts = function (rule, baseFormula) {
	# parse and convert the Camera package adduct formula into a chemical formula.
	# The Camera formula prototype:
	# 	"["[Mcount]"M"([adductType][adductSimbol])*"][chargeCount]*"-"
  # E.g.: [M-2H+NH4]-
  print (paste(rule$name,"and",baseFormula))
	output = character()
  #	charge = rule$charge
	Mcount = rule$nmol	
	for (i in 1:Mcount) { output = paste(output,baseFormula,sep="") }
  #	Get the charge number
	formula = rule$name
	# Get the adduct type and symbol
	adduct = sub (x=formula, pattern="\\[.*M(.*)\\].*$", replacement="\\1")
  
	for (i in 1:(length (strsplit (adduct,split="[A-Z][a-z]?\\d*")[[1]]))) {
		# Get adduct type
		adductType = sub (x=adduct, pattern="^([-|+]).*", replacement="\\1")
		# Get the first element
		adductSymb = sub (pattern="^([-|+]\\d?[A-Z][a-z]?\\d*).*$", replacement="\\1",x=adduct)
		# Get the rest
		adduct = sub (adductSymb,"",adduct, fixed=TRUE) # strsplit(adduct, split=adductSymb)[[1]][2]
		# Add the sign
		if (!length(grep("^[-|+]",adduct))) { adduct = paste (adductType,adduct,sep="") }
		if (adductType=='-') { 
			if (adductSymb == "-H") { adductNum=1 } 
			else if (adductSymb == "-2H") { adductNum=2 } 				
			else if (adductSymb == "-3H") { adductNum=3 } 
			else {
				# look for "[notDigit][digit]" pattern
				adductNum = sub (x=adductSymb,pattern="^-\\D+(\\d*)$",replacement="\\1")
			}
			if (adductNum=="") {
				adductNum = 1
			}
			adductSymb = sub (x=adductSymb,pattern="-\\d?(\\D+)\\d*$",replacement="\\1")
			adductSymb = paste ("(",adductSymb,"-",adductNum,")",sep="")
		}
		adductSymb = sub (x=adductSymb, pattern="^[+]", replacement="")
		output = paste(output,adductSymb,sep="")
	}
	output = CHNOSZ::as.chemical.formula(CHNOSZ::makeup(output))
  print (paste("Converted to:",output))
	return (output)	
}

isoMatchFilter = function (decompObj,lmModel,iso2ratio,ci.level=0.95) {

##########################################################
# Get: the observed 2nd/1st isotope ratio, 
#   RDISOP object, prediction model and required CI level. 
# Does:
# 1) Generate all possible formulas using 'Rdisop'
# 2) Remove non-'CH' compounds
# 3) Get nHydrogens,nCarbons for each formula
# 4) Use a model to predict iso2/iso1 ratio + intervals
# 5) Apply the isotope ratio filter
# 6) Apply the nHydrogens/nCarbons 'golden rule' filter
# 7) Return all valid formulas
##########################################################

	CH_compounds = decompObj$formula[grep("^C\\d*H\\d*", decompObj$formula)]
	CH_scores = decompObj$score[grep("^C\\d*H\\d*", decompObj$formula)]
	print (paste ("Found",length(CH_compounds),"CH-compounds"))
	nC = as.numeric(gsub("^C(\\d*).*","\\1",CH_compounds))
	nC[is.na(nC)] = 1
	nH = as.numeric(gsub("^C.*H(\\d*).*","\\1",CH_compounds))
	nH[is.na(nH)] = 1
	iso2ratio.predict = predict( lmModel, newdata= data.frame(
			nH=nH,
			nC=nC
		), interval="prediction", level=ci.level )

	iso2filter = (iso2ratio>iso2ratio.predict[,"lwr"]) & (iso2ratio<iso2ratio.predict[,"upr"])
	HC_ratioFilter = (nH/nC > 0.2) & (nH/nC < 3.1)

	return (cbind (
		validFormulas = CH_compounds[iso2filter & HC_ratioFilter],
		validScores = CH_scores[iso2filter & HC_ratioFilter]
		))
}

run_isoMatchModule = function (pl.list,intCol,baseMass,baseFormula,adductRules,lmModel,ci.level=0.95) {

## 	run the IsoMatchModule.
##	Input: 
##	pl - list of identified isotopic peaks associated with a single ion specie
## 	intCol - the pl column holding relevant intensity information
##	lmModel - the LM model describing the module's linear model
##	baseMass - measured or theoretical mass of main ion specie
##	baseFormula - theoretical chemical formula of main ion specie
##	ci.level - confidence limit for the application of the LM model

	elements.many = Rdisop::initializeElements(c("C","H","N","O","P","S","Na"))
	elements.few = Rdisop::initializeElements(c("C","H","N","O")) 

	print ("###############") 
	print ("IsoMatch Module")
	print ("###############")

  #	The emprical 99% CI level for the true formula occurance (in a random DB pool of 327 compounds) as in:
  #	quantile (testIt.44[,4],0.99,na.rm=TRUE) = 163
  # sd.analyticon.isotopes = sd (testIt.44[,4], na.rm=T) = 47.33
  # Hence we set the maximum allowed rank for a chemical formula:
	HIT_LIMIT = 163	

	## Iterate over the groups of input peaks 
	isoMatches = sapply (pl.list, function(peaks) {	
	# decompose isotopic peaks
      
    # get the possible formula, given adduct rules and mass difference to main ion
    possibleIon = getPossibleIon (adductRules,peaks[1,],baseMass)
    if (nrow(possibleIon)) {
      # Convert the given (CAMERA!) adduct notation to valid chemical formula
      foundFormula = parseCameraAdducts (possibleIon, baseFormula)
      print (paste("Analyzing formula:",foundFormula))
      # Do decomposition and filter detected isotopes
    	massTol = max (round ((peaks["massTol"]),4))
      # mass vector
    	mass = as.numeric (peaks[,"mz"])
      # intensity vector
    	int = as.numeric (peaks[,intCol])
        if (mass[1] < 800) {
          # Decompose isotopes using 'Rdisop' and get the decomposed object
    		  decompObj = Rdisop::decomposeIsotopes (mass,int,mzabs=massTol,maxisotopes=4,elements=elements.many)	# 
        } else {
          # Decompose isotopes using 'Rdisop' and get the decomposed object			
	        decompObj = Rdisop::decomposeIsotopes (mass,int,ppm=5,maxisotopes=4,elements=elements.few) # 
        }
        # the rank before filters
        matchUnfiltered = match(foundFormula,decompObj$formula)
    		# Compute the relative intensity of [M+1] isotope to [M] (used for filtering)
        iso2ratio = int[2]/int[1]*100	
    		# Apply the isoMatch filter on possible chemical formulas
    	  filtered = isoMatchFilter (decompObj,lmModel,iso2ratio)
    	  if (nrow(filtered)) {
     		  print(paste(length(filtered),"filtered compounds"))
      		  ## Check if any of the valid chemical formulas matches
      		  if (!is.na(match(foundFormula,filtered))) 
		        {
			               print (paste("Score of found formula by rank:",match(foundFormula,filtered)))
                	   isoScore = round(fdrtool::phalfnorm (match(foundFormula,filtered),fdrtool::sd2theta(sd=47), lower.tail=FALSE),2)
#	                   return(isoScore)
        	  } 
#            else {
#                return (NULL)
#        	  }
        	}
    	}
	})
return (isoMatches)
}
