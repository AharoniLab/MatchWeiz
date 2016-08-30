############################
##	MM Utility functions 	##
############################

getAtoms = function(formula) {
	# parse chemical formula and return the atomic composition
	output = vector()
	numOfAtoms = length (strsplit(formula, split="([A-Z][a-z]*\\d*)")[[1]])
	for ( i in 1:(numOfAtoms) ) {
		elem = gsub (pattern=".*([A-Z][a-z]*).*",replacement="\\1",formula)
		formula = strsplit (formula, split=elem)[[1]][1]
		output = c(output,elem)
	}
	return(output)
}

adjustFormula = function (chemFormula, z) {
	# protonate/deprotonate, remove Na and charge symbols from chemical formula
	if ( z == 1 ) {
			# set molecule charge 
			
			# check if molecule is already protonated
			if (length(grep("\\+$",chemFormula))) {
				# remove cherge sign 
				chemFormula = sub("\\+$","",chemFormula)
				# get number of Hydrogens in original formula	
				Hnum = as.numeric( sub( ".*H(\\d+).*","\\1",chemFormula) )
			} else if (length(grep("Na",chemFormula))) {
				# replace Sodium adducts with Hydrogens
				Na_num = as.numeric( sub( ".*Na(\\d*).*","\\1",chemFormula) )
				if (is.na(Na_num)) { Na_num = 1 }
				chemFormula = sub ( "Na(\\d*).*","",chemFormula)
				Hnum = as.numeric ( sub( ".*H(\\d+).*","\\1",chemFormula) ) # regsub + 1 + Na_num
			} else {
				# protonate according to the number of Hydrogens in original formula
				Hnum = as.numeric( sub( ".*H(\\d+).*","\\1",chemFormula) ) + 1
			}
	} else if ( z == -1 ) {
		# get number of Hydrogens in original formula, then deprotonate
		if (length(grep("\\+$",chemFormula))) {
			# remove cherge sign 
			chemFormula = sub("\\+$","",chemFormula)
			# positively charged compounds need to be double deprotonated!
			Hnum = as.numeric( sub( ".*H(\\d+).*","\\1",chemFormula) ) - 2
		} else {
			Hnum = as.numeric( sub( ".*H(\\d*).*","\\1",chemFormula) ) - 1
		}
	}
	chemFormula = sub("H(\\d*)",paste("H",Hnum,sep=""),chemFormula)
	return( chemFormula )
}


annotateFeature = function(
	    peak,       	# A vector of mass, rt, and mass tolerance
	    polarity,   	# polarity 
      DBobj,		# DB object
      rtTol      	# allowed rt tol in seconds
) {               

  ## 	Annotate a single peak against the selected database
  ## 	Returns the feature index in the DB or NULL if not found

  ## 	calculate mass-rt differences
  deltam = abs(DBobj$peaks[,"mz.db"] - peak["mz"])
  deltart = abs(DBobj$peaks[,"rt.db"]*60 - peak["rtShifted"])
  ## 	combining the mass tolerances
  ## ...We assume the library mass values has some variability which is added to the peak mass variability
  idm = deltam < abs(DBobj$peaks$massTol + peak["massTol"])
  idt = deltart <= rtTol   
  output <- which(idm & idt)

	if ( !length(output) )  { output = NULL } 

	return (output)
}


getIds = function(
	        pl.annot, # mz&rt matches of pl peaks
      		DBobj,    # library+mass error surface+ranked fragments
      		pl,       # input peak list
	        intCol1,   # pl column of the first intensity values
	        intCol2,   # pl column of the last intensity values
  		    polarity,   # ionization mode
		      rtModel
		) {
    
    	count=TRUE
    	isoMatch=TRUE
    	main=TRUE
    	xrank=TRUE

  	##
  	## load adduct rules and scoring models
  	## 
    
  	data (rulesets)
    if (polarity=="positive") { adductRules = ruleset.pos} else { adductRules = ruleset.neg }
  	data(isoMatchModel)	
  	data(xrankScoremat) # fragMatch
    
  	# get indices for intensity columns - depends on existance of prefix (01) or (02) for MS channels!
  	colsMs1 = grep(".*01$", names(pl))+3
  	colsMs2 = grep(".*02$", names(pl))+3
  	# get all maximum MS1 intensity columns across peak table rows: 
    maxIntCols = grep(".*01$", names(pl))[max.col(pl[,grep(".*01$", names(pl))])]
    # add extra columns accounting for additinal fields
    maxIntCols = maxIntCols+3

	  # gets annotation list object (output of 'annotateFeature')
	  # converts annotations (by DB IDs) into a list sorted by pcgroups (peak clusters)
    annot.table = list()
	  # Relate peak features to DB peaks position
    for (i in seq(along = pl.annot)){
    	if ( !is.null(pl.annot[[i]]) ) { # annot.table[[i]] = NA } else {
		  annot.table[[i]] = cbind("feature" = i,"db_position" = pl.annot[[i]]) }  
    }
  	# Clean empty matches and convert to data frame. Multiple matches are row-concatanated
  	annot.table = do.call(rbind,annot.table)
  	if( !length(annot.table) ) { return (NULL) }
  	annot.table = as.data.frame (annot.table)
  	# Add DB ID
    annot.table[,"ID"] = as.character (DBobj$peaks[annot.table[,2],"ID"])
  	annot.table = annot.table[!is.na(annot.table[,"ID"]),]    
  	##
  	## Build the annotation oject
  	##
		unmatched = data.frame()
		id_list = list()
  	# group all features related to one DB ID
  	for(id in unique(annot.table[,"ID"])) {
  	    print( id )
		    # add to the pseudospectrum the info about the db
    		sameId = annot.table[which(annot.table["ID"]==id),]
        
    #	Apply the CRITICAL PARAM GROUP_TH
		if ( nrow(sameId)<GROUP_TH) {
			 # check if already not in the unmatched list, otherwise add it
		   if(!any(sameId$feature %in% unmatched$feature)) { 
			      whichNotIn = which(!sameId$feature %in% unmatched$feature)
		        unmatched=rbind(unmatched,cbind(sameId[whichNotIn,],pl[sameId$feature[whichNotIn],])) 
		   }
		   next()
		} 
		# add the feature (peak) data and the corresponding library entry
	  id_list[[id]] = list(
				"peaks" = cbind(sameId, pl[sameId$feature,]),
				"db.peaks" = DBobj$peaks[DBobj$peaks[,"ID"]==id,]
    )
		row.names( id_list[[id]]$peaks ) = NULL
		# Max intensity columns for found peaks
    cols = maxIntCols[id_list[[id]]$peaks$feature]
		# decide on a single column according to the hoghest peak
		maxPeak = which.max(sapply(1:length(cols),function(idx) id_list[[id]]$peaks[idx,cols[idx]]))
		intCol = cols[maxPeak]
		## Remove features matching the same library position
 		dbPosition = sameId$db_position
	## Param 1.a: RT integrity of found peaks
	##
	##	Use h.clust to find the number of RT clusters relative to the DB RT cluster	
	d.found = dist (id_list[[id]][["peaks"]]$rt)
	hc.found = hclust (d.found, method="median")      
	d.db = dist (id_list[[id]][["db.peaks"]]$rt.db * 60)
	hc.db = hclust (d.db, method="median")
	cutHeight = max (max(hc.db$height),RT_TH)
	clusters = cutree (hc.found, h=cutHeight)
	id_list[[id]][["rt.clusters"]] = length (unique(clusters))

	valid = FALSE

   for (rtGroup in 1:length(unique(clusters))) {

	found.peaks = id_list[[id]]$peaks[clusters==rtGroup,]
	# filter out low match feature groups
	if ( nrow(found.peaks)<GROUP_TH ) {
	  # check if already not in the unmatched list 
	  next()
	}
  # 	The validity of each RT cluster peak group has to be checked separatly, hence the flag: "valid"
	valid = TRUE       
	annotData = list()
	annotData[["peaks"]] = found.peaks
  	print ("Found a feature group!") 
	##
	## calculate RT shift panalty
	##
	modelShift = coef(rtModel$model)[1]
	## for RT difference I was the unmodified (unshifted) values - since after the correction the true 
	## values can fall on either side of the reference value... so to get a panelty, the peaks have to
	## originally be found "on the other side" of the ref value BEFORE correction:
	rtDiff = median(found.peaks$rt) - median(id_list[[id]][["db.peaks"]]$rt.db*60)
	model.sd = sd(rtModel$model$residuals)
	if (sign(rtDiff) == sign(modelShift)) {
		annotData[["rt.penalty"]] = round(fdrtool::phalfnorm(abs(rtDiff),theta=fdrtool::sd2theta(model.sd)),2)
	} else {
	  # for the contrary trend - double the penalty:
		annotData[["rt.penalty"]] = 2*round(fdrtool::phalfnorm(abs(rtDiff),theta=fdrtool::sd2theta(model.sd)),2)
	}
	##
  ## Param 1.a: calculate a simple coverage score
  ##  
  if (count) {
   	   ## percent of found peaks out of DB entry
   	   annotData[["coverage"]] = round(nrow(found.peaks)/nrow(id_list[[id]][["db.peaks"]]),2)
  }
	##
	## Param 2: find if the main ion is detected, and if coupled to a principal ion
	##
	## Main ion formula (neutral form)
	id_list[[id]][["formula"]] = as.character(DBobj$peaks[(DBobj$peaks$ID==id),"formula"][1])
	##  Main ion mass (neutral ion)
	id_list[[id]][["M.neutral"]] = as.numeric(DBobj$peaks[(DBobj$peaks$ID==id),"M"][1])
  	if (polarity=="positive") { id_list[[id]][["M.parent"]] = id_list[[id]][["M.neutral"]] + protonMass }
  	else if (polarity=="negative") { id_list[[id]][["M.parent"]] = id_list[[id]][["M.neutral"]] - protonMass }
  	# parent ion is matched by mass&rt and it should be confirmed: 
  	# either by a coupled adduct or intensity check.        
  	annotData[["parentIon"]] = 0
  	annotData[["principalIon"]] = 0
  	mainIonRow = rownames(id_list[[id]][["db.peaks"]])[1] # DB main peak 

	if (mainIonRow %in% found.peaks$db_position) {
	   parentIdx = which (found.peaks$db_position == as.numeric(mainIonRow))
	   parentPeak = found.peaks[parentIdx,]
	   isIsotope = grep("M\\+",parentPeak$isotopes)
	   if (!length(isIsotope)) { 
      	     massErr = abs ((id_list[[id]][["M.parent"]] - parentPeak$mz)/id_list[[id]][["M.parent"]]*1e+06)
	     if (nrow(parentPeak)>1) { 
		parentPeak = parentPeak[which.min(massErr),]
		massErr = min(massErr)		
	     }
  if (length(massErr)>1) { massErr = massErr[which.min(massErr)];  }
   ## Use the halfnormal distribution to create a continuous normalized scale
	 ##	ppm = 5	
           annotData[["parentIon"]] = round(1-fdrtool::phalfnorm(massErr, theta=fdrtool::sd2theta(5)),2) 
	  }
  } 

  ## check if a principal ion exists (different from parent ion)
	## get the most abundant peak in the DB entry:
    principalIdx = which.max(id_list[[id]][["db.peaks"]]$lowE)
	## if that peak ain't the parent ion, proceed 
  	if (length(principalIdx) & principalIdx > 1) {
	## find the corresponding peak in the observed (saple) data:
    	dbRow = rownames(id_list[[id]][["db.peaks"]])[principalIdx] 
    	observedIdx = which(found.peaks$db_position==dbRow)      
    	if (length(observedIdx)) {      
      		## remove any multiple hits
      	  principal = id_list[[id]][["db.peaks"]][principalIdx,]
      	  observed = found.peaks[observedIdx,]
	  isIsotope = grep("M\\+",observed$isotopes)
	  if (!length(isIsotope)) { 
	  ## measure the mass error of observed peak to the DB reference peak
	  ## (it is not the theoretical value, but should be a close approximation!)
	    whichAdduct = which(adductRules$name %in% principal$peak.tag)
	    if (length(whichAdduct)) {    	
         	adductRule = adductRules[whichAdduct,]
         	adductForm = parseCameraAdducts (adductRule,id_list[[id]][["formula"]])
         	adductMass = Rdisop::getMolecule(adductForm)$exactmass
         	massErr = abs((adductMass - observed$mz)/adductMass*1e+06)	
       	    } else {
        	massErr = abs((principal$mz - observed$mz)/principal$mz*1e+06)
      	    }
    # eliminate multiple scoring hits by min mass error
	    if (length(massErr)>1) { 
	      observed = observed[which.min(massErr),]
	      massErr = massErr[which.min(massErr)] 
	      }
	    ## use the same error function as for parent peak, with est. SD.massErr = 6ppm
      	    annotData[["principalIon"]] = round(1-fdrtool::phalfnorm(massErr, theta=fdrtool::sd2theta(6)),2) 		
    	   }
  	 }
  	} 
    
#    runMatch("__tests__/JAsig_neg_2ch_PL_new.tsv","__tests__/stdMix","negative","JAsig_test",DBobj)  
  ## Adjust parent ion annotation according to principal ion annotation:
  ##	if principal ion is NOT detected apply intensity filter	on parent peak annotation
  if ( !(annotData[["principalIon"]]) & annotData[["parentIon"]] ) {
	# check intensity ratio between MS channels - but only when the principal ion was not detected
	  if ( mean(as.numeric(parentPeak[colsMs1])) < mean(as.numeric(parentPeak[colsMs2])) ) {
		  annotData[["parentIon"]] = 0	
	  }
  # now, check if principal ion is also valid 
  } else if ((annotData[["principalIon"]])) {
	# in cases the principal ion is a higher mass adduct - do MS channels comparison to lowerreduce FPs
    # in the other cases - where the mass is lower then mol. mass - we skip the intensity comparison 
    # to allow for cases where a true fragment (e.g. loss of H2O) is also the principal ion - and at
    # the same time it's also observed as stronger ion on Ms2 channel due to increased fragmentation
	if  ((observed$mz > id_list[[id]]$M.parent) & 
	     (mean(as.numeric(observed[colsMs1])) < mean(as.numeric(observed[colsMs2]))) ) 
	{
		annotData[["principalIon"]] = 0
	}
  }
  ##
  ## Param 3: find how many peaks match a tagged DB peak (apart from the main/principal ions)
  ##
  annotData[["tagged"]] = 0
	dbPeaks = id_list[[id]][["db.peaks"]][rownames(id_list[[id]][["db.peaks"]]) %in% found.peaks$db_position,]
	dbTags = which (dbPeaks$peak.tag %in% adductRules$name[-1])
	if (length(dbTags)) {
	          adductRule = adductRules[adductRules$name %in% dbPeaks$peak.tag[dbTags],]
      	    for (i in 1:nrow(adductRule)) { 
	      adductForm = parseCameraAdducts (adductRule[i,],id_list[[id]][["formula"]])	  
              tryCatch ({
                adductMass = Rdisop::getMolecule(adductForm)$exactmass
                massErr = min(abs(adductMass - found.peaks$mz))/adductMass * 1e+6
                annotData[["tagged"]] = annotData[["tagged"]] + round(1-fdrtool::phalfnorm(massErr,theta=fdrtool::sd2theta(6)),2)    
              }, error = function(err) {
        	    print (paste("Warning:",err))
              })    
      	    }	          
	}
  
  ##
  ## Param 4: apply the decompose isotopes module				
  ##    
	if (isoMatch) {
		# chkIsotopes function
	 annotData[["isotopes"]] = getIsotopes (found.peaks,intCol)
	 if (!is.null(unlist(annotData[["isotopes"]])))  {
	   annotData[["iso.confirmed"]] = sum(unlist (run_isoMatchModule (
					annotData[["isotopes"]], 
					intCol, 
					id_list[[id]][["M.neutral"]],
					id_list[[id]][["formula"]],
					adductRules,
					isoMatchModel)))
     	   if (is.null(annotData[["iso.confirmed"]])) { annotData[["iso.confirmed"]] = 0 }
	} else {
		annotData[["isotopes"]] = NA
		annotData[["iso.confirmed"]] = 0
	}
}       

  ##
  ## Param 5: perform highE fragment matching, using my implementation of X-rank algorithm
  ##
	if (xrank) {
  ## first check if the library has ranked fragments:
	db.ranks = DBobj$ranked[[which(names(DBobj$ranked)==id)]]
	  if (!is.null(db.ranks)) {
	  	frgMatch = which (found.peaks[,"mz"] < id_list[[id]][["M.neutral"]]-8*protonMass)
      		## define the fragments in the observed sample
		if (length(frgMatch)) {
		  ## rank fragments by peak intensity
		  frgs = found.peaks[frgMatch,]
		  ## remove peaks labelled as isotopes 
		  ## I do not trust the camera isotopes annotations - use isotopes found by the internal function
      if (!any(is.na(annotData[["isotopes"]]))) {
      	whichIso = unlist(sapply (annotData[["isotopes"]], function(x) 
				as.numeric(x[2:nrow(x),"feature"])))
      	 if (sum(frgs$feature %in% whichIso)) {
			      frgs = frgs[-which(frgs$feature %in% whichIso),] 
      	 }
      }
		  ## if any peaks left after isotope removal - do the matching:
		  if (nrow(frgs)) {        
		    	db.ranks = DBobj$ranked[[which(names(DBobj$ranked)==id)]]
		    	db.ratio = id_list[[id]][["db.peaks"]][1,"lowE"] > db.ranks[1,"lowE"]
			annotData[["ranked.frgs"]] = round(
				frgs[order(frgs[,intCol+1],decreasing=T),"mz"][1:min(nrow(frgs),MAX_RANK)],4) 
		  ## filter (optional, on DB side) peaks with less than INT_TH relative intensity			
			if (INT_TH) {
			    whichTooLow = which(db.ranks$highE < max(db.ranks$highE)*INT_TH/100)
			    if (length(whichTooLow)) {
			    	db.ranks = db.ranks[-whichTooLow,]
			    }			 			    
			}     
			id_list[[id]][["ranked.db.frgs"]] = db.ranks[,1] 			  	    
			sampleScore = scoreSpectra (
				    annotData[["ranked.frgs"]],
				    id_list[[id]][["ranked.db.frgs"]],
				    xrankScoremat )
        if (sampleScore>=0) {
          # get the max positive score
          fullScore = scoreSpectra (
            id_list[[id]][["ranked.db.frgs"]],
            id_list[[id]][["ranked.db.frgs"]],
            xrankScoremat )
        } else {
          # compute the max negative score
          fullScore = scoreSpectra (
            rep(0,length(id_list[[id]][["ranked.db.frgs"]])),
            id_list[[id]][["ranked.db.frgs"]],
            xrankScoremat ) * (-1)         
        }
    ## Using the "native" X-rank score without normalization - this will give priority to higher 
    ##	fragment hits, which is the original idea of Roman Milonas...
			annotData[["xrank"]] = round(sampleScore,1)
    ## Alternatively using the normalized score - for compatibility with otehr scores:
    #		      annotData[["xrank"]] = round(sampleScore/fullScore,2)              
			 } 
		} else { 
      # In case of no fragments for library match with fragments - apply a penalty:
  		db.ranks = db.ranks[,1]
  		mismatch.score = scoreSpectra( rep(0,length(db.ranks)), db.ranks, xrankScoremat )
  		annotData[["xrank"]] = round(mismatch.score,1)
		}		  
	 } else {
		 annotData[["xrank"]] = 0
	 }
	} 
  ##	build annotation entry (for each rt peak group):
	  id_list[[id]][[paste("annot",rtGroup,sep=".")]] = annotData
	} # end rtGroups loop
	## remove non-valid annotations (nPeaks<GRP_TH)
	if (!valid) {  
	  id_list[id] = NULL
	} 
	} # end "for(id in unique(annot.table[,"ID"]))"
	output = list( 
		"groups" = id_list,
		"singles" = unmatched )
	return ( output )
}

##################################################################
##################################################################

xrefIonModes = function(foundPeaks, ref.list, id, xref.table) {
  # calculate compatible masses between ionizaton modes
	# find ID match
  findID = which(names(ref.list$groups) %in% id)
	xref.score = 0
	if (length(findID)) {
    
    filterOther = sapply (grep("annot",names(ref.list$groups[[findID]])),function(i) {
      X=ref.list$groups[[findID]][[i]];
      if(as.numeric(X["parentIon"])==0 & as.numeric(X["principalIon"])==0) { return(NULL) }
      else { return(1) }      
      })
    if (all(is.null(unlist(filterOther)))) { return(0) }
	   ## compare group RTs - difference should be very similar (here it's 10 seconds)	
	   refPeaks = ref.list$groups[[findID]]$peaks
	   if (abs(median(refPeaks$rtShifted)-median(foundPeaks$rtShifted)) < 10) {
	     ## find and sum matching pairs scores
	     xref.score = sum (apply (foundPeaks,1,function(found) {
		    mDiff.df = data.frame (
		      massDiff = as.numeric(found["mz"]) - refPeaks$mz,
		      massTol = as.numeric(found["massTol"]) + refPeaks$massTol )
		    # check ref table for corresponding mass diff match within massTol bounds
		    max(as.numeric(unlist(apply(mDiff.df,1,function(x) {
		      # here, do the mass cross matching 
		      xref.match = which (abs(xref.table-x["massDiff"]) < x["massTol"])
		      if (length(xref.match)) {
			      # find the minimal mass diff
			      bestMatch = which.min(abs(xref.table[xref.match]-x["massDiff"]))
			      # score match using mass error 
			      massErr = abs(xref.table[xref.match[bestMatch]]-x["massDiff"])/as.numeric(found["mz"])*1e+6
			      # There is a logical bug: PI mode gets higher scores because of lower PPM values for mass errors!
            # need to convert scoring to absolute units (Dalton)
            xref.score = round(1-fdrtool::phalfnorm(as.numeric(massErr),theta=fdrtool::sd2theta(massSD)),2)
		      } else {
		        xref.score = 0
		      }}   
        ))))
	     }))
	   }
	} else {
	   ## Optional: look at "single" peaks discarded earlier - more risky
	   ##
	   ## rtMatch = which(ref.list$singles$rt-median(foundPeaks$rt)<RT_TH)
	   ## mzMatch = foundPeaks$mz - ref.list$singles$mz[rtMatch]
	   ## browser()	
	}
  return(xref.score)
}

## Deprecated: use msLib 'convertResults' instead...
summarizeResults = function(mainList,DBpeaks,mode,refList=NULL,xref.table=NULL) {
  # Convert the IdList 'getIds' results into text table
  # if both ion modes available - cross check both modes
  # output present results in nicely readable table format
  data (glmModel)
  # Calc total intensity values per sample/channel
  id.list = mainList$groups
  # Again -reading correct intensity columns depends on (01) and (02) postfix for MS channels
	intCols1 = grep(".*01$", names(id.list[[1]]$peaks))
	intCols2 = grep(".*02$", names(id.list[[1]]$peaks))
  param.cols = c("coverage","parentIon","principalIon","adducts","isoConfirmed","xrank","xref","rtPenalty")
  results = data.frame(row.names = names(id.list),stringsAsFactors = FALSE)
  
  for (i in 1:length(id.list)) {																																												
	## loop over all possible annotation groups ("annot.1","annot.2",etc.) 
	  annotations = grep ("annot",names(id.list[[i]]))
	  for (j in 1:length(annotations)) {
     # print(paste(i,j))
	    annot.data = id.list[[i]][[annotations[j]]]      
      # Molecular ion filter - skip annotations with no parent or principal ion
	    if (!annot.data$parentIon & !annot.data$principalIon) {	next() }
      # Debug: this should not happen:
    	if (length(annot.data$parentIon)>1 | length(annot.data$principalIon)>1) { browser() }
      ID = as.character(names(id.list[i]))
      if (!is.null(xref.table)) {
        # validate annotation with the other ion mode
	      if (mode=="negative") {
	        annot.data$xref = xrefIonModes(annot.data$peaks,refList,ID,xref.table)	
    
	      } else if (mode=="positive") {
	        # invert xref table for PI mode:
	        annot.data$xref = xrefIonModes(annot.data$peaks,refList,ID,t(xref.table)*(-1))
	      }} else {
        annot.data$xref = 0
      }
    # Data of library entry
    DBdata = DBpeaks[which(DBpeaks$ID == ID)[1],] 	
    DBdata$peak.table = basename(as.character(DBdata$peak.table))
	  # fill output list
  	intCols = grep(".*0[1|2]$", names(annot.data$peaks))
	  meanInt = round(colMeans(annot.data$peaks[,intCols1]))
    results = rbind(results, cbind(
	      ID = ID,
	      polarity = mode,	
	      rt.cluster = j,
        nPeaks = nrow(annot.data$peaks),
	      coverage = annot.data$coverage,
	      parentIon = annot.data$parentIon,
	      principalIon = annot.data$principalIon,
	      adducts = annot.data$tagged,
	      isoConfirmed = annot.data$iso.confirmed,
	      xrank = annot.data$xrank,        
	      xref = annot.data$xref, 
	      rtPenalty = annot.data$rt.penalty,
          DBdata = "###",
    	      DBdata,
	      intData = "@@@",
	      "found.rt" = median(round(annot.data$peaks$rt/60,2)),        
	      "found.mz" = paste(round(annot.data$peaks$mz,4),collapse=";"),
		    "main.mass" = annot.data$peaks$mz[which.max(rowMeans(annot.data$peaks[,intCols1]))],
		    "max.samp" = names(which.max(colMeans(annot.data$peaks[,intCols1]))),
		    t(meanInt),
		    "pcgroup1" = names(which.max(table(annot.data$peaks$pcgroup)))        
		  ))
	  }
  }
  rownames(results) = NULL
  # evaluate hits using a predefined glm model
  results$Scoring = "@@@"
  results$score = predict(glmModel,newdata=results[,5:12],type="response")
  results = results[order(results$score,decreasing=TRUE),]
  return (as.data.frame(results))
}

