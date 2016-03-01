#################################################################
# X-rank module
#################################################################s

scoreSpectra = function( ranked1, ranked2, scoreMatrix, massTol=50 ) {
	# return the symmetric match score of input spectra 
	# according to the X-rank method by Mylonas et. al.

  print ("############") 
  print ("Xrank Module")
  print ("############")
  
	# If ranks match, record the rank of the match

	len1 = length(ranked1[!is.na(ranked1)])
	len2 = length(ranked2[!is.na(ranked2)])
  
	print ("Xrank scoring module")
	matches = rep(NA,len1)
	for( rank in 1:len1 ) {
		match = which(abs(ranked2-ranked1[rank])/ranked1[rank]*1e06 < massTol)
		if (length(match)>1) {  
		   match = match[which.min(abs(ranked2[match]-ranked1[rank]))]
		}
#		print( match )
		if (length(match)) { matches[rank] = match }
	}
	# Convert matche and mismatches to Xrank scores
	scores = sapply( 1:len1, function(rank) {
			if (!is.na(matches[rank])) {
			# score a match by its expected (column) and observed (row) rank
				scoreMatrix[matches[rank],rank]
			} else {
			# score a mismatch
				scoreMatrix[31,rank]
			}
		})
	# Compute the inverse match of target spectra to query
	matches = rep(NA,len2)
	for( rank in 1:len2 ) {
		  match = which(abs(ranked1-ranked2[rank])/ranked2[rank]*1e06 < massTol)
	  	if (length(match)>1) {
		     match = match[which.min(abs(ranked2[match]-ranked1[rank]))]
		  }
		  if (length(match)) { matches[rank] = match }
	}
	# Convert to Xrank scores
	invScores = sapply( 1:len2, function(rank) {
			if (!is.na(matches[rank])) {
				score = scoreMatrix[matches[rank],rank]	# score a match
			} else {
				score = scoreMatrix[31,rank]		# score a mismatch
			}
		})
	
	print( mean(c(sum(scores),sum(invScores))) )
# return the overall symmetric X-rank similarity score
	return( mean(c(sum(scores),sum(invScores))) )
}

