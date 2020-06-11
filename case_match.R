
#Martin  -- we'll try and get this function SelectMatchedControls written quickly.  Luke, could you have a go at this?
#		Can I check with you the exact specification for this function:
#		
#		The function will be called once for each SpecimenDate value
#		
#		It will take three arguments:
#		
#		1. cases - a data table with two integer columns upi, mkey.  upi values are unique, mkey values may be duplicated.  Typically this table will have up to 1000 rows.
#		
#		2. popextract - a data table with the same two integer columns.  upi values are unique, mkey values may be duplicated.  Typically this table will have about 40,000 records
#		
#		3. first.stratum.number - a positive integer
#		
#		For each mkey value in cases, count the number n[k] of cases with this mkey value, and sample (with replacement) up to
#       10 * n[k] records (if < 10 * n[k], as many as are available) with this mkey value from popextract
#		
#		Append cases and controls to a new table, with a third column named stratum with integer values incremented for each mkey value. starting at first.stratum.number
#		
#		For each value of SpecimenDate, Martin should call this function with the following arguments
#		
#		1. cases =  subset of all cases with that SpecimenDate
#		
#		2. popextract = all popextract alive on that SpecimenDate with an mkey value matching an mkey value in cases. 
#		
#		3.  first.stratum.number: next unused stratum number, starting from 1. 
#		
#		Martin -- if you fix your code to work with this function, we will get it written.
#' Given data frames containing cases and a background population, randomly sample up to 10 records from
#'   the background population for each case, with strata matching on the mkey column. A data.frame containing both cases and matched
#'   controls are returned, along with a flag identifying which are cases/controls.
#' @param cases data.frame - columns: upi, mkey: upi unique character, mkey integer duplicated - up to 1000 rows
#' @param popextract data.frame - same columns as cases
#' @param first.stratum.number +ve int
#' @return A data.frame of mixed cases and controls. Columns upi (character),
#'   mkey (integer), strata (integer) and is.case (boolean). This should be no more than 11*size(cases) in size.
SelectMatchedControls <- function(cases, popextract, first.stratum.number) {

	cases <- data.table::data.table(cases)
	popextract <- data.table::data.table(popextract)
	
    # Make a copy of cases so do not alter original by reference
	cases.mkey <- data.table::copy(cases)
	
	# For each mkey value in cases, count number n[k] of cases with this mkey value
	cases.mkey[, num.cases := .N, by=mkey]
	# For each mkey in population how many cases are there
	pop.m <- merge(popextract, unique(cases.mkey[, .(mkey, num.cases)]), how="left", on="mkey")
	# For each mkey in population how many potential controls are there
	pop.m[, num.pop.mkey := .N, by=mkey]
	
	# All cases are returned - then controls are appended
	ret.table <- cases
    ret.table$is.case <- TRUE
	
	# If <10*n[k] matched in population then take all available (no replacement)
    ret.pop.lt.10 <- pop.m[(num.cases > 0) & (num.pop.mkey < 10 * num.cases), .(upi, mkey)]
    ret.pop.lt.10$is.case <- FALSE
	ret.table <- rbind(ret.table, ret.pop.lt.10)
	
	
	
	# For >=10*n[k] in each mkey in popextract then randomly sample 10 * n[k] without replacement
	pop.m.gt.10 <- pop.m[(num.cases > 0) & (num.pop.mkey >= 10 * num.cases)]
	
	# Randomly order the data.table for sampling
	rand.order <- sample(nrow(pop.m.gt.10))
	ret.pop.gt.10 <- pop.m.gt.10[rand.order]
	
	# Select first 10*num.cases from each subgroup (.SD)
	ret.pop.gt.10 <- ret.pop.gt.10[, .SD[1:(10*unique(.SD$num.cases))], by=mkey]
	
    ret.pop.gt.10$is.case <- FALSE
	ret.table <- rbind(ret.table, ret.pop.gt.10)
	
	# Add stratum - incrementing integer for each mkey starting at first.stratum.number
	mkey.vals <- unique(ret.table$mkey)
	stratum.vals <- seq(first.stratum.number,first.stratum.number+length(mkey.vals) - 1)
	ret.table[, stratum := as.integer(as.character(factor(mkey, mkey.vals, stratum.vals)))][]
	
	return(as.data.frame(ret.table))
}
	
	

