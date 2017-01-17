### GENERATES A VECTOR OF PARTIAL HAPLOTYPES AT EACH HAPLOTYPE POSITION FROM A READ MATRIX ###
# input data is vector of strings containing 0s, 1s, and -'s representing each row of the matrix. In other words, this is the format of the input matrices in Michael's examples. 
# partial haplotypes for columns with less than a certain amount of reads are represented with NA
get_haplotype <- function(data) {
  require(stringr)
  
  # splits up the reads (which are strings) into a character matrix consisting of 0s, 1s, and -'s. Note that this requires the stringr package.
  data <- str_split_fixed(data[,1], '', n = nchar(data[1,]))
  mec = NULL
  
  # find optimal haplotype for longest read starting at each column in matrix. Last column is excluded because reads of length 1 are trivial.
  for (x in 1:(ncol(data)-1)) {
    
    group <- data[!(data[,1] == '-'),] #obtain all reads that begin at indicated position
    
    # matrix with one row/column automatically converted to vector, so convert back
    if (is.vector(group)) {
      group <- matrix(group,nrow = 1)
    }
   
     # if no reads at particular column, then add NA
    if (length(group) == 0) {
      mec <- c(mec, NA)
      data <- data[,-1]
    }
    
    else {
      data <- data[-(1:nrow(group)),-1] #crop read matrix so next column becomes first column
      
      # if -'s present, remove dashes and paste characters together to individual strings
      if ('-' %in% group) {
        group <- do.call(paste0, as.data.frame(group)) 
        group <- sapply(group, function (x) substr(x,0,regexpr('-',x)[1]-1))
      }
      else {
        group <- do.call(paste0, as.data.frame(group))
      }
      names(group) <- NULL
      num <- max(nchar(group)) #length of longest read
  
      # permute through all combinations of 0s and 1s for longest read and determine which permutation results in the least amount of total errors among all the reads.
      # if less than 2 reads present at particular position, then skip 
      # creates a vector 'mec' that contains optimal haplotype for each column position in read matrix
      if (num > 0 && length(group) >= 2) {
      
        perms <- expand.grid(rep(list(0:1),num)) #obtain all permutations of 0s and 1s
        perms <- perms[1:(2^num/2),] #throw out all permutations that are complements of each other
        perms <- do.call(paste0, as.data.frame(perms)) #paste to string
        errors = NULL
      
        # compare each permutation with reads
        for (i in 1:length(perms)) {
          totalerror = 0;
        
          # for each permutation, compare with each read
          for (j in 1:length(group)) {
            perm1 <- substr(perms[i],0,nchar(group[j])) # subset permutation to match length of read
            perm2 <- paste(1-as.numeric(strsplit(perm1,'')[[1]]), collapse = '') #obtain complement of permutation
            error1 <- sum(mapply(function(x,y) sum(x!=y),strsplit(group[j],''),strsplit(perm1,''))) #compare read to permutation and record the error
            error2 <- sum(mapply(function(x,y) sum(x!=y),strsplit(group[j],''),strsplit(perm2,''))) #compare read to complement of permutation
            totalerror = totalerror + min(error1,error2) #between permutation and its complement, add the lower error to the total error count for the particular permutation
          }
          errors = c(errors, totalerror) #add total error for particular permutation to vector containing errors for all permutations
        }
        names(errors) <- perms
        mec <- c(mec, names(which.min(errors))) #find which permutation has the least amount of errors
      }
    
    # if fewer than 2 reads present, add a NA to the vector
      else {
        mec <- c(mec, NA) 
      }
    }
  }
  hap <- assemble_fragments(mec) #assembles haplotype. see below 
  hap <- c(hap, paste(1-as.numeric(strsplit(hap,'')[[1]]), collapse = '')) #obtain complement of haplotype
  return (hap)
}



### ASSEMBLES THE FULL HAPLOTYPE OUT OF THE PARTIAL HAPLOTYPES
# This function assembles a haplotype based on partial haplotypes. Also does error correction in case partial haplotypes do not match up and accounts for NAs. 
# Input is a vector of strings containing one partial haplotype for each column position in the read matrix (excluding the last column)
# If a particular partial haplotype is missing from a position, then it must be indicated with a NA
assemble_fragments <- function(mec) {
  
  # if first partial haplotype is NA, then add a '-' to the beginning of the haplotype. 
  # This only applies to the very first partial haplotype because if later haplotypes are NA, they can simply be skipped over without affecting the overall haplotype
  if (is.na(mec[1])) {
    hap = NA
  }
  
  # initialize haplotype with first partial haplotype
  else {
    hap = mec[1]
  }
  
  # assemble all remaining snippets
  for (i in 2:(length(mec))) {
    
    # if first partial haplotype is missing, then add a '-' to the beginning 
    if (is.na(hap)) {
      hap <- paste0('-',mec[i])
    }
    
    # partial haplotypes that are NA are simply skipped over
    else if (!is.na(mec[i])) {
      compare <- substring(hap,i) #subset haplotype to start with position matching the position of the partial haplotype  
      add <- mec[i] #partial haplotype
      add2 <- paste(1-as.numeric(strsplit(add,'')[[1]]), collapse = '') #complement of partial haplotype
      
      # if snippet is shorter or equal length to the previous snippet, it will not contribute additional information to the haplotype and can be discarded
      if (nchar(compare) < nchar(add)) {
        compare2 <- substr(add,0,nchar(compare)) #subset partial haplotype to match overlapping section with overall haplotype
        compare3 <- substr(add2,0,nchar(compare)) #same as above, but with complement
        
        # if partial haplotype matches haplotype, then append extra 0s or 1s from partial haplotype to haplotype
        if (compare == compare2) {
          hap = paste0(hap,substring(add,nchar(compare)+1))
        }
        
        # same as above but with complement of partial haplotype
        else if (compare == compare3) {
          hap = paste0(hap,substring(add2,nchar(compare)+1))
        }
        
        # If partial haplotype does not match with overall haplotype, then perform error correction. Mismatches occur in approximately 1% of all partial haplotypes.
        # There are several possible scenarios for errors to occur. If there is an error in the very first position of the whole haplotype, then it is unfortunately impossible to detect.
        # If an errors occurs anywhere among the first to second to last position in the partial haplotype, then the partial haplotype can be discarded without issue.
        # However, because snippets are assembled one by one, if an error occurs in the last nonoverlapping position(s) that doesn't overlap with the haplotype, then the snippet will still be appended.
        # The error will be only be detected when subsequent snippets don't match up with this location. The code will account for this case.
        else {
          
          # assemble haplotype for subsequent 5 partial haplotypes. This will reveal information about the nature of the error in the errant snippet. 
          result <- assemble_fragments(mec[(i+1):(i+6)])
          result <- substr(result, 0, nchar(mec[i])-1) #subset this temporary haplotype to match overlapping sections with errant partial haplotype
          
          # if the errant snippet is followed immediately by an NA (extremely rare, but happens)
          if (substr(result,0,1) == '-') {
            result <- substring(result,2) #remove '-' from beginning of temporary haplotype 
            result2 <- paste(1-as.numeric(strsplit(result,'')[[1]]), collapse = '') # rest of code is the same as below
            comp <- substring(mec[i],3)
            if (comp == result || comp == result2) {
              if (i == 2) {
                hap <- paste0(substring(hap,0,1), mec[i])
              }
              else {
                difference <- nchar(mec[i-1])-nchar(mec[i-2]) + 1
                hap <- substring(hap, 0, nchar(hap)-difference) 
              }
            }
          }
          
          # if the errant snippet is followed by a normal snippet
          else {
            result2 <- paste(1-as.numeric(strsplit(result,'')[[1]]), collapse = '') #find complement of temporary haplotype
            comp <- substring(mec[i],2) #subset errant snippet to match overlapping sections with temporary haplotype 
            
            # If the errant snippet matches the snippets following it, then that means the final 0(s) or 1(s) of the preceding snippet was where the error was located.
            # If the errant snippet does not match the snippets following it, then the error was located in the snippet. In this case, the snippet is simply discarded.
            if (comp == result || comp == result2) {
              
              # In case the error was located in the very first partial haplotype, then the error is not guaranteed to be the in the last position(s).
              # In this case, simply discard the entire first snippet except for the first position.
              if (i == 2) {
                compare1 <- substring(hap,2)
                compare2 <- paste(1-as.numeric(strsplit(compare1,'')[[1]]), collapse = '')
                compare3 <- substring(mec[i],0,nchar(hap)-1)
                error1 <- sum(mapply(function(x,y) sum(x!=y),strsplit(compare1,''),strsplit(compare3,''))) #compare read to permutation and record the error
                error2 <- sum(mapply(function(x,y) sum(x!=y),strsplit(compare2,''),strsplit(compare3,''))) #compare read to complement of permutation
                if (error2 < error1) {
                  hap <- paste(1-as.numeric(strsplit(hap,'')[[1]]), collapse = '')
                }
                hap <- paste0(substring(hap,0,1), mec[i])
              }
              
              # Remove all final nonoverlapping positions from previous snippet
              else {
                position = i-1
                previous = mec[position]
                while (is.na(previous) & position > 0) {
                  position = position - 1 
                  previous = mec[position]
                }
                
                if (position <= 1) {
                  hap <- substring(hap, 0, nchar(hap)-1)
                }
                
                else {
                  previous_position = position - 1
                  previous_previous = mec[previous_position]
                 
                  while (is.na(previous_previous) && previous_position > 0) {
                    previous_position = previous_position - 1 
                    previous_previous = mec[previous_position]
                  }
                  
                  if (previous_position == 0) {
                    compare1 <- substring(previous,2)
                    compare2 <- paste(1-as.numeric(strsplit(compare1,'')[[1]]), collapse = '')
                    compare3 <- substring(mec[i],0,nchar(previous)-1)
                    error1 <- sum(mapply(function(x,y) sum(x!=y),strsplit(compare1,''),strsplit(compare3,''))) #compare read to permutation and record the error
                    error2 <- sum(mapply(function(x,y) sum(x!=y),strsplit(compare2,''),strsplit(compare3,''))) #compare read to complement of permutation
                    if (error1 < error2) {
                      hap <- paste0(substring(hap,0,1), substring(previous,0,1),compare3)
                    }
                    else {
                      hap <- paste0(substring(hap,0,1), substring(paste(1-as.numeric(strsplit(previous,'')[[1]]), collapse = ''),0,1),compare3)
                    }
                  }
                  
                  else {
                    difference <- nchar(previous)-nchar(previous_previous) + 1
                    hap <- substring(hap, 0, nchar(hap)-difference) 
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return (hap)
}
