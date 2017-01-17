### GENERATE A RANDOM READ MATRIX ####
# generates a random read matrix in the exact format of the examples provided by Michael
# user can specify the length of the haplotype, the number of reads, the range of the lengths of reads, and the error rate
# return type is a list containing the two haplotypes and the read matrix
get_reads <- function(hap_length, num_reads, read_length_range, error_rate) {
  require(stringr)
  # because reads are to be "cut off" at the end (as in the example read matrices), extend the haplotype length for now
  hap_length <- (hap_length + max(read_length_range))
  
  # generate random binary string for haplotype
  hap1 <- sample(c(0,1),hap_length, replace = T)
  names(hap1) <- 1:hap_length
  
  # generate the complementary hapotype
  hap2 <- 1-hap1
  
  # randomly determine the locations of the reads
  starting_locations <- sample(c(1:(hap_length-max(read_length_range))),num_reads, replace = T)
  
  # randomly determine the lengths of the reads
  lengths <- sample(read_length_range,num_reads,replace = T)
  
  # randomly determine which haplotype each read came from
  hap1_or_hap2 <- sample(c(1,2),num_reads,replace=T)
  
  # create a reference table which has the location, length, and which haplotype for each read
  table <- as.data.frame(cbind(starting_locations,lengths,hap1_or_hap2))
  table <- table[order(table$starting_locations),]
  rownames(table) <- 1:num_reads
  reads = NULL
  
  # generate the reads
  for (i in 1:num_reads) {
    
    # if read comes from the first haplotype
    if (table[i,3] == 1) {
      
      # subset the haplotype to create the read of proper length
      read <- hap1[table[i,1]:(table[i,1] + table[i,2] - 1)]
      read <- paste(read, collapse = '')
      
      # add "-"s to the beginning and ends of the read to indicate the position of the read
      beginning <- paste(rep('-',table[i,1]-1), collapse = '')
      ending <- paste(rep('-',(hap_length-table[i,2]-table[i,1]+1)), collapse = '')
      reads <- c(reads, paste0(beginning,read,ending))
    } 
    
    # if the read comes from the second haplotype
    else {
      
      # subset the haplotype to create the read of proper length
      read <- hap2[table[i,1]:(table[i,1] + table[i,2] - 1)]
      read <- paste(read, collapse = '')
      
      # add "-"s to the beginning and ends of the read to indicate the position of the read
      beginning <- paste(rep('-',table[i,1]-1), collapse = '')
      ending <- paste(rep('-',(hap_length-table[i,2]-table[i,1]+1)), collapse = '')
      reads <- c(reads, paste0(beginning,read,ending))
    }
  }
  
  reads <- str_split_fixed(reads,'',n = nchar(reads[1]))
  
  # subset the read matrix so it is the proper length. Remember that we extended the haplotype length at the beginning so the final reads would be "cut off".
  reads <- reads[,1:(ncol(reads)-max(read_length_range))]
  count = 0
  
  # insert errors into the reads by switching random 0s with 1s and vice versa. 
  while (count < num_reads*(sum(read_length_range)/length(read_length_range))*error_rate) {
    error <- sample(1:(num_reads*ncol(reads)), 1)
    row <- ceiling(error/hap_length)
    column <- error-(ceiling(error/ncol(reads))-1)*ncol(reads)
    if (reads[row,column] == 1) {
      reads[row,column] = 0
      count = count + 1
    }
    else if (reads[row,column] == 0) {
      reads[row,column] = 1
      count = count + 1
    }
  }
  
  # subset the haplotypes to the proper length
  hap1 <- hap1[1:(hap_length-max(read_length_range))]
  hap2 <- hap1[1:(hap_length-max(read_length_range))]
  hap1 <- paste(hap1, collapse = '')
  hap2 <- paste(hap2, collapse = '')
  return (list(hap1,hap2,reads))
}
