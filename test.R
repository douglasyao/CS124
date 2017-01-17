### CODE FOR TESTING ###
source(generate_read_matrix.R)
source(assemble_haplotype.R)

setwd('~/Downloads/Final_Project_Data/Final_2/medium/test/input/')
data <- read.table('reads_high_error.txt', stringsAsFactors = F)
hap <- get_haplotype(data)
write.table(hap, file = 'medium_difficulty_test_haplotypes_high_error_Douglas_Yao_304478204.txt', sep = '\t', quote = F, row.names = F, col.names = F)

reference <- scan('../output/haplotypes_low_error.txt', what = 'character')
if (identical(reference[1],hap[1]) || identical(reference[1],hap[2])) {
  print ('No errors')
} else {print ("Errors present")}

# check accuracy of haplotype assembly after varying error rate
error_readings <- NULL
reads <- get_reads(10,100,4:6, 0.3)
testhap <- get_haplotype(reads[[3]])
a <- strsplit(reads[[1]],'')[[1]]
b <- strsplit(testhap[1],'')[[1]]
length(a[a!=b])/10




