######################################
# the reference code of program2 
######################################
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library(Biostrings)
######################################
# initial
######################################
# read parameters
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: Rscript pro2_<your student ID>.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}

# parse parameters
i<-1 
while(i < length(args))
{
  if(args[i] == "--input"){
    i_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--score"){
    s_f<-args[i+1]
    i<-i+1
  }else if(args[i] == "--aln"){
    aln_mode <- args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_open"){
    g_o<-args[i+1]
    i<-i+1
  }else if(args[i] == "--gap_extend"){
    g_e<-args[i+1]
    i<-i+1    
  }else if(args[i] == "--output"){
    o_f<-args[i+1]
    i<-i+1
  }else{
    stop(paste("Unknown flag", args[i]), call.=FALSE)
  }
  i<-i+1
}

#print("PARAMETERS")
#print(paste("input file         :", i_f))
#print(paste("output file        :", o_f))
#print(paste("score file         :", s_f))
#print(paste("gap open penalty   :", g_o))
#print(paste("gap extend penalty :", g_e))
######################################
# main
######################################
# your code
file_fasta <- readBStringSet(i_f)
seq_name = names(file_fasta)
sequence = paste(file_fasta)
score_schem <- read.table(s_f, header = T, skip = 9)
score_schem_matrix <- t(matrix(unlist(score_schem), ncol = 24, byrow = TRUE))
rownames(score_schem_matrix) <- rownames(score_schem)
colnames(score_schem_matrix) <- rownames(score_schem)
g_o_int <- strtoi(g_o)
g_e_int <- strtoi(g_e)
pA <- pairwiseAlignment(sequence[1], sequence[2], substitutionMatrix = score_schem_matrix, gapOpening = g_o_int, gapExtension = g_e_int, type = aln_mode)
r = BStringSet( c(toString(pattern(pA)), toString(subject(pA))) )
names(r) <- seq_name
writeXStringSet(r, o_f, format="fasta")
print("generate finish")
