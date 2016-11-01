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
  stop("USAGE: Rscript hw2_105753006.R --input test.fasta --score pam250.txt --aln global --gap_open -10 --gap_extend -2 --output result.fasta", call.=FALSE)
}
aln_mode <- "global"
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
#read fasta, and create var to save id,seq
file_fasta <- readBStringSet(i_f)
seq_name <- names(file_fasta)
sequence <- paste(file_fasta)

#read score file, convert to matrix and give it row and col names 
score_schem <- read.table(s_f, header = T, skip = 9)
score_schem_matrix <- t(matrix(unlist(score_schem), ncol = 24, byrow = TRUE))
rownames(score_schem_matrix) <- rownames(score_schem)
colnames(score_schem_matrix) <- rownames(score_schem)

#convert string g_o and g_e to int
g_o_int <- strtoi(g_o)
g_e_int <- strtoi(g_e)

#convert seq to list
s1_list <- strsplit(sequence[1], "")[[1]]
s2_list <- strsplit(sequence[2], "")[[1]]

#create matrix to compute score with gap-penalty
align_matrix <- matrix(-Inf, nrow = length(s1_list) + 1, ncol = length(s2_list) + 1)
ix <- matrix(-Inf, nrow = length(s1_list) + 1, ncol = length(s2_list) + 1)
iy <- matrix(-Inf, nrow = length(s1_list) + 1, ncol = length(s2_list) + 1)

if(aln_mode == "global")
{
	align_matrix[1,1] <- 0
	for(i in 2:(length(s1_list)+1))
	{
		ix[i,1] <- g_o_int + (i-2) * g_e_int
	}	
	for(j in 2:(length(s2_list)+1))
	{
		iy[1,j] <- g_o_int + (j-2) * g_e_int
	}
	for(i in 2:(length(s1_list)+1))
	{
		for(j in 2:(length(s2_list)+1))
		{
			mat <- align_matrix[i - 1, j - 1] + score_schem_matrix[s1_list[i - 1], s2_list[j - 1]]
			temp_ix <- ix[i - 1, j - 1] + score_schem_matrix[s1_list[i - 1], s2_list[j - 1]]
			temp_iy <- iy[i - 1, j - 1] + score_schem_matrix[s1_list[i - 1], s2_list[j - 1]]
			align_matrix[i, j] <- max(mat, max(temp_ix, temp_iy))
			ix[i, j] <- max((align_matrix[i - 1, j] + g_o_int), (ix[i - 1, j] + g_e_int))
			iy[i, j] <- max((align_matrix[i, j - 1] + g_o_int), (ix[i, j - 1] + g_e_int))
		}
	}
	#trace-back
	pattern <- ""
	subject <- ""
	pointer <- c((length(s1_list)+1),(length(s2_list)+1))
	while( !((pointer[1] == 1) && (pointer[2] == 1)) )
	{
		if(pointer[1] == 1)
		{
			pattern <- paste(pattern, "-", sep = "")
			subject <- paste(subject, s2_list[pointer[2] - 1], sep = "")
			pointer <- c(pointer[1], pointer[2] - 1)
		}
		else if(pointer[2] == 1)
		{
			pattern <- paste(pattern, s1_list[pointer[1] - 1], sep = "")
			subject <- paste(subject, "-", sep = "")
			pointer <- c(pointer[1] - 1, pointer[2])
		}
		else
		{
			mat <- align_matrix[pointer[1], pointer[2]]
			top <- ix[pointer[1], pointer[2]]
			left <- iy[pointer[1], pointer[2]]
			choose_path = max(mat, max(top, left))
			if(choose_path == mat)
			{
				pattern <- paste(pattern, s1_list[pointer[1] - 1], sep = "")
				subject <- paste(subject, s2_list[pointer[2] - 1], sep = "")
				pointer <- c(pointer[1] - 1, pointer[2] - 1)
			}
			else if(choose_path == top)
			{
				pattern <- paste(pattern, s1_list[pointer[1] - 1], sep = "")
				subject <- paste(subject, "-", sep = "")
				pointer <- c(pointer[1] - 1, pointer[2])
			}
			else
			{
				pattern <- paste(pattern, "-", sep = "")
				subject <- paste(subject, s2_list[pointer[2] - 1], sep = "")
				pointer <- c(pointer[1], pointer[2] - 1)
			}
		}
	}
}

if(aln_mode == "local")
{
	align_matrix[1:(length(s1_list)+1), 1] <- 0
	align_matrix[1, 1:(length(s2_list)+1)] <- 0
	for(i in 2:(length(s1_list)+1))
	{
		for(j in 2:(length(s2_list)+1))
		{
			mat <- align_matrix[i - 1, j - 1] + score_schem_matrix[s1_list[i - 1], s2_list[j - 1]]
			temp_ix <- ix[i - 1, j - 1] + score_schem_matrix[s1_list[i - 1], s2_list[j - 1]]
			temp_iy <- iy[i - 1, j - 1] + score_schem_matrix[s1_list[i - 1], s2_list[j - 1]]
			align_matrix[i, j] <- max(mat, max(temp_ix, max(temp_iy, 0)))
			ix[i, j] <- max((align_matrix[i - 1, j] + g_o_int), (ix[i - 1, j] + g_e_int))
			iy[i, j] <- max((align_matrix[i, j - 1] + g_o_int), (ix[i, j - 1] + g_e_int))
		}
	}
	#trace-back
	pattern <- ""
	subject <- ""
	count_max = (length(which(align_matrix == max(align_matrix), arr.ind = TRUE)) / 2)
	tempx = which(align_matrix == max(align_matrix), arr.ind = TRUE)[1]
	tempy = which(align_matrix == max(align_matrix), arr.ind = TRUE)[count_max+1]
	temp_pointer = c(tempx, tempy)
	for(i in 2:count_max)
	{
		choose_pointer = c(which(align_matrix == max(align_matrix), arr.ind = TRUE)[i], which(align_matrix == max(align_matrix), arr.ind = TRUE)[i+count_max])
		if(((temp_pointer > choose_pointer)[1] == FALSE) && ((temp_pointer > choose_pointer)[2] == FALSE))
		{
			temp_pointer <- choose_pointer
		}
		else if(((temp_pointer > choose_pointer)[1] == FALSE) && ((temp_pointer > choose_pointer)[2] == TRUE))
		{
			temp_pointer <- temp_pointer
		}
		else if(((temp_pointer > choose_pointer)[1] == TRUE) && ((temp_pointer > choose_pointer)[2] == FALSE))
		{
			temp_pointer <- choose_pointer
		}
		else
		{
			temp_pointer <- temp_pointer
		}
	}
	pointer <- temp_pointer
	while( !((pointer[1] == 1) || (pointer[2] == 1)) )
	{
		mat <- align_matrix[pointer[1], pointer[2]]
		top <- ix[pointer[1], pointer[2]]
		left <- iy[pointer[1], pointer[2]]
		choose_path = max(mat, max(top, left))
		if(choose_path == mat)
		{
			pattern <- paste(pattern, s1_list[pointer[1] - 1], sep = "")
			subject <- paste(subject, s2_list[pointer[2] - 1], sep = "")
			pointer <- c(pointer[1] - 1, pointer[2] - 1)
		}
		else if(choose_path == top)
		{
			pattern <- paste(pattern, s1_list[pointer[1] - 1], sep = "")
			subject <- paste(subject, "-", sep = "")
			pointer <- c(pointer[1] - 1, pointer[2])
		}
		else
		{
			pattern <- paste(pattern, "-", sep = "")
			subject <- paste(subject, s2_list[pointer[2] - 1], sep = "")
			pointer <- c(pointer[1], pointer[2] - 1)
		}
	}
}
pattern <- reverse(pattern)
subject <- reverse(subject)
r <- BStringSet( c(pattern, subject) )
names(r) <- seq_name
writeXStringSet(r, o_f, format="fasta")
