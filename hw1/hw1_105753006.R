#require pakage expm
require(expm)

#fix round
round2 = function(x, n) {
  posneg = sign(x)
  z = abs(x)*10^n
  z = z + 0.5
  z = trunc(z)
  z = z/10^n
  z*posneg
}

# read PAM1 from data, get header and skip first line
pam1 <- read.table("pam1.txt", header = T, skip = 1)

# construct PAM250 from PAM1
pam1_matrix <- t(matrix(unlist(pam1), ncol = 20, byrow = TRUE))
pam1_config <- pam1_matrix / 10000
pam250_noScale <- pam1_config %^% 250 * 100
for(x in c(1:20))
{
	sum_temp <- sum(pam250_noScale[,x])
	for(y in c(1:20))
	{
		pam250_noScale[y,x] <- pam250_noScale[y,x] / sum_temp * 100
	}
}
pam250 <- round2(pam250_noScale, 0)
# output PAM250 as a file
rownames(pam250) <- rownames(pam1)
colnames(pam250) <- rownames(pam1)
write.table(pam250, "pam250.txt", row.names = TRUE, col.names = NA, quote = FALSE, sep = '\t')
