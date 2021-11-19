# Calculation of F stat 

# Chr	Position	Marker	        A1	A2	A1_Freq	    MAF	Test	    N	   Beta	     SE	         t	        P	        rsid
# 21	48025097	21:48025097:A:G	A	G	0.307479	0.307479	ADD	769	-0.461949	0.0537992	-8.58655	4.96E-17	rs8128872


# effect allele freq: 0.307479
# minor allele freq: 0.307479
MAF <- 0.307479
beta <- -0.461949
N <- 769
k <- 1


r2 <- 2 * MAF * (1-MAF) * beta^2

# r2 = 0.09

F <- ((N - k - 1) / k) * (r2 / (1 - r2))

# F = 76.67



