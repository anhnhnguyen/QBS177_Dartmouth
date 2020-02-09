# Transform the subpopulation data to 21 countries:
chn <- (fulldat[,5] + fulldat[,6])/2
ind <- (fulldat[,11] + fulldat[,14])/2
nga <- (fulldat[,25] + fulldat[,8])/2
usa <- (0.777*fulldat[,4] + 0.132*fulldat[,2] + 0.053*(chn + ind + fulldat[,16] + fulldat[,15])/4)/(0.777+0.132+0.053)
fulldat <-fulldat[,c(3,1,5,7,9,12,11,24,15,17,19,25,21,20,22,18,13,23,10,4,16)]
fulldat[,3] <- chn
fulldat[,7] <- ind
fulldat[,12] <- nga
fulldat[,20] <- usa

# Read in the smoking data.
yy <- read.delim("smoking_outcome.txt")

# Try to write a R code to calculate the correlation between yy and the first 10000 loci. Use proc.time() before and after your code to report the computation time.
y <- yy[,2]
plong <- corlong <- NULL

proc.time()
for (i in 1:10000){
  if(var(fulldat[i,])){   # remove loci with all minor allele frequncy 0
    fit <- lm(y~fulldat[i,]) # run a linear regression 
    corlong[i] <- cor(fulldat[i,], y)
    plong[i] <- summary(fit)$coef[2,4] #output slope from linear regression
  }
}
proc.time()

# Estimate the total computation time if you would apply the code to all loci’s in Chromosome 22.
(20.020-8.952)*dim(fulldat)[1]/10000

# Try to use vector operation to rewrite the code and speed it up. You can use the following code or write your own.
n <- dim(yy)[1] # n is the sample size
pvalue <- matrix(0, length(chr), 1) # initial a variable to store p values
corvalue <- matrix(0, length(chr), 1) # initial a variable to store correlation
colnames(pvalue) <- c("prevave")
colnames(corvalue) <- c("prevave")

proc.time()
asd <- suppressWarnings(apply(t(fulldat), 2 , cor , y)) # obtain correlation using apply
asd1 <- 1 - asd^2 # mid step to calculate p value
asd2 <- sqrt(n-2)*asd/sqrt(asd1) # this is the T statistics
pvalue[,1] <- 2*(1-pt(abs(asd2), (n-2))) # transform T statistics to p-value
corvalue[,1] <- asd # pass correlation to corvalue
proc.time()

# Check if you got the same p value and correlation using the two methods.
abs(sum(plong - pvalue[1:10000], na.rm = T))
abs(sum(corlong - corvalue[1:10000], na.rm = T))

# Draw a mahatten Plot using “manhattan” function in “qqman”. You will need specify the SNP name, chromosome number and location for the plot. 
install.packages("qqman")
library(qqman)
out <- data.frame(chr, location, pvalue)
names(out) <- c("CHR", "BP", "P")
SNP <- paste("SNP", 1:dim(out)[1])
out <- cbind(SNP, out)
manhattan(out[,1:4], ylim = c(0,10), xlim = c(0, max(out[,3])))





