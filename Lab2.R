library(tidyverse)

load('dbp.R')
ls()
dbp[1:5,]

## Run a logistic regression analysis of the affection status regressed on the genotype
## of marker rs1112, using the data in the data frame dbp. Assign the results from 
## the regression analysis to the new object result.snp12:

result.snp12 <- glm(affection ~ rs1112, family = binomial(link = 'logit'), data = dbp)

# Print the results of the regression analysis
print(result.snp12)
print(class(result.snp12))
print(summary(result.snp12))

# Likilihood-ratio test (LRT)
dev.geno = anova (result.snp12, test="Chi")
lrt.pvalue = pchisq(dev.geno[dim(dev.geno)[1],"Deviance"], 
                    df=2, ncp=0, FALSE) 
print ( lrt.pvalue )

# Calculate the odds ratios for the genotypes and their confidence intervals
print ( summary(result.snp12)$coefficients )
snp.beta = summary(result.snp12)$coefficients[2:3,1]
print ( snp.beta )
print ( exp(snp.beta) )

ci = confint (result.snp12)
print (ci) 
print ( exp(ci) )

# Run regression with rs1112 as numeric
snp.data = dbp[,c("affection", "rs1112")]
summary(snp.data)

snp.data[,"rs1112"] <- as.numeric(snp.data[,"rs1112"]) - 1
summary(snp.data)
result.all = glm (affection ~ rs1112, family=binomial("logit"), 
                  data=snp.data)
dev.all    = anova (result.all, test="Chi")
lrt.pvalue = pchisq(dev.all[dim(dev.all)[1],"Deviance"], 
                    df=1, ncp=0, FALSE) 
print ( lrt.pvalue )
summary(result.all)   
print(dev.all) 

# Adjust regression 

snp.data = dbp[,c("affection", "trait","sex", "age", "rs1112", "rs1117")]
summary(snp.data)

snp.data[,"rs1112"] <- as.numeric(snp.data[,"rs1112"]) - 1
snp.data[,"rs1117"] <- as.numeric(snp.data[,"rs1117"]) - 1

result.adj = glm (affection ~ sex + rs1112      , family=binomial("logit"), 
                  data=snp.data)
summary(result.adj)

result.adj = glm (affection ~ age + rs1112      , family=binomial("logit"), 
                  data=snp.data)
summary(result.adj)

result.adj = glm (affection ~ sex + age + rs1112, family=binomial("logit"), 
                  data=snp.data)
summary(result.adj)

result.adj = glm (affection ~ rs1117 + rs1112, family=binomial("logit"), 
                  data=snp.data)
summary(result.adj)
anova (result.adj, test="Chi")

result.adj = glm (affection ~ rs1112 + rs1117, family=binomial("logit"), 
                  data=snp.data)
summary(result.adj)
anova (result.adj, test="Chi")

result.adj = lm (trait ~ rs1112      , data=snp.data)
summary(result.adj)

result.adj = lm (trait ~ sex + rs1112, data=snp.data)
summary(result.adj)

result.inter = glm (affection ~ sex * rs1112, family=binomial("logit"), 
                    data=snp.data)
summary(result.inter)

result.inter = glm (affection ~ age * rs1112, family=binomial("logit"), 
                    data=snp.data)
summary(result.inter)

result.inter = glm (affection ~ rs1112 * rs1117, family=binomial("logit"), 
                    data=snp.data)
summary(result.inter)


