# Download and read in the melanoma dataset into R:
mel <- read.table("PCA.example.melanoma.txt",header=TRUE)

# Fill the missing values using column means.   
for(i in 1:ncol(mel)){
  mel[is.na(mel[,i]), i] <- round(mean(mel[,i], na.rm = TRUE))
}
sum(is.na(mel))

# Perform principle component analysis using spectral decomposition of covariance matrix.
mel.pc.eigen.cor <- prcomp(mel)

# Check the output object and identify the principle components.
names(mel.pc.eigen.cor)
mel.scores.eigen <- mel.pc.eigen.cor$x

# Reconstruct the first principle component using loadings.
pc1 <- mel.pc.eigen.cor$rotation[, 1] %*% t(mel)
cor(t(pc1), mel.scores.eigen[, 1])

# Reconstruct the covariance matrix using loadings and eigenvalues.
eigenvalue <- diag(mel.pc.eigen.cor$sdev^2)
newdata <- (mel.pc.eigen.cor$rotation) %*% eigenvalue %*% t(mel.pc.eigen.cor$rotation)
covdata <- cov(mel)
sum(abs(newdata - covdata))

# Check the variance for the first 15 principle components.
screeplot(mel.pc.eigen.cor, type = "l", npcs = 15, main = "Screeplot of the first 15 PCs")
abline(h = 1, col="red", lty=5)
legend("topright", legend=c("Eigenvalue = 1"),
       col=c("red"), lty=5, cex=0.6)

# Plot cumulative variance plot for the first 15 principle components.
cumpro <- cumsum(mel.pc.eigen.cor$sdev^2 / sum(mel.pc.eigen.cor $sdev^2))
plot(cumpro[0:15], xlab = "PC #", ylab = "Amount of explained variance", main = "Cumulative variance plot")

# Visualize the first two principle components.
plot(mel.scores.eigen[,1],mel.scores.eigen[,2], xlab="PCA 1", ylab="PCA 2")

# The last 505 samples have ethnicity information. Replot with color coding.
plot(mel.scores.eigen[1:2830,1],mel.scores.eigen[1:2830,2],type="p",col="black",pch=1,xlab="PC1", ylab="PC2", xlim = range(mel.scores.eigen[,1]), ylim = range(mel.scores.eigen[,2]));
points(mel.scores.eigen[2831:2995,1],mel.scores.eigen[2831:2995,2],col="red",pch=20);
points(mel.scores.eigen[2996:3132,1],mel.scores.eigen[2996:3132,2],col="green",pch=20)
points(mel.scores.eigen[3133:3335,1],mel.scores.eigen[3133:3335,2],col="blue",pch=20)
title(main="Principal Components Analysis using princomp in R", col.main="black", font.main=1,outer=T)

# Perform PCA using the singular value decomposition
mel1 <- scale(mel, scale = FALSE)
mel.pc.scale.svd <- svd(mel1)

# Check the output object and compare with decomposition result. 
names(mel.pc.scale.svd)
head(mel.pc.scale.svd$d^2/(nrow(mel) - 1))
head(mel.pc.eigen.cor$sdev^2)

mel.pc.scale.svd$v[1:5,1:5]
mel.pc.eigen.cor$rot[1:5,1:5]

# Obtain principle components
pc1 <- (mel.pc.scale.svd$v[,1]%*%t(mel))
cor(t(pc1), mel.scores.eigen[,1])

# Calculate the centroids for these populations based on first 2 principle components. 
mel.euro <- apply(mel.scores.eigen[2831:2995,1:2], 2, mean)
mel.asian <- apply(mel.scores.eigen[2996:3132,1:2], 2, mean)
mel.african <- apply(mel.scores.eigen[3133:3335,1:2],2,mean)

# For the first 2830 samples with unknown ancestry, calculate the distance between each centroid and cluster each subject to the nearest cluster.
sample.vector <- rep("0",2830) 
for(i in 1:2830){
  mel.euro.dist <- sqrt((mel.scores.eigen[i,1] - mel.euro[1])^2 + (mel.scores.eigen[i,2] - mel.euro[2])^2)
  mel.asian.dist <- sqrt((mel.scores.eigen[i,1] - mel.asian[1])^2 + (mel.scores.eigen[i,2] - mel.asian[2])^2)
  mel.african.dist <- sqrt((mel.scores.eigen[i,1] - mel.african[1])^2 + (mel.scores.eigen[i,2] - mel.african[2])^2)
  sample.index <- which.min(c(mel.euro.dist, mel.asian.dist, mel.african.dist))
  if (sample.index == 1){
    sample.vector[i] <- "red"}
  if (sample.index == 2){
    sample.vector[i] <- "green"}
  else if (sample.index == 3){
    sample.vector[i] <- "blue"}
}
  
# Repeat 10 on 2830 sample color coded.
plot(mel.scores.eigen[1:2830,1],mel.scores.eigen[1:2830,2],type="p",col=sample.vector,pch=1,xlab="PC1", ylab="PC2", xlim = range(mel.scores.eigen[,1]), ylim = range(mel.scores.eigen[,2]));
points(mel.scores.eigen[2831:2995,1],mel.scores.eigen[2831:2995,2],col="red",pch=20);
points(mel.scores.eigen[2996:3132,1],mel.scores.eigen[2996:3132,2],col="green",pch=20)
points(mel.scores.eigen[3133:3335,1],mel.scores.eigen[3133:3335,2],col="blue",pch=20)
title(main="Principal Components Analysis using princomp in R", col.main="black", font.main=1,outer=T)

