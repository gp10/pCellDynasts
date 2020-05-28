library(foreach)
library(iterators) #if not installed
library(parallel) #if not installed
library(doMC)
library(modes)

## A set of convenience functions to simplify running multimode in batch

registerDoMC(4) #Change this to match the number of cores on your machine

library(multimode)

# Calling a single specific modality test to be applied on given data (with a particular significance level):
modesearch <- function(dataset,method="HH", p=0.05) {
workingData <- data.frame(dataset)
positive <- foreach(d=workingData,.combine=rbind) %dopar% {
r <- modetest(log2(d), method=method)
s <- if (r$p<p) { 1 } else { 0 }
s
}
result <- sum(positive)/ncol(dataset)
print(sprintf("%f (%d/%d)", result, sum(positive),ncol(dataset)))
return(result)
}

# Calling a single specific modality test to be applied on all data:
searchEverything <- function(method="HH",nmode=1) {
workingData <- data.frame(everything)
result <- foreach(d=names(workingData), .combine=rbind) %dopar% {
r <- modetest(log2(workingData[[d]]), method=method, mod0=nmode)
sprintf("%s with %s = %f", d, method, r$p)
}
result
}

# Calling a given set of modality tests to be applied on given data:
matrixSearch <- function(data,p=0,nmode=1) {
workingData <- data.frame(data)
result <- c() 
if (nmode==1) {
methods <- c("HH","SI","CH","HY","FM","ACR")
} else {
methods <- c("SI","FM","ACR")
}
for (m in methods){
	print(m)
	positive <- foreach(d=workingData,.combine=rbind) %dopar% {
		r <- modetest(log2(d), method=m, mod0=nmode) # modality test is called
		s <- if (p) { if (r$p<p) { 1 } else { 0 } } else { r$p }
		#s <- data.frame(s)
		#rownames(s) <- methods
		}
result <- cbind(result,positive)
}
result <- data.frame(result)
colnames(result) <- methods
rownames(result) <- names(workingData)
result
}

# Calling the calculation of bimodality coefficient for given data:
bimCSearch <- function(data,threshold=0) {
workingData <- data.frame(data)
result <- foreach(d=workingData,.combine=rbind) %dopar% {
raw <- bimodality_coefficient(log2(na.exclude(d))) # bimodality coefficient is calculated
r <- if (threshold) {if (raw>threshold) { 1 } else { 0 } } else { raw }
}
rownames(result) <- names(workingData)
result
}
