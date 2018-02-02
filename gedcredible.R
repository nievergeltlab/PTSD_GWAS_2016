#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/r_2.11.0/bin/Rscript
args <- commandArgs(TRUE)

R2_threshold <- 0.4

file <- args[1]
#file <-"locus.txt"
dat <- read.table(file=file, header=T)

getCredibleSNP <- function(snp, logProb, threshold=0.99){

 prob <- exp(logProb) 
 prob_normed <- prob/sum(prob)

 prob_cumsum <- cumsum(sort(prob_normed, decreasing=TRUE))
 nSNP <- as.numeric(which.max(prob_cumsum>threshold ))

 if(nSNP<=0){ nSNP=1 }

 credible_set <- snp[order(prob_normed, decreasing=TRUE)[1:nSNP]]
 select <- rep(FALSE, length(snp))
 select[order(prob_normed, decreasing=TRUE)[1:nSNP]] <- T 
 ret <- list(nSNP = nSNP, 
	prob_normed = prob_normed, 
	prob_cumsum = prob_cumsum[rank(-prob_normed, ties.method="random")], 
	credible_set=credible_set, 
	select=select 
 )
}

select <- !is.na(dat$R2) & dat$R2 > R2_threshold

max_SCZ_P <- min(dat$P[select])
max_SCZ_snp <- as.character(dat$SNP[select][which.min(dat$P[select])])
print(max_SCZ_P)
print(max_SCZ_snp)
print(qchisq(dat$P[select], 1, low=F))
ret <- getCredibleSNP(as.character(dat$SNP[select]), qchisq(dat$P[select], 1, low=F)/2)

inCredible <- rep(NA, length(dat$P))
inCredible[match( dat$SNP[select], dat$SNP )] <- 0
inCredible[match( dat$SNP[select][ret$select], dat$SNP )] <- 1
prob_norm <- rep(NA, length(dat$P))
prob_norm[match(dat$SNP[select], dat$SNP )] <- ret$prob_normed
prob_cumsum <- rep(NA, length(dat$P))
prob_cumsum[match(dat$SNP[select], dat$SNP )] <- ret$prob_cumsum

result <- cbind(dat, inCredible=inCredible, probNorm=prob_norm, cumSum=prob_cumsum)

write.table(result, file=paste(file,"credible.txt",sep="_"), sep="\t", quote=F, col.names=T, row.names=F)
