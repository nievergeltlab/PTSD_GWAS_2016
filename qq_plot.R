args <- commandArgs(trailingOnly = TRUE)
pvals <- args[1]
outfilename <- args[2]
errorbars <- args[3]

data <- read.table(pvals,header=F,stringsAsFactors=F,skip=1,nr=20000000,colClasses="numeric")

print(head(data))

unadj_filtered <- sort(data[,1])
UNADJ <- -log(unadj_filtered,10)
QQ <- -log(ppoints(length(UNADJ)),10)



GCfactor= round(median(qchisq(data[,1],1,lower.tail=F),na.rm=T)/.455,3)

png(paste(outfilename,'.png', sep=''))
par(bty='l')

plot(c(0,max(QQ)), c(0,max(UNADJ)), xlab='Expected -log10(p)', ylab='Observed -log10(p)', col='white', cex=1.3, cex.axis=1.2, cex.lab=1.5,pch=20)

if(errorbars == 1)
{
    #code for error bars
    ranks <- c(1:length(QQ))
    CIlower <- qbeta(.025, ranks, length(QQ)-ranks +1)
    CIupper <- qbeta(.975, ranks, length(QQ)-ranks +1)
    plotCIlower <- -log(CIlower,10)
    plotCIupper <- -log(CIupper,10)
    segments(x0=QQ,x1=QQ, y0=plotCIlower,y1=plotCIupper,lwd=2,col='grey')
}

abline(0,1,col='red', lwd=2)
points(QQ, UNADJ, ,pch=20,col='blue')

legend('topleft', paste('GC Lambda =', GCfactor),  bty='n', cex=1.5, xjust=1)

dev.off()