png('esophagus_dilution.png')
hist(log2(eso00$mouse6c),xlim=c(-6,8),freq=F,col=rgb(1,0,0,0.25),main="Esophageal histone dilution",xlab="log2(Histone Intensity)",border=NA)
hist(log2(eso00$mouse6d),xlim=c(-6,8),freq=F,col=rgb(1,0,0,0.25),add=T,border=NA)
hist(log2(eso00$mouse1k),xlim=c(-6,8),freq=F,col=rgb(1,0,0,0.25),add=T,border=NA)

hist(log2(eso07$mouse1c),xlim=c(-6,8),freq=F,col=rgb(0,1,0,0.25),add=T,border=NA)
hist(log2(eso07$mouse1i),xlim=c(-6,8),freq=F,col=rgb(0,1,0,0.25),add=T,border=NA)
hist(log2(eso07$mouse1j),xlim=c(-6,8),freq=F,col=rgb(0,1,0,0.25),add=T,border=NA)

hist(log2(eso12$mouse6e),xlim=c(-6,8),freq=F,col=rgb(0,0,1,0.25),add=T,border=NA)
hist(log2(eso12$mouse6g),xlim=c(-6,8),freq=F,col=rgb(0,0,1,0.25),add=T,border=NA)
hist(log2(eso12$mouse1a),xlim=c(-6,8),freq=F,col=rgb(0,0,1,0.25),add=T,border=NA)

hist(log2(eso18$mouse1d),xlim=c(-6,8),freq=F,col=rgb(1,1,0,0.25),add=T,border=NA)
hist(log2(eso18$mouse1e),xlim=c(-6,8),freq=F,col=rgb(1,1,0,0.25),add=T,border=NA)
legend("topleft", c("0 day", "7 day", "12 day", "18 day"), col=c("red","green","blue","yellow"), lwd=10)
dev.off()
