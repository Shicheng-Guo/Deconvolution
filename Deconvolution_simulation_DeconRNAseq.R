
library("DeconRNASeq")
k=11   # 11 reference tissue
s=3    # each tissue own 3 samples
p=30   # number of tissue specific probe for each tissue
signatures<-c()
tissue=c("Brain","CCT","Colon","Esophagus","Heart","Intestine","Kidney","Liver","Lung","Stomach","WBC")

signatures<-abs(matrix(rnorm(p*k*s,0.1,0.1),p*k,k*s))
for(i in 1:k){
  signatures[(p*(i-1)+1):(p*i),(s*(i-1)+1):(s*i)]=matrix(rnorm(p*s,0.9,0.1),p,s) 
}
signatures[signatures>1]<-1
colnames(signatures)=rep(tissue,each=3)
rownames(signatures)=paste("cg00",1:nrow(signatures),sep="")

barplot(colMeans(signatures))
barplot(rowMeans(signatures))
heatmap(signatures,Rowv=NA,Colv=NA)
? heatmap
gsi<-data.frame(GSI(signatures))
group<-as.character(unique(gsi$group))

# tissue specific MHL selectin
rlt<-c()
rank<-c(rep(20,length(group)))
for (i in 1:length(group)){
  subset=gsi[which(gsi$group==group[i]),]
  subset=subset[order(subset[,3],decreasing=T)[1:rank[i]],]
  rlt<-rbind(rlt,subset)
}

# prepare signatures for QR
signaturesQR<-AverageWithGroup(signatures)  # QR
# prepare signatures for cibersort
signaturesCS<-(signatures)  # QR

# Simulation mixture for 2 componment
VirtualMatrix<-c()
for(F1 in seq(0,1,by=0.05)){
  VirtualSample<-F1*(signaturesQR[,grep("CCT",colnames(signaturesQR))])+(1-F1)*(signaturesQR[,grep("WB",colnames(signaturesQR))])
  VirtualMatrix<-cbind(VirtualMatrix,VirtualSample)
}
VirtualMatrix<-data.frame(VirtualMatrix)
colnames(VirtualMatrix)=paste("VM",seq(0,1,by=0.05),sep="")
dim(VirtualMatrix)

DeconData<-data.frame(VirtualMatrix,signaturesQR)
VirtualMatrix=data.frame(DeconData[,grep("VM",colnames(DeconData))])
Signatures=data.frame(DeconData[,-grep("VM",colnames(DeconData))])

library("car")
library("DeconRNASeq")
plot(density(as.numeric(unlist(VirtualMatrix))))
plot(density(as.numeric(unlist(logit(VirtualMatrix)))))

VirtualMatrix=logit(VirtualMatrix)
Signatures=logit(Signatures)
Rlt<-DeconRNASeq(VirtualMatrix, Signatures,checksig=FALSE,known.prop = F, use.scale = TRUE, fig = TRUE)
rlt<-Rlt$out.all
output<-data.frame(CCTInput=seq(0,1,by=0.05),rlt)
plot(output[,1],output[,3])
lines(x = c(0,1), y = c(0,1))
plot(output[,1],output[,3]-output[,1],xlab="Simulated contribution (%)",ylab="Prediction error (%)",pch=19,cex=2,col="blue",lty=2)


# Figure parepartion
library("ggplot2")
Fig<-data.frame(seq(0,1,by=0.05),Rlt$out.all[,grep("CCT",colnames(Rlt$out.all))],Rlt$out.all[,grep("WB",colnames(Rlt$out.all))])
colnames(Fig)<-c("Simulated","CCT","WB")
Fig<-100*Fig
c <- ggplot(Fig, aes(Simulated,CCT))
c <- c + xlab("Simulated contribution (%)")+ylab("Predicted contribution (%)")
c <- c + stat_smooth(se = TRUE,n=10,size=1.5) + geom_point(size=3)
c <- c + stat_smooth(se = TRUE,n=10,size=1.5) + geom_point(size=3)
c <- c + theme_bw()
c <- c + theme(axis.text=element_text(size=10), axis.title=element_text(size=14,face="bold"))
print(c)
ggsave(c,file="coloncancer-deconvolution.simultaion-2.pdf")
abline()
