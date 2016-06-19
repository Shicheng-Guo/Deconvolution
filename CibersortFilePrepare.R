options("scipen"=100, "digits"=4)
args <- commandArgs(trailingOnly = TRUE)
pure=args[1]
mix=args[2]
class=args[3]
A = read.table(pure,T,row=1)
B = read.table(mix,T,row=1)
C = read.table(class,F,row=1)
A = A # *100
B = B # *100
for( i in 1:nrow(C) ){
        cur_ref = which(C[i,]==1)
        average_values = apply(A[,cur_ref],1,function(x) mean(x,na.rm=T))
        for( j in 1:nrow(A) ){
                na_pos = which(is.na(A[j,]) & C[i,]==1)
                A[j, na_pos] = average_values[j]
        }
}
write.table(cbind(data.frame(probe_id=rownames(A)), A), file=paste(pure,".mod",sep=""), quote=F,sep="\t",row=F)
write.table(cbind(data.frame(probe_id=rownames(B)), B), file=paste(mix,".mod",sep=""), quote=F,sep="\t",row=F)
