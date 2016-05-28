args <- commandArgs(TRUE)


process_gene=function(gene_name,fusion_location_in_gene,sample_name, output_file)
{
# Read input file for the given gene
file_in=paste(gene_name,"_data.tsv", sep="")
x<-read.table(file_in,header=F, sep="\t", row.names=1)
# Get all sample names from TCGA file
y=read.delim("HiSeqV2_exon",sep="\t", header=F, nrows=1, row.names=1)
z=c()
for(i in 1:length(y[1,]))
{z<-c(z,as.character(y[1,i]))}
colnames(x)<-z
if(sum(colnames(x)==sample_name)==0) # if sample doesn't exist
{
  print("Invalid sample name")
  dev.off()
  write.table("",args[6], append=TRUE,sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\n")
  quit(save="no")
}

# Add extra column for exon number and find the breakpoint exon/intron
colnames(x)[length(x)]<-"Exon_number"
z<-c()
exon_breakpoint<-0
for(i in rownames(x))
{z<-c(z,as.numeric(strsplit(i,"[:-]")[[1]][2:3]))}
for(i in 1:(length(z)-1))
{
  if(as.numeric(fusion_location_in_gene)>z[i]&&as.numeric(fusion_location_in_gene)<=z[i+1])
  {
    exon_breakpoint=(i+1)/2
  }
}
if(exon_breakpoint==0) # means the breakpoint is outside exonic regions
{
  print("break position out of exon boundaries. Try a different method")
  gene_name=paste(gene_name,"outbreak",sep="_")
  if(as.numeric(fusion_location_in_gene)<z[1]) # 5' to gene
   {
     exon_breakpoint=1
   }
  else if(as.numeric(fusion_location_in_gene)>z[length(z)]) # 3' to gene
   {
     exon_breakpoint=length(x[,1])
   }

}
if(strsplit(rownames(x)[1],"[:]")[[1]][3]=='-' ) # Is the gene on the positive or negative strand? If yes, reverse exon numbers. Then calculate log2 value
{
  x[,length(x)]<-c(length(x[,1]):1)
  log2_diff<-log2(mean(x[1:floor(exon_breakpoint),colnames(x)==sample_name])/mean(x[ceiling(exon_breakpoint):length(x[,1]),colnames(x)==sample_name]))
} else {  
  x[,length(x)]<-c(1:length(x[,1]))
log2_diff<-log2(mean(x[ceiling(exon_breakpoint):length(x[,1]),colnames(x)==sample_name])/mean(x[1:floor(exon_breakpoint),colnames(x)==sample_name]))
}
# Calculate mean absolute expression
mean_abs_expression=mean(rowMeans(x[,-length(x)]))

# Normalise exons
for(i in 1:(length(x[,1])-1)){
  x[i,-length(x)]<-x[i,-length(x)]-rowMeans(x[,-length(x)])[i]
}

# Add empty datapoints to pass to xlim such that the fusion is always in the centre of the plot
if(x[exon_breakpoint,length(x)]-1 > length(x[,1])-x[exon_breakpoint,length(x)])
{
  xl=c(1,2*x[exon_breakpoint,length(x)]-1) #xlim
  xat=c(1:length(x[,1])) # where to put tickmarks on the plot
  labx=c(1:length(x[,1])) # labels for tickmarks
}
else
{
  xl=c(2*x[exon_breakpoint,length(x)]-length(x[,1]),length(x[,1]))
  xat=c(1:length(x[,1]))
  labx=c(1:length(x[,1]))
}

# Build plot area with axes names and limits
plot(x=x[,length(x)],y=x[,1], type='n', xlim=xl, ylim=c(min(x[,-length(x)]),1.4*max(x[,-length(x)])),  cex.axis=0.8, xlab="exon number", ylab="relative exon counts", main=gene_name,xaxt='n')
axis(1,at=xat,labels=labx,cex.axis=0.8)

# Plot cohort expressions
for(i in 1:(length(x[1,])-1))
{lines(x[,length(x)],x[,i], type='l', col="grey")}

# Plot mean expression of cohort
lines(x[,length(x)],rowMeans(x[,-length(x)]), type='l', col='black', lwd=4)
# Plot area before breakpoint
lines(x[1:floor(exon_breakpoint),length(x)],x[1:floor(exon_breakpoint),colnames(x)==sample_name], type='l', col='red', lwd=4)
# Plot area after breakpoint
lines(x[ceiling(exon_breakpoint):length(x[,1]),length(x)],x[ceiling(exon_breakpoint):length(x[,1]),colnames(x)==sample_name], type='l', col='blue', lwd=4)

# fill in blank spaces around breakpoint
lines(c(x[floor(exon_breakpoint),length(x)],(x[floor(exon_breakpoint),length(x)]+x[ceiling(exon_breakpoint),length(x)])/2),c(x[floor(exon_breakpoint),colnames(x)==sample_name],(x[floor(exon_breakpoint),colnames(x)==sample_name]+x[ceiling(exon_breakpoint),colnames(x)==sample_name])/2),type='l', col='red', lwd=4)
lines(c((x[floor(exon_breakpoint),length(x)]+x[ceiling(exon_breakpoint),length(x)])/2,x[ceiling(exon_breakpoint),length(x)]),c((x[floor(exon_breakpoint),colnames(x)==sample_name]+x[ceiling(exon_breakpoint),colnames(x)==sample_name])/2,x[ceiling(exon_breakpoint),colnames(x)==sample_name]),type='l', col='blue', lwd=4)
# Plot vertical line at the beakpoint
abline(v=((x[floor(exon_breakpoint),length(x)]+x[ceiling(exon_breakpoint),length(x)])/2),lwd=2,col='black')

# Perform t tests on areas before and after the breakpoint. If it's an error, set values to NA
cc<-NULL
cc<-try(a<-t.test(x[1:floor(exon_breakpoint),colnames(x)==sample_name],rowMeans(x)[1:floor(exon_breakpoint)], paired=TRUE),silent=TRUE)
if(length(cc)==1)
{
  a<-list(p.value="NA",estimate="NA")
}
cc<-NULL
cc<-try(b<-t.test(x[ceiling(exon_breakpoint):length(x[,1]),colnames(x)==sample_name],rowMeans(x)[ceiling(exon_breakpoint):length(x[,1])], paired=TRUE),silent=TRUE)
if(length(cc)==1)
{
  b<-list(p.value="NA",estimate="NA")
}


# Insert text onto the plot for p values (before and after breakpoint). Fill output table with data on gene
if(strsplit(rownames(x)[1],"[:]")[[1]][3]=='-' )
{
  text(xl[1],1.3*max(x[,-length(x)]), paste("Upstream     p=",b$p.value, sep=""), pos=4, cex=0.6)
  text(xl[1],1.15*max(x[,-length(x)]),paste("Downstream p=",a$p.value, sep=""), pos=4, cex=0.6)
  write.table(cbind(gene_name,fusion_location_in_gene,b$estimate,b$p.value,a$estimate,a$p.value, log2_diff, mean_abs_expression),output_file, append=TRUE,sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\t")
} else {  
  text(xl[1],1.3*max(x[,-length(x)]), paste("Upstream     p=",a$p.value, sep=""), pos=4, cex=0.6)
  text(xl[1],1.15*max(x[,-length(x)]),paste("Downstream p=",b$p.value, sep=""), pos=4, cex=0.6)
  write.table(cbind(gene_name,fusion_location_in_gene,a$estimate,a$p.value,b$estimate,b$p.value, log2_diff,mean_abs_expression),output_file, append=TRUE,sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\t")
}

# Add text to plot for Log2 ratio and MAE
mtext(paste("log2(Downstream/Upstream) = ",signif(log2_diff,4)," and Mean Abs. Expression = ",signif(mean_abs_expression,4),sep=""),side=3, cex=0.75)
}
# END OF FUNCTION



# Add newline character to output file. This is important if the previous write operation terminated improperly
write.table("",args[6], append=TRUE,sep="\t",quote=FALSE, row.names=FALSE, col.names=FALSE, eol="\n")
sample_name=args[5]
# Produce a jpeg device to plot to
jpeg(paste(args[1],"_",args[2],"_",args[3],"_",args[4],"_",sample_name,".jpeg", sep=""))
par(mfrow=c(2,1))
process_gene(gene_name=args[1],fusion_location_in_gene=args[2],sample_name=sample_name, output_file=args[6])
process_gene(gene_name=args[3],fusion_location_in_gene=args[4],sample_name=sample_name,output_file=args[6])
write.table(sample_name,args[6], append=TRUE,sep="",quote=FALSE, row.names=FALSE, col.names=FALSE, eol="")
dev.off()

