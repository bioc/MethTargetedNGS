###### Sequence alignment and creating matrix for further analysis#############
methAlign <-function(sequence_fasta,ref_seq,sub_mat=FALSE,align_type="local")
{  t1 <- Sys.time()
   seq_list <- read.fasta(sequence_fasta,as.string=TRUE,set.attributes=FALSE)
   ref_list <-read.fasta(ref_seq,as.string=TRUE,set.attributes=FALSE)
   cglist <- unlist(str_locate_all(toString(ref_list[1]),"c(-)*g"))
   final_mat <- matrix(nrow=length(seq_list),ncol=length(cglist)/2,dimnames=list(c(names(seq_list)),c(1:(length(cglist)/2))),"-")
   if(sub_mat==FALSE)
   {
     mat <- nucleotideSubstitutionMatrix(match = 6, mismatch = -18, baseOnly = TRUE)
     mat[2,4] <- 3
     mat[4,2] <- 3
     mat <- rbind(mat,3)
     mat <- cbind(mat,3)
     colnames(mat)[5] <- "N"
     rownames(mat)[5] <- "N" }
   else
   { mat <- sub_mat}
   for (i in 1:length(seq_list))
   { p_align1 <- pairwiseAlignment(DNAString(toString(seq_list[i])),DNAString(toString(ref_list[1])),type=align_type,gapOpening=0,gapExtension=-1,substitutionMatrix=mat)
     p_align2 <- pairwiseAlignment(reverseComplement(DNAString(toString(seq_list[i]))),DNAString(toString(ref_list[1])),type=align_type,gapOpening=0,gapExtension=-1,substitutionMatrix=mat)
     if (pairwiseAlignment(DNAString(toString(seq_list[i])),DNAString(toString(ref_list[1])),type=align_type,gapOpening=0,gapExtension=-1,substitutionMatrix=mat,scoreOnly=TRUE) >= pairwiseAlignment(reverseComplement(DNAString(toString(seq_list[i]))),DNAString(toString(ref_list[1])),type=align_type,gapOpening=0,gapExtension=-1,substitutionMatrix=mat,scoreOnly=TRUE)){  p_align=p_align1 }
     if (pairwiseAlignment(DNAString(toString(seq_list[i])),DNAString(toString(ref_list[1])),type=align_type,gapOpening=0,gapExtension=-1,substitutionMatrix=mat,scoreOnly=TRUE) < pairwiseAlignment(reverseComplement(DNAString(toString(seq_list[i]))),DNAString(toString(ref_list[1])),type=align_type,gapOpening=0,gapExtension=-1,substitutionMatrix=mat,scoreOnly=TRUE)) {  p_align=p_align2 }
     cglist <- unlist(str_locate_all(toString(subject(p_align)),"C(-)*G"))
     final_cglist <- unlist(cglist[1:(length(cglist)/2)])    
     cur_seq <- unlist(strsplit((toString(toString(pattern(p_align)))),""))
     if (length(final_cglist) > 0){
       for (j in 1:length(final_cglist))
       {final_mat[i,j] <- cur_seq[final_cglist[j]] }}}
   t2 <- Sys.time()
   print (round(t2-t1,digits=2))
   return (final_mat)}

#######################################################################
#######################################################################
############Calculating Methylation average for every position########
######################################################################
methAvg <- function(Sample,plot=FALSE)
{ sumt_average <- list()
  for (j in 1:ncol(Sample))
  { sum=length(which(Sample[,j]=="C"))
    sumt=length(which(Sample[,j]=="T"))
    sumt_average[j] <- ((sum)/(sum+sumt))*100   }
  if (plot=="TRUE")
    plot (unlist(sumt_average),type="l",main="Methylation Average",xlab="CpG",ylab="Percentage")
  return (round(unlist(sumt_average),digits=2))  }
########################################################################
###########Function for Counting average and Heatmap####################
########################################################################
methHeatmap <- function(Sample,yl="",plot=TRUE,title="")
{ heat_map <- Sample
  result_percentage <- (length(which(heat_map == "C")) + length(which(heat_map == "T"))) / (nrow(heat_map) * ncol(heat_map)) * 100
  heat_map[heat_map == "A"] <- ""
  heat_map[heat_map == "G"] <- ""
  heat_map[heat_map == "-"] <- ""
  heat_map[heat_map == "C"] <- strtoi(1)
  heat_map[heat_map == "T"] <- strtoi(0)
  suppressWarnings(class(heat_map) <- "numeric")
  csum <- list()
  for (i in 1:nrow(heat_map))
    csum[i] <- length(which(heat_map[i,]==1))
  heat_map <- cbind(heat_map,unlist(csum))
  heat_map2 <- heat_map[order(heat_map[,ncol(heat_map)]),]
  heat_map <- heat_map2[,-ncol(heat_map2)]
  if (plot==TRUE)
    heatmap.2(heat_map,Rowv=NA,Colv=NA,col=c("yellow","blue"),xlab="",ylab=yl,main=title,margins=c(2,2),dendrogram="none",trace="none",key=FALSE,xaxt="n",yaxt="n")
  cat (paste("Percentage Result is",round(result_percentage,digits=2)))
  return (heat_map) }
###########################################################################
#############################END###########################################
###########################################################################

###########################################################################
############################### Entropy ##################################
###########################################################################
methEntropy <- function(x)
{
  t1 <- Sys.time()
  entropy <- list()
  ni_list <- list()
  for (i in 1:(ncol(x)-3))
  {
    sum=0
    entropy_list <- unlist(1:nrow(x))
    for (j in 1:(nrow(x)-1))
    {
      ni_new=0
      for (k in (j+1):length(entropy_list))
      {
        if (identical(x[j,i:(i+3)],x[entropy_list[k],i:(i+3)]))
        { 
          entropy_list <- entropy_list[-k]
          ni_new=ni_new + 1
        }
      }
      
      if (ni_new > 0)
      {
        ni_list[j] <- ni_new
        ni_for_sum <- (-((ni_new/nrow(x))*(log(ni_new/nrow(x)))))
        sum= sum + ni_for_sum
      }
    }
    entropy[i] <- sum/4
  }
  t2 <- Sys.time()
  print (round(t2-t1,digits=2))
  return(round(unlist(entropy),digits=2))
}


#####################Fisher testn#########
fishertest_cpg <- function(healthy,tumor,plot=TRUE,main="Fisher Exact Test")
{
  final_mat <- healthy
  final_mat2 <- tumor
  con_mat <- list()
  for (j in 1:ncol(final_mat))
  {
    a <- length(which(final_mat[,j]=="C"))
    c <- length(which(final_mat[,j]=="T"))
    b <- length(which(final_mat2[,j]=="C"))
    d <- length(which(final_mat2[,j]=="T"))
    M=rbind(c(a, c),c(b,d))
    rownames(M) <- c("C","T")
    colnames(M) <- c("Normal","Tumor")
    con_mat[j] <- fisher.test(M)$p
  }
  con_mat <- -log10(p.adjust(unlist(con_mat),method="fdr"))
  if (plot==TRUE){
    barplot(con_mat,xlab="CpG",ylab="-log10(FDR)",main=main)
    abline(h=-log10(0.05),col="Red",lty=2)}
  return (round(con_mat,digits=2))
}

odd_ratio <- function(SampA,SampB,plot=TRUE,main="Log Odd Ratio")
{
  odd_ratio <- list()
  for (j in 1:ncol(SampA))
  {
    a <- length(which(SampA[,j]=="C"))
    c <- length(which(SampA[,j]=="T"))
    b <- length(which(SampB[,j]=="C"))
    d <- length(which(SampB[,j]=="T"))
    M=rbind(c(a, c),c(b,d))
    rownames(M) <- c("C","T")
    colnames(M) <- c("Normal","Tumor")
    odd_ratio[j]=log((a*d)/(b*c))
  }
  if (plot==TRUE)
    plot(unlist(odd_ratio),xlab="CpG",ylab="log(Odd Ratio)",type="l",main=main)
  return (round(unlist(odd_ratio),digits=2))
}
### compare samples ############

compare_samples <- function(healthy,tumor){
  par(mfrow=c(2,2))
  healthy_average <- methAvg(healthy)
  healthy_entropy <- methEntropy(healthy)
  tumor_average <- methAvg(tumor)
  tumor_entropy <- methEntropy(tumor)
  plot(healthy_average,type="l",col="Green",ylim=c(1,100),main="Methylation Average",ylab="Average in Percentage",xlab="CpGs")
  lines(tumor_average,col="Red")
  plot(healthy_entropy,type="l",col="Green",ylim=c(0,1),main="Methylation Entropy",ylab="Entropy",xlab="CpG")
  lines(tumor_entropy,col="Red")
  ft <- fishertest_cpg(healthy,tumor,plot=TRUE)
  odr <- odd_ratio(healthy,tumor,plot=TRUE)
  par(mfrow=c(1,1))
}
