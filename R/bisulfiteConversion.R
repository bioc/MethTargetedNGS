################bisulfite conversion##############
bconv <- function(fasta_file,out_file="output.fasta"){
  ref_seq <- read.fasta(fasta_file,as.string=TRUE,set.attributes=FALSE)
  for (i in 1:length(ref_seq))
  {
    final_seq <- ""
    sub_seq <- toString(ref_seq[i])
    cur_seq <- unlist(strsplit((toString(sub_seq)),""))
    for (j in 1:(length(cur_seq)-1))
    {
      if ((cur_seq[j]=="c") && (!(cur_seq[j+1]=="g")))
        final_seq <- paste(final_seq,"t",sep="")
      else
        final_seq <- paste(final_seq,cur_seq[j],sep="")
      
    }
    if (out_file!="output.fasta"){
      out_file = fasta_file
      write.fasta(final_seq,names(ref_seq)[i],file.out=paste(out_file,"_bisulfiteconverted.fasta"),open="a")
    }
  }
}
