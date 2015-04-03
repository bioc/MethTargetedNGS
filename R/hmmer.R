hmmbuild <- function(file_seq,file_out,pathHMMER="")
{
  if (pathHMMER == "")
  {
    cat ("Please install HMMER from http://hmmer.janelia.org/")
  }
  else
    system(paste(pathHMMER,"/hmmbuild --informat afa ",file_out," ",file_seq,sep=""))
}


nhmmer <- function(file_hmm,file_seq,pathHMMER="")
{
  if (pathHMMER == "") 
  {
    cat ("Please install HMMER from http://hmmer.janelia.org/")
    return (NULL)
  }
  else 
    system(paste(pathHMMER,"/nhmmer --tblout hmm_table.out ",file_hmm," ",file_seq,sep=""))
  tab_read <- read.table("hmm_table.out")
  e_value <- cbind(data.frame(tab_read$V1))
  e_value <- cbind(e_value,data.frame(tab_read$V13))
  e_value <- cbind(e_value,data.frame(tab_read$V12))
  e_value <- cbind(e_value,data.frame(tab_read$V14))
  colnames(e_value) <- c("Sequences","E.value","Forward.Reverse","score")
  res <- list(e_value,sum(e_value$score))
  res <- setNames(res,c("HMM.Table.Out","Total.Likelihood.Score"))
  return (res)
}
