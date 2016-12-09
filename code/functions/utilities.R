shifter <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}
score <- function(ref,vec,bpUnit,df){
  
  svec <- vec[2:(length(vec)-1)]
  seqL <- max(df[,3])
  L <- length(svec)
  lvec <- vec[1:L]
  rvec <- vec[3:(L+2)]
  ldist <- diag(adist(append(paste0(ref[svec,3],ref[svec,1]),paste0(ref[lvec,2],ref[lvec,4])),ignore.case=TRUE)[1:L,(L+1):(2*L)])
  rdist <- diag(adist(append(paste0(ref[svec,4],ref[svec,2]),paste0(ref[rvec,1],ref[rvec,3])),ignore.case=TRUE)[1:L,(L+1):(2*L)])
  return(sum(ldist+rdist)/(2*sum(1:bpUnit)*(seqL-bpUnit+1)))
}
