#' Internal function
#' 
#' @name lichee2edge
#' @docType package
#' @import igraph,ape
#' @author Yusuke Matsui & Teppei Shimamura
#' @export
#' 
lichee2edge <- function(input=NULL,licheeParamIO=NULL,licheeParamFilter=NULL,licheeParamPhy=NULL,outdir=NULL,idx=1){
  options(scipen = 9)
  if(is.null(unlist(licheeParamIO))){
    normal <- 0
    nsave <- 1
    net <- FALSE
    ntree <- 0
  }else{
    normal <- licheeParamIO$normal;if(is.null(normal)){normal <- 0}
    nsave <- licheeParamIO$save;if(is.null(nsave)){nsave <- 1}
    net <- licheeParamIO$net;if(is.null(net)){net <- FALSE}
    ntree <- licheeParamIO$tree;if(is.null(ntree)){ntree <- 0}
  }
  
  if(is.null(unlist(licheeParamFilter))){
    present <- 0.05
    absent <- 0.05
    maxVAFValid <- 0.6
    minProfileSupport <- 2
  }else{
    present <- licheeParamFilter$present;if(is.null(present)){present <- 0.05}
    absent <- licheeParamFilter$absent;if(is.null(absent)){absent <- 0.05}
    maxVAFValid <- licheeParamFilter$maxVAFValid;if(is.null(maxVAFValid)){maxVAFValid <- 0.6}
    minProfileSupport <- licheeParamFilter$minProfileSupport;if(is.null(minProfileSupport)){minProfileSupport <- 2}
  }
  
  if(is.null(unlist(licheeParamPhy))){
    minClusterSize <- 2
    minPrivateClusterSize <- 1
    minRobustNodeSupport <- 2
    maxClusterDist <- 0.2
    completeNetwork <- FALSE
    e <- 0.1
    nTreeQPCheck <- 0
  }else{
    minClusterSize <- licheeParamPhy$minClusterSize;if(is.null( minClusterSize)){minClusterSize <- 2}
    minPrivateClusterSize <- licheeParamPhy$minPrivateClusterSize;if(is.null(minProfileSupport)){minPrivateClusterSize <- 1}
    minRobustNodeSupport <- licheeParamPhy$minRobustNodeSupport;if(is.null(minRobustNodeSupport)){minRobustNodeSupport <- 2}
    maxClusterDist <- licheeParamPhy$maxClusterDist;if(is.null(maxClusterDist)){maxClusterDist <- 0.2}
    completeNetwork <- licheeParamPhy$completeNetwork;if(is.null(completeNetwork)){completeNetwork <- FALSE}
    e <- licheeParamPhy$e;if(is.null(e)){e <- 0.1}
    nTreeQPCheck <- licheeParamPhy$nTreeQPCheck;if(is.null(nTreeQPCheck)){nTreeQPCheck <- 0}
  }
  
  if(is.null(outdir)){
    if(!dir.exists("out")){dir.create("out");outdir <- "out"}
  }else{
    if(!dir.exists(outdir)){dir.create(outdir)}
  }
  #output <- paste0(outdir,"/o",idx,"_p",present,"_a",absent,"_mv",maxVAFValid,"_mps",minProfileSupport,"_mcs",minClusterSize,"_mpcs",minPrivateClusterSize,"_mrns",minRobustNodeSupport,"_mcd",maxClusterDist,"_e",e,"_m.txt")
  output <- paste0(outdir,"/o",idx,"_p",present,"_a",absent,"-minClusterSize",minClusterSize,".txt")
  
  cmd <- paste("java -jar lichee.jar -build -i", input, "-o","tmp.txt","-minVAFPresent", present, "-maxVAFAbsent",absent,"-n 0", "-s",nsave,"-showTree",ntree,"-minClusterSize",minClusterSize,"-minPrivateClusterSize", minPrivateClusterSize,"-minRobustNodeSupport",minRobustNodeSupport,"-maxClusterDist",maxClusterDist,"-e",e,"-nTreeQPCheck",nTreeQPCheck)  
  
  if(net){paste(cmd,"-showNetwork")}
  if(completeNetwork){paste(cmd,"-c")}
  
  system(cmd,intern = F)
  
  #cmd <- paste("grep","\'>\'","tmp.txt")
  cmd <- paste("grep","\\->","tmp.txt")
  temp <- system(cmd,intern=T)
  temp <- sapply(temp,function(x)gsub(" -> ","\t",x))
  temp <- sapply(temp,function(x)strsplit(x,"\t"))
  edgelist <- t(sapply(temp,function(x)as.numeric(x))) + 1
  
  maxind <- max(edgelist)
  cmd <- paste("head","-n",maxind,"tmp.txt")
  temp <- system(cmd,intern=T)[-1]
  temp <- sapply(temp,function(x)unlist(strsplit(x,"]"))[2])
  temp <- sapply(temp,function(x)unlist(strsplit(x,"\t"))[-1])
  len <- sapply(temp,length)
  edge <- cbind(edgelist,len)
  colnames(edge) <- c("v1","v2","length")
  write.table(edge,file=output,row.names=F,col.names=T,quote=F,sep="\t")
  
}
