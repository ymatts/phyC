setwd("/Users/matsui/Documents/lichee/LICHeE/data/")
fname <- list.files("mimoriken_CRC/")
dat <- vector("list",length(fname))
for(i in seq_along(fname)){
  dat[[i]] <- read.table(paste0("mimoriken_CRC/",fname[i]))
  cat(i,"\n")
}

dat2 <- dat
for(i in seq_along(dat)){
  info <- strsplit(rownames(dat[[i]]),split = "::")
  chr <- sapply(info,function(x)x[1])
  pos <- sapply(info,function(x)as.numeric(x[2]))
  des <- sapply(info,function(x)x[5])
  temp <- cbind(chr = chr,pos=pos,description=des,normal=rep(0,nrow(dat[[i]])),dat[[i]])
  dat2[[i]] <- temp[order(temp$pos,decreasing = F),]
  write.table(dat2[[i]],paste0("mimoriken_CRC/s",i,".txt"),quote = F,sep = "\t",row.names = F,col.names = T)  
}

