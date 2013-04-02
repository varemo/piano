print.GSC <- function(x, ...){
   n <- 50
   if(length(x$gsc) < 50) n <- length(x$gsc)
   m <- 50
   if(length(unique(unlist(x$gsc))) < 50) m <- length(unique(unlist(x$gsc)))
   cat(paste("First ",n," (out of ",length(x$gsc),") gene set names:\n",sep=""))
   tmp <- names(x$gsc)[1:n]
   for(i in 1:n) {
      if(nchar(tmp[i])>15) tmp[i] <- paste(substr(tmp[i],1,15),"...",sep="")
   }
   print(tmp)
   cat(paste("\nFirst ",m," (out of ",length(unique(unlist(x$gsc))),") gene names:\n",sep=""))
   print(unique(unlist(x$gsc))[1:m])
   cat("\nGene set size summary:\n")
   print(summary(unlist(lapply(x$gsc,length))))
   if(class(x$addInfo) == "data.frame") {
      k <- 10
      if(nrow(x$addInfo) < 10) k <- nrow(x$addInfo)
      cat(paste("\nFirst",k,"gene sets with additional info:\n"))
      tmp <- x$addInfo[1:k,] 
      colnames(tmp) <- c("Gene set","Additional info")
      rownames(tmp) <- NULL
      for(i in 1:k) {
         if(nchar(tmp[i,1])>20) tmp[i,1] <- paste(substr(tmp[i,1],1,20),"...",sep="")
         if(nchar(tmp[i,2])>60) tmp[i,2] <- paste(substr(tmp[i,2],1,60),"...",sep="")
      }
      print(tmp)
   } else {
      cat("\nNo additional info available.\n")  
   }
}