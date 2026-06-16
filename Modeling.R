
library(pheatmap)
# Parameters

# Define a function for generations
gen <- function(p,q,gens,Ne,z = 0){
  for(i in 1:gens){
    if (p == 1 || p == 0){
      break
    }
    random.matrix <- matrix(runif(Ne,0,1), nrow= 2*Ne, ncol=1) 
    allele.m <- random.matrix >= q
    psum<-sum(allele.m==TRUE)
    qsum<-sum(allele.m==FALSE)
    
    p <- psum/(2*Ne)
    q <- qsum/(2*Ne)
    
    PP <- p*p
    PQ <- 2*p*q
    QQ <- q*q*(1-z)
    S <- PP+PQ+QQ
    
    p <- PP/S+PQ/S/2
    q <- 1-p
  }
  
  return(c(p,q))
}


NeSeq <- c(1e02)
ZSeq <- c(0,0.05,0.1,0.15,0.2)
m <- matrix(0, nrow = length(NeSeq), ncol = length(ZSeq))
rownames(m) <- NeSeq
colnames(m) <- ZSeq
for(i in 1:length(NeSeq)){
  for(j in 1:length(ZSeq)){
    runs <-1000
    Ne <- NeSeq[i]
    z <- ZSeq[j]
    test <- matrix(0, nrow=2, ncol=runs)
    
    initial.pfreq<- 0.1
    initial.qfreq <- 1 - initial.pfreq
    generations <- 1000		#Number of gens to cycle through for allele frequencies
    for(k in seq(1,runs,1)){
      test[,k] <- gen(initial.pfreq, initial.qfreq,generations,Ne, z)
    }
    m[i,j] <- sum(test[1,]==1)/runs
    print(m)
  }
}
print(m)

pheatmap(m, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, breaks = seq(0,1,0.01))

