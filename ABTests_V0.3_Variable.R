rm(list=ls())
library(exact2x2)

AB_fixedN <- function(out_df, n, pB, log_dif){
  #n = trunc(runif(1, 100, 10000))
  #pB = runif(1, 0.0001, 0.03)
  #log_dif = runif(1, -0.18, 0.18)
  pA <- exp(log_dif + log(pB))
  #print(log(pA/pB))
  nA <- rbinom(1, n, .5)
  nB <- n-nA
  
  sA = pA * nA
  sB = pB * nB
  fA = nA - sA
  fB = nB - sB
  
  data <- matrix(c(sA,fA,sB,fB), nrow = 2, byrow = TRUE, dimnames = list(Row= c("A","B"), Col = c("Success","Failure")))
  
  #print(data)
  
  #Pearson Chi Test
  pearson_chi <- chisq.test(data, correct = F)$p.value
  #print(pearson_chi)
  
  #N-1 chi squared Tes
  n1_chi_stat <- chisq.test(data, correct = F)$statistic * ((n-1)/n)
  n1_chi <- pchisq(n1_chi_stat, df=1, lower.tail=FALSE)[["X-squared"]]
  #print(n1_chi)
  
  #Yates Chi Test
  yates_chi <- chisq.test(data, correct = T)$p.value
  #print(yates_chi)
  
  #Fisher Test
  ft <- fisher.test(data, alternative = 'greater', conf.int = TRUE)$p.value
  #print(ft)
  
  
  
  ft2 <- fisher.exact(data,tsmethod="minlike")$p.value
  #print(ft2)
  
  blaker_val <- blaker.exact(data)$p.value
  #print(blaker_val)
  
  ex_val <- exact2x2(data)$p.value
  #print(ex_val)
  
  #unc_val <- uncondExact2x2(sA, n, sB, n)$p.value
  ##print(unc_val)
  
  Howell_OneMargin <- function(origdata,tot){
    
    rowmarg <- rowSums(origdata)
    nrow <- length(rowmarg)
    colmarg <- colSums(origdata)
    ncol <- length(colmarg)
    colprob <- colmarg/tot
    
    obtchisq <- chisq.test(origdata, correct = F)$statistic
    
    # Now start the resampling 
    nreps <- 1000
    chisquare <- numeric(nreps)
    results <- numeric(nreps)
    extreme <- 0
    
    cum <- numeric()
    cum[1] <- colprob[1]
    for (j in 2:ncol) {cum[j] <- cum[j-1] + colprob[j]}
    
    # Now resample nreps = 10,000 times 
    count <- c()
    countx <- c()
    for (i in 1:nreps) {
      
      randomtable <- c()
      for (j in 1:nrow) {
        randnum <- runif(rowmarg[j],0,1) 
        for (j in 1:ncol) {count[j] <- length(randnum[randnum <= cum[j]])}
        countx[1] <- count[1]
        for (k in 2:ncol) {countx[k] <- count[k] - count[k-1]} 
        randomtable <- rbind(randomtable, countx)
      }
      chisquare[i] <- chisq.test(randomtable, correct = FALSE)$statistic
      
    }
    signif <- length(chisquare[chisquare >= obtchisq]) /nreps
    
    return(signif)
  }
  
  how_1M <- Howell_OneMargin(data,2*n)
  #print(how_1M)
  
  Howell_NoMargin <- function(origdata) {
    N <- sum(origdata)
    nr <- nrow(origdata)
    nc <- ncol(origdata)
    ncells <- nr*nc
    # Establish row and column marginal totals
    rowtots <- rowSums(origdata)
    coltots <- colSums(origdata)
    # Calculate expected cell probabilities based on original marginal totals
    # and convert to vector
    expected <- matrix(data = NA, nrow = nr, ncol = nc, byrow = TRUE)
    for (i in 1:nr) {
      for (j in 1:nc) {
        expected[i,j] <- rowtots[i] * coltots[j] /N^2
      }
    }
    # Convert to vector  
    expectedvector <- as.vector(t(expected))
    # This can now be sampled
    
    # Calculate chi-square on original data
    obtchisq <- chisq.test(origdata, correct = F)$statistic    
    
    nreps <- 1000                    
    xx <- numeric(N)
    results <- numeric(nreps)
    greater <- 0                
    cells <- numeric(ncells)
    cells <- c(1:ncells) 
    # Now resample nreps times
    
    for (i in 1:nreps) {
      xx <- sample(cells,N,replace = TRUE, prob = expectedvector)
      vv <- matrix(0,1,ncells)
      for (j in 1:N) {
        b <- xx[j]
        vv[b] = vv[b] + 1
      }
      x <- matrix(vv,nrow = nr, byrow = TRUE)
      # x is now the contingency table for ith random sample
      
      results[i]<- chisq.test(x, correct = FALSE)$statistic
      greater <- ifelse(results[i] >= obtchisq, greater + 1, greater + 0)  
    }
    signif <- length(results[results >= obtchisq]) /nreps
    
    return(signif)
  }
  
  how_NM <- Howell_NoMargin(data)
  #print(how_NM)
  
  out_df <- rbind(out_df,data.frame(n,nA, nB, pA, pB, pearson_chi, n1_chi, yates_chi, ex_val, ft2, blaker_val, how_1M, how_NM))
  
  return(out_df)
}

output <- data.frame(N = double(),nB = double(), nA = double(), pA = double(), pB = double(),Pearson_chi_sq = double(), n1_chi_sq = double(), yates = double(),exact = double(), fisher = double(), blaker = double(), howell_onemargin = double(), howell_nomargin = double(), stringsAsFactors = TRUE)
#new_row = list(n,pA, pB, pearson_chi, n1_chi, yates_chi, ex_val, ft2, blaker_val, ex_val, how_1M, how_NM)

#out_df <- rbind(out_df2,data.frame(n,pA, pB, pearson_chi, n1_chi, yates_chi, ex_val, ft2, blaker_val, ex_val, how_1M, how_NM))

#size = c(100, 1000, 10000, 100000)
#prob = c(0.0001, 0.0005, 0.001, 0.005, 0.03)
#diff = c(0, 0.03, 0.06, 0.09, 0.18)

size = c(200, 2000,20000,200000)
prob = c(0.0001, 0.0005, 0.001, 0.03 )
diff = c(0, 0.03, 0.09, 0.18)

for (i in 1:4) {
  for (j in 1:4) {
    for (k in 1:4){
      print( (i-1)*16 + (j-1)*4 + (k-1))
      output <- AB_fixedN(output,size[i],prob[j],diff[k])
    }
  }
}

write.csv(output,"Output_VariableN.csv")



