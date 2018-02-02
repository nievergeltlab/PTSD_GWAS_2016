#For each dataset, randomly add to pool A or B
#Increment total count for A and B
#If out of balance by more than 5000, add to smaller dataset
#With a simple biased coin design: because UKBB has 30000 people, there will always be roughly 25k in the smaller sample
#Adaptive design not much better until gamma very low!
dat <- read.csv('pgc_ptsd_study_order_v8_eur.csv',stringsAsFactors=F)
replication_count <- 10
diffz=1000

gamma = 0.00005

p_adj = .33333

for (repz in 1:replication_count)
{
 a_size <- 0
 b_size <- 0
 #Null value for assignments
 dat$assignment <- NA

 #For every study, flip a biased coin and assign to a dataset
 for (i in 1:dim(dat)[1])
 {
  #Flip coin 
  #If datasets are approximately equal, set p = 0.5
  
  # Biased coin design
  # if ( abs(a_size - b_size) <= diffz )
  # {
   # p_coin <- .5
  # } else if ( a_size - b_size > diffz )
  # {
   # p_coin <- p_adj
  # } else if ( a_size - b_size < -diffz )
  # {
   # p_coin <- 1-p_adj
  # }
  
  #Bayesian adaptive design (A Atksinon)
    
   probnum = (1 + b_size / ((a_size + b_size)*a_size))^ (1/gamma)
   probdnom1 = (1 + b_size / ((a_size + b_size)*a_size))^ (1/gamma)
   probdnom2 = (1 + a_size / ((a_size + b_size)*b_size))^ (1/gamma)

   p_coin = probnum / (probdnom1 + probdnom2) 
    if(is.na(p_coin))
    {
     p_coin = .5 #Set to default probability if it gets fucked up
     }

  trial <- rbinom(1,1,p=p_coin)
  #print(p_coin)
  #print(trial)
  if (trial == 1)
  {
   dat$assignment[i] <- "A"
   a_size <- a_size + dat$neff[i]
  }
  
  if (trial == 0)
  {
   dat$assignment[i] <- "B"
   b_size <- b_size + dat$neff[i]
  }
 }
 studies_a <- subset(dat,assignment == "A",select=c("descriptor","neff","study"))
 studies_b <- subset(dat,assignment == "B",select=c("descriptor","neff","study"))
 print (sum(studies_a$neff))
 write.csv(studies_a,paste("groupa_",repz,".csv",sep=''),row.names=F)
 write.csv(studies_b,paste("groupb_",repz,".csv",sep=''),row.names=F)

}