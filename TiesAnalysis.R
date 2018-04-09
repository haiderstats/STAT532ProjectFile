######################################################################################################
#This R script was created by Humza Haider on March 8th, 2018.
#Below is the code used to replicate the output of survfit(coxph())$cumhaz in the presence of ties.
######################################################################################################
#References:

#Some survfit source code
#https://github.com/cran/survival/blob/master/R/agsurv.R

#Cagsurv5 code source:
#mtweb.cs.ucl.ac.uk/mus/lib64/R/library/survival/doc/sourcecode.pdf
######################################################################################################

#Libraries
library(survival)
library(MASS)

#Build the survival model
motocox=coxph(Surv(time,cens)~temp, data=motors,ties="efron")

#Indicator vector of deaths/censors
deaths = motors$cens

#Number of deaths per time step (censoring included).
nevent = as.vector(rowsum(deaths,motors$time))

#The risk for each individual.
risk = predict(motocox, type = "risk")

#Reverse cumulative sum function
rcumsum <- function(x) rev(cumsum(rev(x)))

#The cumulative risk moving forward, starting at t = 408
nrisk = rcumsum(rowsum(risk, motors$time))

# This is the variable of interest, that causes a change between lecture and the R output. 
# erisk is the total risk accumulated at each time step. E.g. the risk for t = 408
# is 2*31.323 + 2*1.991 = 66.630. This is essentially allowing us to 'average' the risk and not 
# have to deal with the difference in ordering.
erisk = rowsum(risk*deaths, motors$time)

#A vector of 0s to store the hazards
haz = double(length(nevent))

#We will iterate through each unique event time (either a death or a censor).
#The first if handles censors. In this case nothing happens.
#The second if handles only 1 death, i.e. no ties occured. In this case we divide 1 by the sum of the risk.
#The third if is the event we have a tie. 
#Let d be the number of deaths in a time point and i be the current time point.
#Since we have summed all the risk for each time point (erisk), instead of forcing multiple sums over all the risks, we can take
#1/nrisk[i] + 1/(nrisk[i] - (1/d)* erisk[i]) + 1/(nrisk[i] - (2/d)* erisk[i])  + ... +1/(nrisk[i] - ((d-1)/d)* erisk[i]). 
#This essentially averages across the risks and thus does not have to deal with different values according to the different orders we consider
#the observations.
for (i in 1:length(nevent)) {
  d = nevent[i];
  if(d==0) next
  else if (d==1){
    temp = 1/nrisk[i];
    haz[i] = temp;
  }
  else {
    for (j in 1:d) {
      temp = 1/(nrisk[i] - erisk[i]*(j-1)/d);
      haz[i] = haz[i] + temp;
    }
  }
}

#Our cumulative hazard:
cumsum(haz)
#Showing that this loop is equivalent to the survfit model, to 12 decimal places.
round(cumsum(haz),12) == round(survfit(motocox)$cumhaz,12)



