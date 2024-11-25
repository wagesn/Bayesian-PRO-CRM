###install required R packages
#install.packages("nnet")
#install.packages("binom")
#install.packages("dfcrm")
library(nnet)
library(binom)
library(dfcrm)

###Load the function 'procrm' 
procrm<-function(truthc,truthp,skeletonc,skeletonp,sc,sp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl){
 
 bcrmh<-function(a,p,y,n,s2){
 lik=exp(-0.5*a*a/s2)
 for(j in 1:length(p)){
 pj=p[j]**exp(a)
 lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
 }
 return(lik);
     }
 
 bcrmht<-function(a,p,y,n,s2){
 lik=a*exp(-0.5*a*a/s2)
 for(j in 1:length(p)){
 pj=p[j]**exp(a)
 lik=lik*pj^y[j]*(1-pj)^(n[j]-y[j]);
 }
 return(lik);
     }

    ###simulate a trial 	
    ndose = length(skeletonc);   #number of combos
    yc=yp=n=rep(0,ndose);  #number of toxicities and number of treated patients at each dose level
    curr = start;  # current dose level	 
    ptox.hatc = ptox.hatp = numeric(ndose); # estimate of toxicity prob
    dose.select=rep(0,ndose); # a vector of indicators for dose selection
    stopc=stopp=0; #indicate if trial stops early
    i=1	
while(i <= ncohort)
    {
	# generate data for a new cohort of patients
		yc[curr] = yc[curr] + rbinom(1,cohortsize,truthc[curr]);
		yp[curr] = yp[curr] + rbinom(1,cohortsize,truthp[curr]);
		n[curr] = n[curr] + cohortsize;

		if(any(n>n.stop)){
			stop<-0
			break
		}
 	##Model-based estimation of DLT probabilities
 	marginalc = integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletonc, y=yc,n=n, s2=sc,abs.tol = 0)$value;
 	estc=integrate(bcrmht,lower=-10,upper=10, skeletonc, yc, n, sc,abs.tol = 0)$value/marginalc

 	marginalp = integrate(bcrmh,lower=-Inf,upper=Inf, p=skeletonp, y=yp,n=n,s2=sp,abs.tol = 0)$value;
 	estp=integrate(bcrmht,lower=-10,upper=10, skeletonp, yp, n, sp, abs.tol = 0)$value/marginalp

 	ptox.hatc=skeletonc**exp(estc)
 	ptox.hatp=skeletonp**exp(estp)

     #########stopping rules
     safetyc=binom.confint(yc[1],n[1],conf.level=cl,methods="agresti-coull")$lower
     if(safetyc>targetc){
      stopc=1
      break
     }
    
      safetyp=binom.confint(yp[1],n[1],conf.level=cl,methods="agresti-coull")$lower
      if(safetyp>targetp){
        stopp=1
        break
     }
     
  	##Allocation algorithm
 	distancec=abs(ptox.hatc-targetc)
 	distancep=abs(ptox.hatp-targetp)

 	bestc=which.is.max(-distancec)
 	bestp=which.is.max(-distancep)
 	best = min(bestc,bestp)
 	curr=min(best,curr+1)

	i=i+1
	}
	if(stopc==0 & stopp==0){
		dose.select[curr]=dose.select[curr]+1;
		}
	return(list(dose.select=dose.select,tox.datac=yc,tox.datap=yp,pt.allocation=n,stopc=stopc,stopp=stopp))
}
##########'procrm' end here

###Load the function 'procrm.sim' 
procrm.sim<-function(truthc,truthp,skeletonc,skeletonp,sc,sp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl,ntrial){
	ndose=length(truthc)
	
	dose.select<-yc<-yp<-n<-matrix(nrow=ntrial,ncol=ndose)
	nstopc=nstopp=0
	
	for(i in 1:ntrial){
		result<-procrm(truthc,truthp,skeletonc,skeletonp,sc,sp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl)
		dose.select[i,]=result$dose.select
		yc[i,]=result$tox.datac
		yp[i,]=result$tox.datap
		n[i,]=result$pt.allocation
		nstopc=nstopc+result$stopc
		nstopp=nstopp+result$stopp
	}

 cat("Simulation results for the Bayesian PRO-CRM design (Lee et al, 2020)\n");
 cat("targeting a NCI-CTCAE DLT rate of", targetc,"\n");
 cat("and targeting a PRO-CTCAE DLT rate of", targetp,"\n\n");

	cat("True NCI-CTCAE DLT probability:\n");
	cat(round(truthc,3), sep="\t",  "\n");
	cat("True PRO-CTACE DLT probability:\n");
      cat(round(truthp,3), sep="\t",  "\n");
	cat("MTD selection percentage:\n");
	cat(formatC(colMeans(dose.select)*100, digits=1, format="f"), sep="\t",  "\n");
	cat("Average number of NCI-CTCAE DLTs:\n");
      cat(formatC(colMeans(yc), digits=1, format="f"), sep="\t",   "\n");
	cat("Average number of PRO-CTCAE DLTs:\n");
	cat(formatC(colMeans(yp), digits=1, format="f"), sep="\t",   "\n");
	cat("Average number of patients treated:\n");
	cat(formatC(colMeans(n), digits=1, format="f"), sep="\t",   "\n");
	cat("Percentage of trials stopped for NCI-CTCAE safety:\n");
	cat(nstopc/ntrial*100, "\n");
	cat("Percentage of trials stopped for PRO-CTCAE safety:\n");
	cat(nstopp/ntrial*100, "\n");
}
##########'procrm.sim' end here



####################################################
#
# Phase I study evaluating adjuvant hypofractionated
# whole pelvis radiation therapy (WPRT) in
# endometrial cancer (NCT04458402)
#
####################################################
start=1  		##starting dose
targetc=0.20      ##target c toxicity rate 
targetp=0.55      ##target p toxicity rate 

##Specify a set of skeleton values
skeletonc<-c(0.20,0.30)
skeletonp<-c(0.55,0.65)

##specify the prior variance for each model
sc <- 1.6**2
sp <- 1.58**2

##True DLT probability scenarios
c1<-c(0.05,0.15)
p1<-c(0.18,0.35)

c2<-c(0.20,0.40)
p2<-c(0.18,0.35)

c3<-c(0.10,0.20)
p3<-c(0.35,0.55)

c4<-c(0.08,0.15)
p4<-c(0.50,0.65)

c5<-c(0.08,0.15)
p5<-c(0.65,0.75)

c6<-c(0.40,0.45)
p6<-c(0.25,0.35)

cohortsize=3    	##cohort size for each inclusion
ncohort=5      	##number of cohorts
n.stop=25       	##Number of patients needed on one combination to stop the trial
ntrial=10000      ##number of simulated trials 
cl=0.7       	##confidence level for the confidence interval 

truthc<-c5
truthp<-p5
set.seed(34895)
procrm.sim(truthc,truthp,skeletonc,skeletonp,sc,sp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl,ntrial)



##############################################
#
# Comparison to Lee et al (Stat Med, 2020)
#
##############################################
start=1  		##starting dose
targetc=0.25      ##target c toxicity rate 
targetp=0.35      ##target p toxicity rate 

###Specify a set of skeleton values
skeletonc<-c(0.08,0.16,0.25,0.35,0.46)
skeletonp<-c(0.16,0.25,0.35,0.45,0.54)

##specify the prior variance for each model
sc <- 0.79**2
sp <- 0.74**2

###True DLT probability scenarios
c1<-c(0.05,0.05,0.25,0.40,0.55)
p1<-c(0.17,0.18,0.35,0.50,0.65)

c2<-c(0.05,0.25,0.40,0.55,0.70)
p2<-c(0.10,0.15,0.35,0.50,0.65)

c3<-c(0.01,0.02,0.05,0.10,0.25)
p3<-c(0.04,0.09,0.17,0.20,0.35)

c4<-c(0.02,0.05,0.10,0.25,0.40)
p4<-c(0.09,0.17,0.20,0.35,0.50)

c5<-c(0.05,0.10,0.16,0.25,0.40)
p5<-c(0.05,0.20,0.35,0.50,0.65)

c6<-c(0.05,0.18,0.20,0.25,0.40)
p6<-c(0.17,0.35,0.50,0.65,0.80)

c7<-c(0.01,0.05,0.10,0.16,0.25)
p7<-c(0.04,0.05,0.20,0.35,0.50)

cohortsize=1 	##cohort size for each inclusion
ncohort=18   	##number of cohorts
n.stop=25    	##Number of patients needed on one combination to stop the trial
ntrial=10000   	##number of simulated trials 
cl=0.9999         ##confidence level for the confidence interval 

truthc<-c7
truthp<-p7
set.seed(234667)   ##random seed
procrm.sim(truthc,truthp,skeletonc,skeletonp,sc,sp,targetc,targetp,cohortsize,ncohort,n.stop,start,cl,ntrial)
