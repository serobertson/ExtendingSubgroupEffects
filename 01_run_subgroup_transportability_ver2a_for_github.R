
#set your working directory here
source('00_subgroup_generalizability_source_code_ver2_for_github.R')

#create a simulated dataset to illustrate the methods
n=2000
DF_sim = data.frame( prevmi=rbinom(n, size=1,prob=0.50),
                     age=runif(n),angina=rbinom(n, size=1,prob=0.50),ejecfr=runif(n),lvscor=runif(n),
                     prxlad31=runif(n),prxves31=runif(n))
DF_sim$A<-rbinom(n, 1, 0.5)
DF_sim$S <- rbinom(n, 1, plogis(-1 + 0.5*( DF_sim$prevmi + DF_sim$age + DF_sim$angina + DF_sim$ejecfr)))
DF_sim$Y <- rbinom(n, 1, plogis(-1 + 0.5*( DF_sim$prevmi + DF_sim$age + DF_sim$angina + DF_sim$ejecfr)*DF_sim$A))

head(DF_sim)


run_all<-function(DF=DFcass,bootstrap_sample=bootstrap_sample, modeling="parametric"){
  
  if(bootstrap_sample==F){DFbs<-DF}
  else if(bootstrap_sample==T) {
      index <- sample(1:nrow(DF), replace=T)
      DFbs<-DF[index, ]
    }
  
  if(modeling=="parametric"){

     model=fit_nuisance_functions_parametric(data=DFbs, Smod=S~age+angina+ejecfr+lvscor+prevmi+prxlad31+prxves31,
                                             Amod=A~age+angina+ejecfr+lvscor+prevmi+prxlad31+prxves31,
                                             Ymod=Y~age+angina+ejecfr+lvscor+prevmi+prxlad31+prxves31)
  }
  
  if(modeling=="gam"){
    library(mgcv) 
    model=fit_nuisance_functions_GAM(data=DFbs, Smod=S~s(age)+angina+s(ejecfr)+lvscor+prevmi+prxlad31+prxves31,
                                     Amod=A~age+angina+ejecfr+lvscor+prevmi+prxlad31+prxves31,
                                     Ymod=Y~s(age)+angina+s(ejecfr)+lvscor+prevmi+prxlad31+prxves31)
  }
 
DFfit=model$dat

#fit among the entire target population
OM=OM_est(data=DFfit, subgroup=DFfit$prevmi)
IPW1=IPW1_est(data=DFfit, subgroup=DFfit$prevmi)
IPW2=IPW2_est(data=DFfit, subgroup=DFfit$prevmi)
AIPW1=AIPW1_est(data=DFfit, subgroup=DFfit$prevmi)
AIPW2=AIPW2_est(data=DFfit, subgroup=DFfit$prevmi)

all_estimates<-data.frame(cbind(OM, IPW1, IPW2, AIPW1, AIPW2))
tall_estimates<-t(all_estimates)


#fit among S=0
OM_S0=OM_est_S0(data=DFfit, subgroup=DFfit$prevmi)
IOW1=IOW1_est_S0(data=DFfit, subgroup=DFfit$prevmi)
IOW2=IOW2_est_S0(data=DFfit, subgroup=DFfit$prevmi)
AIOW1=AIOW1_est_S0(data=DFfit, subgroup=DFfit$prevmi)
AIOW2=AIOW2_est_S0(data=DFfit, subgroup=DFfit$prevmi)

all_estimates_S0<-data.frame(cbind(OM_S0, IOW1, IOW2, AIOW1, AIOW2))
tall_estimates_S0<-t(all_estimates_S0)

subgroup=DFfit$prevmi
return(data.frame(cbind(subgroup=c(1:dim(table(subgroup))), all_estimates, all_estimates_S0)))

}

set.seed(123)
result_original_parametric=run_all(DF=DFcass, bootstrap_sample=F, modeling="parametric")
result_original_gam=run_all(DF=DFcass, bootstrap_sample=F, modeling="gam")

set.seed(123)
result_onesample=run_all(DF=DFcass, bootstrap_sample=T, modeling="parametric")


#nReps=10000
nReps=10
set.seed(123)
sim<-do.call(rbind, replicate(nReps, as.matrix(run_all(DF=DFcass, bootstrap_sample=T, modeling="parametric")), simplify=FALSE)) 
simresults_parametric<-as.data.frame(sim)
 
file.name<- paste('simresults_parametric', '_nReps', nReps,'_v', format(Sys.time(), "%d%b%Y"), '.Rdata', sep='')
save(simresults_parametric, file=file.name)

#nReps=10000
nReps=10
set.seed(123)
sim<-do.call(rbind, replicate(nReps, as.matrix(run_all(DF=DFcass, bootstrap_sample=T, modeling="gam")), simplify=FALSE)) 
simresults_gam<-as.data.frame(sim)

file.name<- paste('simresults_gam', '_nReps', nReps,'_v', format(Sys.time(), "%d%b%Y"), '.Rdata', sep='')
save(simresults_gam, file=file.name)


