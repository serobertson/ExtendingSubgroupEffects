

#------------------------------------------------#
#functions to fit the working models for participation, treatment, and outcome
fit_nuisance_functions_parametric<-function(data, Smod, Amod, Ymod){
  S1data<-subset(data, S==1)
  w_reg<-glm(Smod, family="binomial", data=data)
  ps<- predict(w_reg,newdata=data, type="response") 
  w_reg2<-glm(Amod, family="binomial", data=S1data)
  pa<- predict(w_reg2,newdata=data, type="response") 
  
  #weights 
  ipweight= (data$A*data$S )/(ps*pa) + ((1 -data$A) *data$S) /(ps*(1-pa))
  
  ioweight= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))
  
  #Fit outcome model in trial, in each treatment arm
  S1data_A1<-subset(data, S==1 & A==1)
  OM1mod<-glm(Ymod, family="binomial", data=S1data_A1)
  g1<- predict(OM1mod,newdata=data, type="response")
  
  S1data_A0<-subset(data, S==1 & A==0)
  OM0mod<-glm(Ymod, family="binomial", data=S1data_A0)
  g0<- predict(OM0mod,newdata=data, type="response")
  
  #add to dataset
  data$ps<-ps
  data$pa<-pa
  data$ipweight<-ipweight
  data$ioweight<-ioweight
  data$g1<-g1
  data$g0<-g0
  
  list<-list(dat=data, Smod=w_reg, Amod=w_reg2, Ymod1=OM1mod, Ymod0=OM0mod)
  return(list)
}

fit_nuisance_functions_GAM<-function(data, Smod, Amod, Ymod){
  S1data<-subset(data, S==1)
  w_reg<-gam(Smod, family="binomial", data=data)
  ps<- predict(w_reg,newdata=data, type="response") 
  w_reg2<-glm(Amod, family="binomial", data=S1data)
  pa<- predict(w_reg2,newdata=data, type="response") 
  
  #weights 
  ipweight= (data$A*data$S )/(ps*pa) + ((1 -data$A) *data$S) /(ps*(1-pa))

  ioweight= (data$A*data$S*(1-ps) )/(ps*pa) + ((1 -data$A) *data$S*(1-ps) ) /(ps*(1-pa))

  #Fit outcome model in trial, in each treatment arm
  S1data_A1<-subset(data, S==1 & A==1)
  OM1mod<-gam(Ymod, family="binomial", data=S1data_A1)
  g1<- predict(OM1mod,newdata=data, type="response")

  S1data_A0<-subset(data, S==1 & A==0)
  OM0mod<-gam(Ymod, family="binomial", data=S1data_A0)
  g0<- predict(OM0mod,newdata=data, type="response")

  #add to dataset
  data$ps<-ps
  data$pa<-pa
  data$ipweight<-ipweight
  data$ioweight<-ioweight
  data$g1<-g1
  data$g0<-g0
  
  list<-list(dat=data, Smod=w_reg, Amod=w_reg2, Ymod1=OM1mod, Ymod0=OM0mod)

  return(list)
}

#------------------------------------------------#
# estimators for the entire target population
OM_est<-function(data,subgroup){
  
  #calculate in each level of the subgroup
  OM_1<-c()
  OM_0<-c()
  OM_ATE<-c()
  for(v in 1:dim(table(subgroup))){
    S0sub<-subset(data,subgroup==v-1)
    OM_1[v]<-mean(S0sub$g1)
    OM_0[v]<-mean(S0sub$g0)
    OM_ATE[v]<-mean(S0sub$g1)-mean(S0sub$g0)
  }
  
  results=data.frame(cbind(OM_1, OM_0, OM_ATE))
  print(results)
  return(results)
  
}

IPW1_est<-function(data,subgroup){ 
  A<-data$A
  S<-data$S
  weight<-data$ipweight
  Y<-data$Y
  
  #calculate in each level of the subgroup
  IPW1_1<-c()
  IPW1_0<-c()
  IPW1_ATE<-c()
  for(v in 1:dim(table(subgroup))){
    IPW1_1[v] <-(sum((subgroup==v-1))^-1)* sum(A*S*weight*(subgroup==v-1)*Y)
    IPW1_0[v] <-(sum((subgroup==v-1))^-1)* sum((1-A)*S*weight*(subgroup==v-1)*Y)
    IPW1_ATE[v] = IPW1_1[v] - IPW1_0[v]
  }
  
  results=data.frame(cbind(IPW1_1, IPW1_0, IPW1_ATE))
  print(results)
  return(results)
  
}

IPW2_est<-function(data,subgroup){
  
  
  #calculate in each level of the subgroup
  IPW2_1<-c()
  IPW2_0<-c()
  IPW2_ATE<-c()
  
  for(v in 1:dim(table(subgroup))){
    
    data$weight<-data$ipweight*(subgroup==v-1)
    S0data<-subset(data,subgroup==v-1)  
    
    S1data_A1<-subset(data, S==1 & A==1 & subgroup==v-1)
    IPW1mod<-glm(formula=Y~1, data=S1data_A1, weights=weight)
    p1<- predict(IPW1mod,newdata=S0data, type="response") 
    S1data_A0<-subset(data, S==1 & A==0 & subgroup==v-1)
    IPW0mod<-glm(formula=Y~1, data=S1data_A0, weights=weight)
    p0<- predict(IPW0mod,newdata=S0data, type="response") 
    IPW2_1[v]<-mean(p1)
    IPW2_0[v]<-mean(p0)
    IPW2_ATE[v]<-mean(p1)-mean(p0)
  }
  
  results=data.frame(cbind(IPW2_1, IPW2_0, IPW2_ATE))
  print(results)
  return(results)
}

AIPW1_est<-function(data,subgroup){  
  A<-data$A
  S<-data$S
  Y<-data$Y
  g1<-data$g1
  g0<-data$g0
  w<-data$ipweight
  
  #calculate in each level of the subgroup
  AIPW1_1<-c()
  AIPW1_0<-c()
  AIPW1_ATE<-c()
  for(v in 1:dim(table(subgroup))){
    AIPW1_1[v]<-(sum((subgroup==v-1))^-1)* sum(S*A*w*(subgroup==v-1)*(Y-g1) + (subgroup==v-1)*g1)
    AIPW1_0[v]<-(sum((subgroup==v-1))^-1)* sum(S*(1-A)*w*(subgroup==v-1)*(Y-g0) + (subgroup==v-1)*g0)
    AIPW1_ATE[v]<-AIPW1_1[v]-AIPW1_0[v]
  }
  results=data.frame(cbind(AIPW1_1, AIPW1_0, AIPW1_ATE))
  print(results)
  return(results)
}

AIPW2_est<-function(data,subgroup){
  A<-data$A
  S<-data$S
  Y<-data$Y
  g1<-data$g1
  g0<-data$g0
  w<-data$ipweight
  #calculate in each level of the subgroup
  AIPW2_1<-c()
  AIPW2_0<-c()
  AIPW2_ATE<-c()
  for(v in 1:dim(table(subgroup))){
    
    sum1_AIPW2<-sum(S*A*w*(subgroup==v-1)*(Y-g1)) 
    sum0_AIPW2<-sum(S*(1-A)*w*(subgroup==v-1)*(Y-g0)) 
    norm1<-(sum(S*A*w*(subgroup==v-1)))^-1
    norm0<-(sum(S*(1-A)*w*(subgroup==v-1)))^-1
    AIPW2_1[v]<-norm1*sum1_AIPW2 + (sum(subgroup==v-1)^-1)*sum((subgroup==v-1)*g1)
    AIPW2_0[v]<-norm0*sum0_AIPW2 + (sum(subgroup==v-1)^-1)*sum((subgroup==v-1)*g0)
    AIPW2_ATE[v]<-AIPW2_1[v]-AIPW2_0[v]
  }
  
  results=data.frame(cbind(AIPW2_1, AIPW2_0, AIPW2_ATE))
  print(results)
  return(results)
}

#------------------------------------------------#
# estimators for the non-randomized subset
OM_est_S0<-function(data,subgroup){
  
  #calculate the RD in each level of the subgroup
  OM_1<-c()
  OM_0<-c()
  OM_ATE<-c()
  for(v in 1:dim(table(subgroup))){
    S0sub<-subset(data, S==0 & subgroup==v-1)
    OM_1[v]<-mean(S0sub$g1)
    OM_0[v]<-mean(S0sub$g0)
    OM_ATE[v]<-mean(S0sub$g1)-mean(S0sub$g0)
  }
  
  results=data.frame(cbind(OM_1, OM_0, OM_ATE))
  print(results)
  return(results)
  
}

IOW1_est_S0<-function(data,subgroup){ 
  A<-data$A
  S<-data$S
  weight<-data$ioweight
  Y<-data$Y
  
  #calculate in each level of the subgroup
  IOW1_1<-c()
  IOW1_0<-c()
  IOW1_ATE<-c()
  for(v in 1:dim(table(subgroup))){
  IOW1_1[v] <-(sum((subgroup==v-1 & S==0))^-1)* sum(A*S*weight*(subgroup==v-1)*Y)
  IOW1_0[v] <-(sum((subgroup==v-1 & S==0))^-1)* sum((1-A)*S*weight*(subgroup==v-1)*Y)
  IOW1_ATE[v] = IOW1_1[v] - IOW1_0[v]
  }
  
  results=data.frame(cbind(IOW1_1, IOW1_0, IOW1_ATE))
  print(results)
  return(results)
  
  }

IOW2_est_S0<-function(data,subgroup){
  
  
  #calculate in each level of the subgroup
  IOW2_1<-c()
  IOW2_0<-c()
  IOW2_ATE<-c()
  
  for(v in 1:dim(table(subgroup))){
  
  data$weight<-data$ioweight*(subgroup==v-1)
  S0data<-subset(data, S==0 & subgroup==v-1)  
  
  S1data_A1<-subset(data, S==1 & A==1 & subgroup==v-1)
  IOW1mod<-glm(formula=Y~1, data=S1data_A1, weights=weight)
  p1<- predict(IOW1mod,newdata=S0data, type="response") 
  S1data_A0<-subset(data, S==1 & A==0 & subgroup==v-1)
  IOW0mod<-glm(formula=Y~1, data=S1data_A0, weights=weight)
  p0<- predict(IOW0mod,newdata=S0data, type="response") 
  IOW2_1[v]<-mean(p1)
  IOW2_0[v]<-mean(p0)
  IOW2_ATE[v]<-mean(p1)-mean(p0)
  }
  
  results=data.frame(cbind(IOW2_1, IOW2_0, IOW2_ATE))
  print(results)
  return(results)
}

AIOW1_est_S0<-function(data,subgroup){  
  A<-data$A
  S<-data$S
  Y<-data$Y
  g1<-data$g1
  g0<-data$g0
  w<-data$ioweight
  
  #calculate in each level of the subgroup
  AIOW1_1<-c()
  AIOW1_0<-c()
  AIOW1_ATE<-c()
  for(v in 1:dim(table(subgroup))){
  AIOW1_1[v]<-(sum((subgroup==v-1 & S==0))^-1)* sum(S*A*w*(subgroup==v-1)*(Y-g1) + (subgroup==v-1 & S==0)*g1)
  AIOW1_0[v]<-(sum((subgroup==v-1 & S==0))^-1)* sum(S*(1-A)*w*(subgroup==v-1)*(Y-g0) + (subgroup==v-1 & S==0)*g0)
  AIOW1_ATE[v]<-AIOW1_1[v]-AIOW1_0[v]
  }
  results=data.frame(cbind(AIOW1_1, AIOW1_0, AIOW1_ATE))
  print(results)
  return(results)
}

AIOW2_est_S0<-function(data,subgroup){
  A<-data$A
  S<-data$S
  Y<-data$Y
  g1<-data$g1
  g0<-data$g0
  w<-data$ioweight
  #calculate in each level of the subgroup
  AIOW2_1<-c()
  AIOW2_0<-c()
  AIOW2_ATE<-c()
  for(v in 1:dim(table(subgroup))){
    
  sum1_AIOW2<-sum(S*A*w*(subgroup==v-1)*(Y-g1)) 
  sum0_AIOW2<-sum(S*(1-A)*w*(subgroup==v-1)*(Y-g0)) 
  norm1<-(sum(S*A*w*(subgroup==v-1)))^-1
  norm0<-(sum(S*(1-A)*w*(subgroup==v-1)))^-1
  AIOW2_1[v]<-norm1*sum1_AIOW2 + (sum(subgroup==v-1 & S==0)^-1)*sum((subgroup==v-1 & S==0)*g1)
  AIOW2_0[v]<-norm0*sum0_AIOW2 + (sum(subgroup==v-1 & S==0)^-1)*sum((subgroup==v-1 & S==0)*g0)
  AIOW2_ATE[v]<-AIOW2_1[v]-AIOW2_0[v]
  }
  
  results=data.frame(cbind(AIOW2_1, AIOW2_0, AIOW2_ATE))
  print(results)
  return(results)
}
