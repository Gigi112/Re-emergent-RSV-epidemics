I_ex_model <- function(t,y,parms, time.step='month'){
  
  States<-array(y, dim=dim(parms$yinit.matrix))
  dimnames(States) <- dimnames(parms$yinit.matrix)
  
  if(parms$time.step=='month'){
    period=12
    length.step=30.44 #days
  }else if(parms$time.step=='week'){
    period=52.1775
    length.step=7 #days
  }
  
  
  if(t>=466 & t<= 484){
    I_ex= sum(States)*parms$ex_proportion*ex_increase_new[t-465] # ex_increase_new is a vector indicating how 
    # travel restritions and virus activities in other regions affects virus importation
  }else{
    I_ex <- sum(States)*parms$ex_proportion 
    ## external infections
    ##  I_ex=1/100000*St ; I_ex=5/100000*St ; I_ex=10/100000*St and etc
  }
  
  omega = 1/(parms$DurationMatImmunityDays/length.step)
  
  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(WidthAgeClassMonth*4.345)
  }
  
  gamma1= 1/(parms$dur.days1/length.step)  #converts 1/days to 1/lenth.step
  gamma2= 1/(parms$dur.days2/length.step)  
  gamma3= 1/(parms$dur.days3/length.step)  
  gamma4= gamma3  #??Is this right?? Yes, gamma3 stands for rate of recovery from subsequent infection
  
  #Pull out the states  for the model as vectors
  M <-  States[,'M']
  S0 <-  States[,'S0']
  I1 <-  States[,'I1']
  
  S1 <-  States[,'S1']
  I2 <-  States[,'I2']
  
  S2 <-  States[,'S2']
  I3 <-  States[,'I3']
  
  S3 <-  States[,'S3']
  I4 <-  States[,'I4']
  
  N.ages <- length(M)
  
  ###Check the standardization of beta and overall structure of lambda here
  ##force of transmission##################
  seasonal.txn <- (1+b1*cos(2*pi*(t-phi*period)/period))# seasonality waves
  b <- parms$baseline.txn.rate[t]/ (parms$dur.days1/length.step) # transmission probability per unit time
  beta <-  (b/100)/(sum(yinit.matrix)^(1-q))*c2# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix
  
  beta_a_i <- seasonal.txn * beta/sum(States)
  infectiousN <- I1 + rho1*I2 + rho2*I3 + rho2*I4 + rho2*I_ex # include internal force and external force
  
  lambda <- infectiousN %*% beta_a_i 
  lambda <- as.vector(lambda)
  ##########?????????????????##########################  
  
  dy <- matrix(NA, nrow=N_ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)
  
  period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period #B=annual birth rate
  
  #https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1004591.s016&type=supplementary
  #mu represents aging to the next class
  #um is death rate
  
  Aging.Prop <- c(0,mu[1:(N.ages-1)])
  
  dy[,'M'] <- period.birth.rate*sum(States) - 
    (omega+(mu+um))*M +
    Aging.Prop*c(0,M[1:(N.ages-1)]) 
  
  dy[,'S0'] <- omega*M -
    lambda*S0 -
    (mu + um)*S0 + 
    Aging.Prop*c(0,S0[1:(N.ages-1)]) 
  
  dy[,'I1'] <-   lambda*S0 - 
    (gamma1 + mu + um)*I1 + 
    Aging.Prop*c(0,I1[1:(N.ages-1)]) 
  
  dy[,'S1'] <- gamma1*I1 - 
    sigma1*lambda*S1 - 
    (mu+um)*S1 + 
    Aging.Prop*c(0,S1[1:(N.ages-1)]) 
  
  dy[,'I2'] <- sigma1*lambda*S1 - 
    gamma2*I2-(mu + um)*I2 + 
    Aging.Prop*c(0,I2[1:(N.ages-1)]) 
  
  dy[,'S2'] <- gamma2*I2 - 
    sigma2*lambda*S2 -
    (mu+um)*S2 + 
    Aging.Prop*c(0,S2[1:(N.ages-1)]) 
  
  dy[,'I3'] <- sigma2*lambda*S2 -
    (gamma3 + mu+um)*I3 +  
    Aging.Prop*c(0,I3[1:(N.ages-1)]) 
  
  dy[,'S3'] <- gamma3*I3 +  
    gamma4*I4 -
    sigma3*lambda*S3 -
    (mu + um)*S3 + 
    Aging.Prop*c(0,S3[1:(N.ages-1)]) 
  
  dy[,'I4'] <- sigma3*lambda*S3 - 
    gamma4*I4 - 
    (mu + um)*I4 + 
    Aging.Prop*c(0,I4[1:(N.ages-1)]) 
  
  derivs <- as.vector(dy)
  
  res <- list(derivs)
  
  return(res)
}