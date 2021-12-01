I_ex_model <- function(t,y,parms, time.step='month'){
  
  States<-array(y, dim=dim(parms$yinit.matrix))  # MSIS states and the number of age groups
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
  
  omega = 1/(parms$DurationMatImmunityDays/length.step) # rate of waning of maternal immunity
  
  mu= 1/parms$WidthAgeClassMonth
  if(parms$time.step=='week'){
    mu= 1/(WidthAgeClassMonth*4.345) # aging process
  }
  
  gamma1= 1/(parms$dur.days1/length.step)  #converts 1/days to 1/lenth.step # gamma1 stands for rate of recovery from the first infection
  gamma2= 1/(parms$dur.days2/length.step)  # gamma2 stands for rate of recovery from the second infection
  gamma3= 1/(parms$dur.days3/length.step)  # gamma3 stands for rate of recovery from subsequent infection
  gamma4= gamma3  # gamma3 stands for rate of recovery from subsequent infection
  
 #Pull out the states  for the model as vectors
  M <-  States[,'M'] # protected by maternal immunity
  S0 <-  States[,'S0'] # fully suceptible to infections
  I1 <-  States[,'I1'] # first time being infected
  
  S1 <-  States[,'S1'] # partial immunity after first infection
  I2 <-  States[,'I2'] # second time infections
  
  S2 <-  States[,'S2'] # further boosted immunity after second time infections
  I3 <-  States[,'I3'] # third time infections
  
  S3 <-  States[,'S3'] # further boosted immunity after subsequent infections
  I4 <-  States[,'I4'] # any subsequent infections
  
  N.ages <- length(M)# age groups 
  
  ###Check the standardization of beta and overall structure of lambda here
  ##force of transmission##################
  seasonal.txn <- (1+b1*cos(2*pi*(t-phi*period)/period))# seasonality waves
  b <- parms$baseline.txn.rate[t]/ (parms$dur.days1/length.step) # transmission probability per unit time; baseline.txn.rate is a vector stands for the transmissibility at that time point
  # will be lower during COVID mitigation period
  beta <-  (b/100)/(sum(yinit.matrix)^(1-q))*c2# q depends on transmission type (whether depends on population density or not); c2 is the contact matrix
  
  beta_a_i <- seasonal.txn * beta/sum(States)
  infectiousN <- I1 + rho1*I2 + rho2*I3 + rho2*I4 + rho2*I_ex # include internal force and external force # infectious population times their relative infectiousness
  
  lambda <- infectiousN %*% beta_a_i 
  lambda <- as.vector(lambda)
  ##########?????????????????##########################  
  
  dy <- matrix(NA, nrow=N_ages, ncol=ncol(States))
  colnames(dy) <- colnames(States)
  
  period.birth.rate <- log(parms$PerCapitaBirthsYear[t,]+1)/period #B=annual birth rate
  
  #https://journals.plos.org/plospathogens/article/file?id=10.1371/journal.ppat.1004591.s016&type=supplementary
  #mu represents aging to the next class
  #um is death rate
  
  Aging.Prop <- c(0,mu[1:(N.ages-1)]) # aging in
  
  dy[,'M'] <- period.birth.rate*sum(States) - # birth into protected status
    (omega+(mu+um))*M + # waning immunity + aging out + death/immigration/emmigration (um can be either positive or negative based on the the population growth)
    Aging.Prop*c(0,M[1:(N.ages-1)]) # aging in
  
  dy[,'S0'] <- omega*M - # after the maternal immunity wanes, individuals are susceptible to infections
    lambda*S0 - # getting infections
    (mu + um)*S0 + # aging out + death/immigration/emmigration (um can be either positive or negative based on the the population growth)
    Aging.Prop*c(0,S0[1:(N.ages-1)]) # aging in
  
  dy[,'I1'] <-   lambda*S0 - # being infected
    (gamma1 + mu + um)*I1 + # rate of recovery + aging out + death/immigration/emmigration (um can be either positive or negative based on the the population growth)
    Aging.Prop*c(0,I1[1:(N.ages-1)]) # aging in
  
  dy[,'S1'] <- gamma1*I1 - # after recovered from infections, individuals gain partial immunity
    sigma1*lambda*S1 - # which reduced their risk of re-infection; sigma 1 here represnts the reduced risk 
    (mu+um)*S1 + 
    Aging.Prop*c(0,S1[1:(N.ages-1)]) 
  
  dy[,'I2'] <- sigma1*lambda*S1 - 
    gamma2*I2-(mu + um)*I2 + # second time infections will last shorter than the first time infection
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
