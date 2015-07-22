#Cost-utility analysis of treating all HIV patients
parameters <- c(
  mui = 1/50 , #birth ?dP rate of population change?
  muoS = 1/50 , #death rate for normal people
  muoIa = 1/20 , #death rate for HIV stage I/II
  muoIb = 1/5 , #death rate for HIV stage III/IV
  R0 <- 5, #Sexual contact (2-5) #wikipedia. Treatment could increase by 2- http://cid.oxfordjournals.org/content/44/8/1115.full% , and this does not consider how effect of transmission could be modified by increased longevity (mind you, in so far as transmission rate is increases that MIGHT reduce longevity if life expectancy on therapy is still lower than if not infection)
  #Rx_cov <- .3, #ART treatment coverage  change this, it's currently for those who 'shoudl be on treatment', we want all HIV positive.
  Rx_cov <- 6.2/23.5, #treatment coverage #avert.org
  
  LE		<- 10 #life_expectancy #average SURVIVAL TIME for HIV positive without treametn is 9-11 years (wikipedia)
  
)

time <- 0:50  #years

cost_par <- c(Rx_costyr= 140, #in USD per person per year
              AIDS_cost =    745.27, #cost of treatment of AIDS per patient year in 2015 US$. From page 177 of Confronting AIDS: Public Priorities in a Global Epidemic
              costIa = 745.27*.3, #Rx cost for HIV stage I/II
              costIb = 745.27*.7 #Rx cost for HIV stage III/IV
)


initP <-1*10^9 # http://populationpyramid.net/sub-saharan-africa/2015/
dP <- 1.6*10^9/50 #rate of pop change over 50 years

initI <- 23.5*10^6 #avert.org
initIa <- initI*.3 #HIV stage I/II
initIb <- initI*.7 #HIV stage III/IV
initS <- initP-initI

state <- c(S = initS, Ia=initIa, Ib = initIb, D = 0)



HIV_Rx <- function(t, state, parms) 
{
  with(as.list(c(state, parms)),
       {
         
         P <- S+Ia+Ib
         lam <- R0
         
         #rate of change
         dS <- mui*P-muoS*S-lam*S
         dIa <- -muoIa*Ia+lam*S-trans*Ia
         dIb <- -muoIb*Ib+trans*Ia
         dD <- muoS*S+muoIa*Ia+muoIb*Ib
         
         # return the rate of change
         list(c(dS, dIa, dIb, dY))
         
       }
  ) 
  
}

############################
#Basecase output
############################
out <- ode(y = state, times = time, func = HIV_Rx, parms = parameters)






