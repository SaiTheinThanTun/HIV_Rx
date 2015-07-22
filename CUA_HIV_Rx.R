#Cost-utility analysis of treating all HIV patients
parameters <- c(
  mui =  #birth ?dP rate of population change?
    muoS =  #death rate for normal people
    muoIa = #death rate for HIV stage I/II
    muoIb = #death rate for HIV stage III/IV
    R0 <- 5, #Sexual contact (2-5) #wikipedia. Treatment could increase by 2- http://cid.oxfordjournals.org/content/44/8/1115.full% , and this does not consider how effect of transmission could be modified by increased longevity (mind you, in so far as transmission rate is increases that MIGHT reduce longevity if life expectancy on therapy is still lower than if not infection)
  #Rx_cov <- .3, #ART treatment coverage  change this, it's currently for those who 'shoudl be on treatment', we want all HIV positive.
  Rx_cov <- 6.2/23.5, #treatment coverage #avert.org
  time <- 50, #years
  <- life_expectancy #average SURVIVAL TIME for HIV positive without treametn is 9-11 years (wikipedia)
  
)

cost_par <- c(Rx_costyr= 140, #in USD per person per year
              AIDS_cost =    745.27, #cost of treatment of AIDS per patient year in 2015 US$. From page 177 of Confronting AIDS: Public Priorities in a Global Epidemic
              
)


initP <-1*10^9 # http://populationpyramid.net/sub-saharan-africa/2015/
dP <- 1.6*10^9/50 #rate of pop change over 50 years

initI <- 23.5*10^6 #avert.org
initIa <- initI*.3 #HIV stage I/II
initIb <- initI*.7 #HIV stage III/IV
initS <- initP-initI






