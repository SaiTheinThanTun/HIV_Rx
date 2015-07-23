#Cost-utility analysis of treating all HIV patients

library(deSolve)


parameters <- c(
  mui = 1/50, 				#birth ?dP rate of population change?
  muoS = 1/50, 				#death rate for normal people 12.6 per 1000 person years (ref ???trading economics???)
  muoIa = 1/20, 				#death rate for HIV stage I/II
  muoIb = 1/5, 				#death rate for HIV stage III/IV
  time_without_drugs <- 5,
  resistance <- 1 + time_without_drugs,
  trans = 1/10*resistance, 		#transition from Ia to Ib
  R0 <- 5, #Sexual contact (2-5) 	#wikipedia. Treatment could increase by 2- http://cid.oxfordjournals.org/content/44/8/1115.full% , and this does not consider how effect of transmission could be modified by increased longevity (mind you, in so far as transmission rate is increases that MIGHT reduce longevity if life expectancy on therapy is still lower than if not infection)
  #Rx_cov <- .3, 			#ART treatment coverage  change this, it's currently for those who 'shoudl be on treatment', we want all HIV positive.
  Rx_cov <- 6.2/23.5, 			#treatment coverage #avert.org
  
  LE <- 10 				#life_expectancy #average SURVIVAL TIME for HIV positive without treametn is 9-11 years (wikipedia)
  
  
)
time <- 0:50  				#years
costIa = 140 		#Rx cost for HIV stage I/II
costIb = 745.27 

cost_par <- c(Rx_costyr= 140, 	#in USD per person per year
              AIDS_cost =    745.27, 		#cost of treatment of AIDS per patient year in 2015 US$. From page 177 of Confronting AIDS: Public Priorities in a Global Epidemic
              costIa = 140, 		#Rx cost for HIV stage I/II
              costIb = 745.27 			#Rx cost for HIV stage III/IV
)


#Initial values

initP <-1*10^9 # http://populationpyramid.net/sub-saharan-africa/2015/
dP <- 1.6*10^9/50 #rate of pop change over 50 years

initI <- 23.5*10^6 #avert.org
initIa <- initI*.3 #HIV stage I/II
initIb <- initI*.7 #HIV stage III/IV
init_not_S <- 0.15*initP	#From ???http://populationpyramid.net/sub-saharan-africa/2015/??? assuming that those under 10 not infected (not true)

initS <- initP-initI-init_not_S


state <- c(S = initS, Ia=initIa, Ib = initIb, Y = 0, D = 0)
#Cost_Ia = initIa*140, Cost_Ib = initIb*745.27,

#######################################################################

HIV_Rx <- function(t, state, parms) 
{
  with(as.list(c(state, parms)),
       {
         
         P <- S+Ia+Ib
         lam <- R0/50 *((Ia+Ib)/P) #Force of infection
         
         
         #rate of change
         #dCost_Ia <- Ia*costIa
         #dCost_Ib <- Ib*costIb
         dS <- mui*P-muoS*S-lam*S
         dIa <- -muoIa*Ia+lam*S-resistance*Ia
         dIb <- -muoIb*Ib+resistance*Ia
         dY <- 1
         dD <- muoIa*Ia+muoIb*Ib #To calculate Years life lost from infection
         
         # return the rate of change
         list(c(dS, dIa, dIb, dY, dD))
         
         
         
         
       }
  ) 
  
}

############################
#Basecase output
############################
out.base <- ode(y = state, times = time, func = HIV_Rx, parms = parameters)


total_Ia_base <- sum(out.base[,3]) #sum of no. of Ia cases treated each year
total_Ib_base <- sum(out.base[,4]) #sum of no. of Ib cases treated each year

total_ARTcost_base <- total_Ia_base*costIa+total_Ib_base*costIb #need to incoperate current coverage

DW_Ia <-.221
DW_Ib <- .547
DW_ART <- .053

YL_Dis <- total_Ia_base*Rx_cov*DW_ART+total_Ia_base*(1-Rx_cov)*DW_Ia + total_Ib_base*Rx_cov*DW_ART+total_Ib_base*(1-Rx_cov)*DW_Ib
