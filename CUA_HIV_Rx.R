#Cost-utility analysis of treating all HIV patients

#Donald_bump assumes ‘patchy’ pattern (at least temporaly) of drug taking.

library(deSolve)

#Base case

parameters <- c(
  mui = 1/30, 			#http://e360.yale.edu/content/images/0911-haub-population-chart1.html assuming lifetime of 50 years and 50% women and low projection averaged over 50 years (and its a slightly earlier 50 years)
  muoS = 1/100, 				#Assume crude death rate maintained at about 10 (‘trading economics’)
  muoIa = 1/45, 				#death rate for HIV stage I/II
  muoIb = 1/20, 				#death rate for HIV stage III/IV, using data from ‘avert’ from 2011, assuming 70% were in stage III/IV
  #Rx_cov <- .3, 			#ART treatment coverage  change this, it's currently for those who 'should be on treatment', we want all HIV positive.
  #Rx_cov <- 1, 			#treatment coverage #avert.org 6.2/23.5
  LE <- 10 ,				#life_expectancy #average SURVIVAL TIME for HIV positive without treatment is 9-11 years (wikipedia)
  donald_bump <- 53 #has to equal to Y, years of projection
  
)

#Donald bump version

parameters_db <- c(
  mui = 1/30, 			#http://e360.yale.edu/content/images/0911-haub-population-chart1.html assuming lifetime of 50 years and 50% women and average of middle projection.
  muoS = 1/100, 			#Assume crude death rate maintained at about 10 (‘trading economics’)
  muoIa = 1/45, 				#death rate for HIV stage I/II
  muoIb = 1/20, 				#death rate for HIV stage III/IV, using data from ‘avert’ from 2011, assuming 70% were in stage III/IV
  #Rx_cov <- .3, 			#ART treatment coverage  change this, it's currently for those who 'should be on treatment', we want all HIV positive.
  #Rx_cov <- 1, 			#treatment coverage #avert.org 6.2/23.5
  LE <- 10 ,				#life_expectancy #average SURVIVAL TIME for HIV positive without treatment is 9-11 years (wikipedia)
  donald_bump <- 12
  
)

time <- 0:50  				#years

costIa = 140 		#Rx cost for HIV stage I/II
costIb = 745.27 			#Rx cost for HIV stage III/IV

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
init_not_S <- 0.15*initP	#From “http://populationpyramid.net/sub-saharan-africa/2015/” assuming that those under 10 not infected (not true)

initS <- initP-initI-init_not_S


state <- c(S = initS, Ia=initIa, Ib = initIb, Y = 0, D = 0, D0 = 0, IaDis= 0, IbDis=0, IaCost = 0, IbCost = 0)
#Cost_Ia = initIa*140, Cost_Ib = initIb*745.27,

#######################################################################

HIV_Rx <- function(t, state, parms) 
{
  with(as.list(c(state, parms)),
       {
         
         P <- S+Ia+Ib
         lam <-  1.25*((Ia+Ib)/P) ##Infection rate per infected patient per year. Using a value similar to what it was for most of this othe model  http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3512439/
         
         
         
         resistance <- ((Y<donald_bump) + 10*(Y>=donald_bump))
         

    Rx_cov <- ((Y<donald_bump) + (6.2/23.5)*(Y>=donald_bump))
    
    
    
    
    
    #rate of change
    #dCost_Ia <- Ia*costIa
    #dCost_Ib <- Ib*costIb
     
      dS <- mui*P-muoS*S-lam*S
    dIa <- -muoIa*Ia+lam*S-resistance*Ia
    dIb <- (-muoIb*Ib)-(resistance/200*Ib)+(Ia*resistance)
    dY <- 1
    dD <- muoIa*Ia+muoIb*Ib #To calculate Years life lost from infection
    dD0 <- muoIa*Ia*(1-Rx_cov)+muoIb*Ib*(1-Rx_cov) #Deaths without Rx
    dIaDis <- Ia*(1-Rx_cov)*DW_Ia+Ia*Rx_cov*DW_ART
    dIbDis <- Ib*(1-Rx_cov)*DW_Ib+Ib*Rx_cov*DW_ART
    dIaCost <- Ia*Rx_cov*costIa
    dIbCost <- Ib*Rx_cov*costIb
    
    
    # return the rate of change
    list(c(dS, dIa, dIb, dY, dD, dD0, dIaDis, dIbDis, dIaCost, dIbCost))
    
    
    
    
    
    
       }
    ) 

}

############################
#Basecase output
############################
donald_bump <- length(time)+1
out.base <- ode(y = state, times = time, func = HIV_Rx, parms = parameters)
donald_bump <- 12
out.db <- ode(y = state, times= time, func= HIV_Rx, parms= parameters_db)

total_Ia_base <- sum(out.base[,3]) #sum of no. of Ia cases treated each year
total_Ib_base <- sum(out.base[,4]) #sum of no. of Ib cases treated each year

total_Ia_db <- sum(out.db[,3]) #sum of no. of Ia cases treated each year
total_Ib_db <- sum(out.db[,4]) #sum of no. of Ib cases treated each year

total_ARTcost_base <- out.base[length(time),10]+out.base[length(time),11]

total_ARTcost_db <- out.db[length(time),10]+out.db[length(time),11]

#Disability weights, Salomon et al
#HIV: symptomatic, pre AIDS= 0.221
#HIV/AIDS: receiving ART = 0.053
#AIDS: not receiving ART = 0.547

DW_Ia <-.221
DW_Ib <- .547
DW_ART <- .053

YL_Dis_base <- out.base[length(time),8]+out.base[length(time),9]

YL_Dis_db <- out.db[length(time),8]+out.db[length(time),9]

LYlost <- 20 #average life lost without ART Rx
LYlost_base <- LYlost*out.base[length(time),7] #TO CHANGE IF THE TIME HORIZON CHANGES
LYlost_db <- LYlost*out.db[length(time),7]

DALY_base <- LYlost_base+YL_Dis_base
DALY_db <- LYlost_db+YL_Dis_db


Net_cost <- total_ARTcost_db - total_ARTcost_base
Net_DALY <- (DALY_db-DALY_base)

ICER <- Net_cost/Net_DALY


plot(out.base[,10]+out.base[,11], type="l", col="blue")
points(out.db[,10]+out.db[,11], type="l", col="red")

plot(out.base[,8]+out.base[,9], type="l", col="blue")
points(out.db[,8]+out.db[,9], type="l", col="red")

