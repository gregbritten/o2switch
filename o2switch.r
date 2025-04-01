
#calculate pathway-specific ATP yield
#parameters

df <- data.frame(y=c(2,2,2,2,2,2),  #yields
                 N=c(10,10,10,10,10,10),  #electron acceptor concentration
                 c=c(0.1,0.1,0.1,0.1,0.1,0.1),  #cost of enzyme
                 k=c(5,5,5,5,5,5),   #half-saturation constants
                 Vmax=c(1,1,1,1,1,1))   #Vmax

row.names(df) <- c("aerobic LAf","aerobic HAf","narG","nirB","nosZ","nrfA")

PATP <- df$y*(df$Vmax*df$N/(df$k + df$N)) - df$c

#now run a line search

O2switch <- function(par,N,...){
    #evaluate all terms except aerobic terms
    #find largest non-aerobic term X
    #find O2 where X - R(O2) = 0
}

