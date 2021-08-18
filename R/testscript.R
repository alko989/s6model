
#####################################################################################################################################
## Calculate a and test
# library(remotes)
# install_github("alko989/s6model")
# library(s6model)

# Simulate data, characteristics specified in 'parameters'
# d <- simulate((c("a",'Fm'), c(0.22,0.0001), FALSE), 1e4, binsize = 500) # Simulation function - maybe documented? 

d1 <- simulateData3(samplesize = 1e5, params = parameters(c("a",'Fm'), c(0.22,0.0001), FALSE), binsize = 500) 

# Simulates pop with bin size 500, number =e4, parameters allows to set phys mortality (a) + a lot else - False = not log transformed

# Q: What's difference between this 'simulate' and 'simulateData3' function (which is documented)?

# Q: what does this function do? It's always the biomass which changes a little according to the parameters set, but then what?
plot(parameters(c("a",'Fm'), c(0.22,0.0001), FALSE))
plot(parameters(c("a"), c(0.22), FALSE))
plot(parameters(c("Fm"), c(0.0001), FALSE))

s6model:::plot.Parameters # Takes into account all parameters but if not defined it uses the default

#### Notes  ####
# w/Winf = scaling individual weight to asymptotic weight on a 0 to 1 scale
# Wfs = 50% retainment size: depends not only on mesh size but also on net width, skipper behaviour, vessel speed etc. 
# Can't be estimated based on survey data, makes more sense to do this with commercial data.
# if you want to see selectivity: increase Fm in simulation and don't pre-define Fm in estimate_TMB
# Set growth rates in A (growth is hidden in here) based on growth info from literature
# Do we have to supply a and a.sd - sdloga is the actual sd here? sd.a is from makeAssessment

# No Fm predefined - gives exacty the same result as when you do pre-define it as very low.
Fm_0 <- estimate_TMB(d1$df, fnout = NULL, a = 2, sdloga = 0.5, nsample = 1000,map = list()) 
print(Fm_0)
attributes(Fm_0)

plotFit(parameters(c("a",'Fm'), c(0.22,0.0001), FALSE), d1$df)

#  Plot this with the estimated parameters from estimate_TMB
plotFit(attr(Fm_0, "estpars"), d1$df)


# Fm = 0.0001 - a is estimated with predefined search space 
FM_0001 <- estimate_TMB(d1$df, fnout = NULL, a = 2, sdloga = 0.5, Fm = 0.0001,map = list()) 

# nsample = 1000, doesn't do anything here. 


# Q: Why do we pre-define all a, sds and Fm already? For the Robin Hood approach and to set the search space?

str(FM_0001)
print(FM_0001)
attributes(FM_0001)

plotResiduals(FM_0001)

plotFit(object = attr(FM_0001, "estpars"), data = d1$df)

p <- getParams(attr(FM_0001, "estpars"))
data = d1$df

if(add == FALSE) {
  plot(p$w, p$pdfN.approx(p$w), type="l",
       # main="Fitted pdf and histogram of the simulated data",
       xlab="Weight (g)",
       ylab="Probability")
} else {
  lines(p$w, p$pdfN.approx(p$w), col="grey28")
}
points(data$Weight, data$Freq/sum(data$Freq)/diff(c(data$Weight,tail(data$Weight,1))),
       pch=16,cex=1, col="#7E6148B2") #
# lines(density(rep(data$Weight, data$Freq)), col=2, lty=2, lwd=2)
##hist(rep(data$Weight, data$Freq), breaks = 35, add=T, freq=FALSE)
lines(p$w, p$pdfN.approx(p$w), lwd = 2) #col="blue", 
legend("topright", NULL, c("Fitted Probability Density Function (PDF)"),
       lty=1, lwd=2, seg.len=5, bty = "n") #, "Data kernel density" ,"red" col=c("blue"),
invisible(NULL)




plotGrowth(attr(FM_0001, "estpars")) # Plots growth and selectivity
plotGrowth(parameters()) # Uses all the default parameters and plots the model - Not the data
parameters() # Default parameters
#  weight in grams

p <- getParams(attr(FM_0001, "estpars"))
ylim.min <- max(p$g) / 100
ylim.max <- max(p$g)
plot(p$w / p$Winf, p$g, type="n", 
     xlab="",  log="xy", ylab="",
     xlim=c(0.01, 1), ylim=c( ylim.min- ylim.min*0.01, ylim.max + ylim.max*0.6), yaxt="n",xaxt="n")
polygon(c(p$w/p$Winf, 1 ), c(p$psi_m, 0) * ylim.min * 100 + ylim.min, border="lightgrey",col="lightgrey") 

# c(p$psi_m, 0) * ylim.min * 3 + ylim.min

plotMortality(attr(FM_0001, "estpars")) # Plots natural mortality in log scale according to different sizes
# Solid line is natural, spotted is fishing? 
plot(attr(FM_0001, "estpars")) # Plots biomass according to relative weight
