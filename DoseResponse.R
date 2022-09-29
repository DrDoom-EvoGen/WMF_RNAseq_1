
library(drc)
library(dplyr)
library(tidyverse)

## Import and view the datasets
H_DR <- read_csv("Desktop/HAY_DR.csv")
head(H_DR)

C_DR <- read_csv("Desktop/CDA_DR.csv")
head(C_DR)


## Hayden HWM genotype
# Test with raw dry wt first
H.m1<- drm(drywt ~ dose, data = H_DR, fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
#you don't need the 'names = ' argument but it's useful to label the b, c, d, and e parameters until you're familiar with
plot(H.m1, type="all")

# Convert all the responses into a percent of the control response
H_DR <- H_DR %>% 
  mutate(percent_response = drywt/(mean(H_DR$drywt[H_DR$dose==0]))*100)
#check that we now have a response out of 100
head(H_DR)

H_fixed<- drm(percent_response ~ dose, data = H_DR, 
                  fct=LL.4(fixed=c(NA, 0, 100, NA),
                           names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

plot(H_fixed, main="LL.4(fixed=c(NA, 0, 100, NA))")


H_fixed2<- model_fixed<- drm(percent_response ~ dose, data = H_DR, 
                                 fct=LL.3(fixed=c(NA, 100, NA),
                                          names = c("Slope", "Upper Limit", "ED50")))

plot(H_fixed2, main="LL.3(fixed=c(NA, 100, NA))")

# Testing model fits
H_model.LL3<- drm(percent_response ~ dose, data = H_DR, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(H_model.LL3, fctList = list(W1.3(fixed=c(NA, 100, NA)),W1.4(), W2.3(fixed=c(NA, 100, NA)), W2.4(),  LL.4()),linreg=TRUE) 


H_model.W23 <-  drm(percent_response~dose, data=H_DR, fct=W2.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
H_model.W24 <-  drm(percent_response~dose, data=H_DR, fct=W2.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
H_model.LL4 <-  drm(percent_response~dose, data=H_DR, fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
H_model.W14 <-  drm(percent_response~dose, data=H_DR, fct=W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

plot(H_model.LL3, broken = TRUE, xlab="Concentration", ylab="Percent Response", type='all',lty=1, lwd=2)
plot(H_model.W23, add=TRUE,col="orange",lty=1, lwd=2)
plot(H_model.W24, add=TRUE,col="blue",lty=2, lwd=2)
plot(H_model.LL4, add=TRUE,col="forestgreen",lty=2, lwd=2)
plot(H_model.W14, add=TRUE,col="pink",lty=2, lwd=2)



## CDA EWM genotype
# Test with raw dry wt first
C.m1<- drm(drywt ~ dose, data = C_DR, fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
#you don't need the 'names = ' argument but it's useful to label the b, c, d, and e parameters until you're familiar with
plot(C.m1, type="all")

# Convert all the responses into a percent of the control response
C_DR <- C_DR %>% 
  mutate(percent_response = drywt/(mean(C_DR$drywt[C_DR$dose==0]))*100)
#check that we now have a response out of 100
head(C_DR)

C_fixed<- drm(percent_response ~ dose, data = C_DR, 
              fct=LL.4(fixed=c(NA, 0, 100, NA),
                       names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

plot(C_fixed, main="LL.4(fixed=c(NA, 0, 100, NA))")


C_fixed2<- model_fixed<- drm(percent_response ~ dose, data = C_DR, 
                             fct=LL.3(fixed=c(NA, 100, NA),
                                      names = c("Slope", "Upper Limit", "ED50")))

plot(C_fixed2, main="LL.3(fixed=c(NA, 100, NA))")

# Testing model fits
C_model.LL3<- drm(percent_response ~ dose, data = C_DR, fct=LL.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit", "ED50")))
mselect(C_model.LL3, fctList = list(W1.3(fixed=c(NA, 100, NA)),W1.4(), W2.3(fixed=c(NA, 100, NA)), W2.4(),  LL.4()),linreg=TRUE) 


C_model.W23 <-  drm(percent_response~dose, data=C_DR, fct=W2.3(fixed=c(NA, 100, NA), names = c("Slope", "Upper Limit",  "ED50")))
C_model.W24 <-  drm(percent_response~dose, data=C_DR, fct=W2.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
C_model.LL4 <-  drm(percent_response~dose, data=C_DR, fct=LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
C_model.W14 <-  drm(percent_response~dose, data=C_DR, fct=W1.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))

plot(C_model.LL3, broken = TRUE, xlab="Concentration", ylab="Percent Response", type='all',lty=1, lwd=2)
plot(C_model.W23, add=TRUE,col="orange",lty=1, lwd=2)
plot(C_model.W24, add=TRUE,col="blue",lty=2, lwd=2)
plot(C_model.LL4, add=TRUE,col="forestgreen",lty=2, lwd=2)
plot(C_model.W14, add=TRUE,col="pink",lty=2, lwd=2)



## Both genotypes plotted together using their best fit model
bmp(file="/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/DoseResponse.bmp",
    width=6, height=5, units="in", res=600)
plot(C_model.LL3, xlab="mg/L-1 2,4-D", ylab="Percent Control Dry Weight (g)",lty=1, lwd=2, cex = 1.5)
plot(H_model.W23, add=TRUE,col="black",lty=2, lwd=2, pch = 2, cex = 1.5)
legend(0.7, 123, legend=c("EWM", "HWM"),
       col=c("black", "black"), lty=1:2, pch = 1:2, cex=1)
dev.off()

# Calculate ED50 for the 2 genotypes
ED(H_model.W23, 50, interval="delta")
ED(C_model.LL3, 50, interval="delta")

# Plot with confidence intervals and error bars
bmp(file="/Users/Greg/Documents/GitHub/WMF_RNAseq_1/results/DoseResponse_AllPoints.bmp",
    width=6, height=5, units="in", res=600)
plot(C_model.LL3, xlab="mg/L-1 2,4-D", ylab="Percent Control Dry Weight (g)",lty=1, lwd=2, cex = 1.5, type = "all", ylim = c(0,160))
plot(H_model.W23, add=TRUE,col="black",lty=2, lwd=2, pch = 2, cex = 1.5, type = "all")
legend(0.7, 163, legend=c("EWM", "HWM"),
       col=c("black", "black"), lty=1:2, pch = 1:2, cex=1)
dev.off()

plot(C_model.LL3, broken = FALSE, xlab="mg/L-1 2,4-D", ylab="Percent Response",lty=1, lwd=2, type = "confidence", log = "x", confidence.level = 0.95, ylim = c(0,160))
plot(H_model.W23, add=TRUE,col="black",lty=2, lwd=2, pch = 2, type = "confidence", log = "x", confidence.level = 0.95)
legend(0.8, 120, legend=c("EWM", "HWM"),
       col=c("black", "black"), lty=1:2, cex=1)


plot(C_model.LL3, broken = FALSE, xlab="mg/L-1 2,4-D", ylab="Percent Response",lty=1, lwd=2, type = "bars", ylim = c(0,160))
plot(H_model.W23, add=TRUE,col="black",lty=2, lwd=2, pch = 2, type = "bars")
legend(0.8, 120, legend=c("EWM", "HWM"),
       col=c("black", "black"), lty=1:2, cex=1)


