getwd()
setwd("C:/Users/steff/Documents/R-Projects")
getwd()
######################################################

library("ggplot2")
library("ggpubr")                    # used for publication
library("tidyverse")

theme_set(theme_bw()+theme(legend.position = "top"))  #Setting the theme for the plot

?mtcars()        

head(mtcars)
# 1. For the SCATTER Plot using mtcars dataset

p <- ggplot(mtcars, aes(mpg, wt))+    # Initializes a ggplot using the mtcars dataset
  geom_point()+                       # Layers individual data points to the plot
  geom_smooth(method = lm)+           # Fits and plots a linear regression line
  stat_cor(method= "person", label.x = 20)
print(p)

# Info: the above scatter plot- stat_cor-- is calculates and displays the Person Correlation 
# Coefficient (r) between mpg and wt. Label = 20 is the x- coordinate of the correlation on the plot.
########################################################################################################


# 2. Scatter Plot   (My method scatter plot)
plot(mtcars$wt, mtcars$mpp,
     pch = 19,
     cex =1.5,
     col = "#cc0000",
     main= "MPG as a function of weight of Cars",
     xlab = "Weight(in 1000 pounds)",
     ylab = "MPG")

##########################################################################################################

# 3. Jittered scatter Plot

library(tidyverse)
library(ggpubr)
theme_set(theme_bw()+theme(legend.position= "top"))

# Basic scatter points

r <-ggplot(mpg, aes(cty, hwy))+               #initializing a ggplot thr mtcars dataset
  geom_point(size=0.5)
print(r)

# Jittered points
s <- ggplot(mpg, aes(cty, hwy))+
  geom_jitter(size= 0.5, width = 0.5)
print(s)

head(mtcars)
##################################################################################