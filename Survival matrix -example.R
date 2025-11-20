# Installing package
install.packages("survival")
# detach the package - survival
detach("package:survival", unload=TRUE)
# Loading package
library(survival)

# Dataset information
?lung
data(lung)

length(lung$time)
length(lung$status)


# Dataset from the survival tool/R- info, description details codes
?lung
head(lung)
class(lung)
dim(lung)
View(lung)

sum(is.na(lung$time))
sum(is.na(lung$status))

# Fitting the survival model
Survival_Function = survfit(Surv(lung$time,
                                 lung$status == 2)~1)
Survival_Function

# Plotting the function
plot(Survival_Function)

