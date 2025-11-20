MyMatrix <- matrix(c(1:5))
print(MyMatrix)
MySuperMatrix <- matrix(c(1:12), 3, 4, TRUE)
print(MySuperMatrix)
print(MySuperMatrix[3,4] )
print(MySuperMatrix[1,3])

"LeosMatrix" <- matrix(c(1:20), 4, 5, TRUE) 
print("LeosMatrix")

MyMatrix7 <- matrix(c(5:25), 3, 7)
print(MyMatrix7)
print(MyMatrix7[3,5])

####################################################################
# THIS IS CORRECT METHOD TO CODE Matrix IN R-Studio IN 2025

MyMatrixZ1  <- matrix(c(1:8), nrow = 2, ncol = 4, byrow = TRUE)
rownames(MyMatrixZ1) <- c("R1", "R2")
colnames(MyMatrixZ1)  <- c("C1", "C2","C3", "C4")
cat("The 3x2 matrix:\n")     #matrix:\n - is exprsion in arrg in r&C pattern
print(MyMatrixZ1)            #cat func - conc. and shows the values of 1 or 
                             #more R-objects as a single string 




########################################################################

MyArray <- array(c(1:8), c(2,2,2), list(c("ROW1", "ROW2"), c("Col1,", "Col2"), c("Mat1","Mat2")))
print(MyArray)
dim(MyArray) <- c(2,2,2)
print(MyArray)
dim(MyArray) <- c(4, 2)
print(MyArray)
 ####################################################

dim(LeoArray) <- c(2, 2, 2)
print(LeoArray)                 #This did not work because you need to 1st creat a data sets like MyArray


#########################################################

MyFactor <- factor(c("North", "East", "West", "South", "West"))
print(MyFactor)
levels(MyFactor)
nlevels(MyFactor)
 levels(MyFactor) <- c("North", "West", "South", "East")
 print(MyFactor)
 
 NumFactor <- factor(c(1, 2, 3, 4, 1, 2, 1, 3))
 print(NumFactor)
 StrFactor <- factor(NumFactor, labels= c("North", "East", "West", "South"))
 print(StrFactor)
 Appearances["East"]
 print(Appearances)
 
 ###################For creating Basic Data Frame ########################

Name <- c("Jane", "Shane", "Tim")
 Age <- c(42, 20, 54)
Location <- factor(c("North", "West", "East"))
Empdata <- data.frame(Name = Names, Age = Ages, Location = Locations, stringsAsFactors = FALSE)
print(Empdata)

Empdata[1,2]
Empdata$Locations
str(Empdata)

# Summarizing data frame

summary(Empdata)

 # Obtaining Individual Data Frame Data 
Subframe <- data.frame(Names = Empdata$Names, Ages = Empdata$Ages)
print(Subframe)

# Expanding Data Frame - Rows & Columns

 HireDates <- as.Date(c("2012/06/15", "2012/12/12", "2014/07/30"))
 Empdata$HireDate <-  HireDates
print(Empdata)
library(tidyverse)        
NewEmp <- data.frame(Name = "Job", 
                     Age = 34,
                     Location = "South",
                     HireDate = as.Date("2015/02/15"),
                     stringsAsFactors = FALSE)
Empdata <- rbind(Empdata, NewEmp)
print(Empdata)

#HOW To Remove Duplicate ROW From The Data Frame
Empdata_New <- Empdata[!duplicated(Empdata),]
Empdata_New                                    # You can use this way by clicking run instead of typing print