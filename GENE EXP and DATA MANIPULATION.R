# Manipulate gene expression data 


#Install all required Libraries- BiocManger, Rstudio-shiny
#Here we have installed all the libraries which required for next step

# script to manipulate gene expression data- (GSE183947)
# setwd("C:/Users/steff/Documents/Data Sciences For Biologist/Data Manipulation-R-prog/Scripts")

#Load Libraries
library(dplyr)
library(tidyverse)
install.packages("BiocManager")                       #if required do this step
library(BiocManager)
BiocManager::install("Biobase")
library(Biobase)

 
# read in the data--
#Had to change the directory because of hex error, later we will change.
setwd("C:/Users/steff/OneDrive/Desktop")             #new current directory

stewd("C:/Users/steff/Documents/R-Projects")
getwd()  #to check where is the r-script directory 

dat <- read.csv(file= "C:/Users/steff/OneDrive/Desktop/Rough R-scripts/GSE183947_fpkm.csv" )
dim(dat)

dat <- read.csv(file = "C:/Users/steff/Documents/Data Sciences For Biologist/Data Munupilation-R-prog/GSE183947_fpkm.csv")
dim(dat)                                             #I got it done, check spellings maybe it caused errors

#get metadata---------
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)
Sys.setenv("VROM_CONNECTION_SIZE" = 131073*1000)
gse


metadata <-pData(phenoData(gse[[1]]))                #command to use phenotype data associated with 1st sample within GEO  
head(metadata)

metadata.subset <- select(metadata, c(1,10, 11, 17))       #select used to make a rw&col of those nos.

metadata.modified <-metadata%>%
  select(1, 10,11, 17) %>%           
  rename(tissue = characteristics_ch1)%>%
  rename(metastasis = characteristics_ch1.1)%>%     # mutate used to remove the string in the column, 
  mutate(tissue = gsub("tissue:", "", tissue))%>%             # mutate can add col or change an existing value in data frame 
  mutate(metastasis = gsub("metastasis:", "", metastasis))
  
   
 head(dat) 


 # reshaping data
 
  dat.long <- dat %>%
   rename(gene = X) %>%                   
   gather(key = 'samples', value = 'FPKM', -gene) 
  
# join dataframes = dat.long + metadata.modified
  dat.long <- dat.long %>%
    left_join(., metadata.modified, by = c("samples" = "description"))
  

  
# explore data
  dat.long %>%
    filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
    group_by(gene, tissue) %>%
    summarise(mean_FPKM = mean(FPKM),
              median_FPKM = median(FPKM)) %>%
    arrange(mean_FPKM)%>%
    head()
    