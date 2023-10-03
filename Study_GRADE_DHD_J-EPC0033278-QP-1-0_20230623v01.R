###### Script for adding GRADE to data tables #######
# By David Demaree
# ORISE Participant, US EPA, LASMB, Athens, GA.
# 6/23/2023

# Script for adding GRADES for articles to the main spread sheets.

#Removes data from environment
rm(list = ls(all = T))

# Set working directory
setwd("C:/Users/DDEMAREE/OneDrive - Environmental Protection Agency (EPA)/FI to Pathogen RainRWH/R Stuff/Name experiment")

# Datatale to GRADES are on
GRADES <- read.csv("Study_GRADE_DHD_J-EPC0033278-QP-1-0_20230627v02.csv")

# Datasets the GRADES are being added to
# Coliform
RawData <- read.csv("Rawdata/Coliform_STUDYN_DHD_J-EPC0033278-QP-1-0_20230623_v01.csv")

# E. Coli 
RawData <- read.csv("Rawdata/EColi_STUDYN_DHD_J-EPC0033278-QP-1-0_20230623_v01.csv")

# Enterococci
RawData <- read.csv("Rawdata/Enteriococci_STUDYN_DHD_J-EPC0033278-QP-1-0_20230623_v01.csv")

# Salmonella 
RawData <- read.csv("Rawdata/Salmonella_STUDYN_DHD_J-EPC0033278-QP-1-0_20230623_v01.csv")

# Cryptosporidium 
RawData <- read.csv("Rawdata/Cryptosporidium_STUDYN_DHD_J-EPC0033278-QP-1-0_220230623_v01.csv")

# Giardia 
RawData <- read.csv("Rawdata/Giardia_STUDYN_DHD_J-EPC0033278-QP-1-0_20230623_v01.csv")

# Campylobacter 
RawData <- read.csv("Rawdata/Campylobacter_STUDYN_DHD_J-EPC0033278-QP-1-0_20230623_v01.csv")


# Stays the same each Run
b <- nrow(GRADES)

# Rerun each run
a <- nrow(RawData)
# Works now I need it to go through all the names in taxonomy

RawData[,'GRADES'] <- NA
RawData[,'GRADES_Score'] <- NA

for (j in 1:b) {
  for (i in 1:a) {
    if (RawData$Full.Reference[i] == GRADES$Study[j]) {
      RawData$GRADES[i] <- GRADES$Rating.1[j]
      RawData$GRADES_Score[i] <- GRADES$Rating.[j]
    }
  }
}

# This works as a fist pass there are still a bunch of places that need fixing

write.csv(RawData, "Output/Coliform_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv", row.names = FALSE)
write.csv(RawData, "Output/EColi_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv", row.names = FALSE)
write.csv(RawData, "Output/Enteriococci_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv", row.names = FALSE)
write.csv(RawData, "Output/Salmonella_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv", row.names = FALSE)
write.csv(RawData, "Output/Cryptosporidium_STUDYN_DHD_J-EPC0033278-QP-1-0_220230627_v03.csv", row.names = FALSE)
write.csv(RawData, "Output/Giardia_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv", row.names = FALSE)
write.csv(RawData, "Output/Campylobacter_STUDYN_DHD_J-EPC0033278-QP-1-0_20230627_v03.csv", row.names = FALSE)

# After script is finished reading make sure to go through all the data sets because
# there are a couple areas where there were typos in the names

#######################
#Testing

RawData$Full.Reference[1]
[1] "Cox et al. 2005"
GRADES$Study[1]
[1] "Abulreesh et al. 2014"

GRADES$StudyNum <- as.numeric(as.factor(GRADES$Study))

for (i in 1:a) {
  if (RawData$Full.Reference[i] == GRADES$Study[40]) {
    RawData$GRADES[i] <- GRADES$Rating.1[40]
  }
}

# Ran the Script looks like everything works
