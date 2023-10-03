# Script for adding names to animals

setwd("C:/Users/DDEMAREE/OneDrive - Environmental Protection Agency (EPA)/FI to Pathogen RainRWH/R Stuff/Name experiment")

# File with animal Taxonomy
Taxonomy <- read.csv("Animal Taxonomy DHD_J-EPC0033278-QP-1-0_20230407v03.csv")

Input files
RawData <- read.csv("Salmonella.csv")
RawData <- read.csv("Cryptosporidium.csv")
RawData <- read.csv("Giardia.csv")
RawData <- read.csv("Campylobacter.csv")


a <- nrow(RawData)
b <- nrow(Taxonomy)

Taxonomy$AniNum <- as.numeric(as.factor(Taxonomy$Common.Name.Used))

# Works now I need it to go through all the names in taxonomy

for (j in 1:b) {
  for (i in 1:a) {
    if (RawData$Animal.common.name[i] == Taxonomy$Common.Name.Used[j]) {
      RawData$kingdom[i] <- Taxonomy$kingdom[j]
      RawData$class[i] <- Taxonomy$class[j]
      RawData$order[i] <- Taxonomy$order[j]
      RawData$family[i] <- Taxonomy$family[j]
      RawData$genius[i] <- Taxonomy$genius[j]
    }
  }
}

# This works as a fist pass there are still a bunch of places that need fixing

write.csv(RawData, "AddedAnimalNamesCampy.csv", row.names = FALSE)
write.csv(RawData, "AddedAnimalNamesCryp.csv", row.names = FALSE)
write.csv(RawData, "AddedAnimalNamesGia.csv", row.names = FALSE)
write.csv(RawData, "AddedAnimalNamesSal.csv", row.names = FALSE)
