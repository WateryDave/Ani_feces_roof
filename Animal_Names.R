# Script for adding names to animals

setwd("C:/Users/DDEMAREE/OneDrive - Environmental Protection Agency (EPA)/FI to Pathogen RainRWH/R Stuff/Name experiment")

Taxonomy <- read.csv("Animal Taxonomy DHD_J-EPC0033278-QP-1-0_20230407v03.csv")
RawData <- read.csv("Salmonella.csv")
RawData <- read.csv("Cryptosporidium.csv")
RawData <- read.csv("Giardia.csv")
RawData <- read.csv("Campylobacter.csv")


a <- nrow(RawData)
b <- nrow(Taxonomy)

RawData$Animal.common.name[1]
[1] "Acadian flycatcher"
Taxonomy$Common.Name.Used[1]
[1] "African fat-tailed gecko"

Taxonomy$AniNum <- as.numeric(as.factor(Taxonomy$Common.Name.Used))


for (i in 1:a) {
  if (RawData$Animal.common.name[i] == Taxonomy$Common.Name.Used[215]) {
  RawData$kingdom[i] <- Taxonomy$kingdom[215]
  RawData$class[i] <- Taxonomy$class[215]
  RawData$order[i] <- Taxonomy$order[215]
  RawData$family[i] <- Taxonomy$family[215]
  RawData$genius[i] <- Taxonomy$genius[215]
}
}

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
