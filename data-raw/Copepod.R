## code to prepare `Copepod` dataset goes here

# Here I process the data off the coast of Japan to make a simplified dataset

MayCop <- read.csv("C:/Users/goods/Desktop/WSBW Synch/Copepod/copepod-2012__biomass-fields/data/copepod-2012__wetmass-m05-qtr.csv")
JunCop <- read.csv("C:/Users/goods/Desktop/WSBW Synch/Copepod/copepod-2012__biomass-fields/data/copepod-2012__wetmass-m06-qtr.csv")
JulCop <- read.csv("C:/Users/goods/Desktop/WSBW Synch/Copepod/copepod-2012__biomass-fields/data/copepod-2012__wetmass-m07-qtr.csv")

LongLatdf1 <- MayCop[, c(1:2)]
LongLatdf2 <- JunCop[, c(1:2)]
LongLatdf3 <- JulCop[, c(1:2)]

LongLatdf <- merge(LongLatdf1, LongLatdf2, by = c("Longitude", "Latitude"))
LongLatdf <- merge(LongLatdf, LongLatdf3, by = c("Longitude", "Latitude"))

# restricting to the region around Japan
LongLatdf <- subset(LongLatdf, Longitude > 120)
LongLatdf <- subset(LongLatdf, Latitude > 20)
LongLatdf <- subset(LongLatdf, Latitude < 60)

MayCop <- merge(MayCop, LongLatdf, by = c("Longitude", "Latitude"), all.y = TRUE)
JunCop <- merge(JunCop, LongLatdf, by = c("Longitude", "Latitude"), all.y = TRUE)
JulCop <- merge(JulCop, LongLatdf, by = c("Longitude", "Latitude"), all.y = TRUE)

colnames(MayCop)
MayCop <- MayCop[, c(1, 2, 4)]
colnames(MayCop) <- c("Longitude", "Latitude", "wtmass5")

colnames(JunCop)
JunCop <- JunCop[, c(1, 2, 4)]
colnames(JunCop) <- c("Longitude", "Latitude", "wtmass6")

colnames(JulCop)
JulCop <- JulCop[, c(1, 2, 4)]
colnames(JulCop) <- c("Longitude", "Latitude", "wtmass7")

Copepod <- merge(MayCop, JunCop, by = c("Longitude", "Latitude"))
Copepod <- merge(Copepod, JulCop, by = c("Longitude", "Latitude"))

usethis::use_data(Copepod, overwrite = TRUE)
