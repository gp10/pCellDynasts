## Create a set of variables by reading the files

#Experimental data
eso00 <- read.table('eso_0day.txt',  sep='\t', header=T)
eso12 <- read.table('eso_12day.txt', sep='\t', header=T)
eso18 <- read.table('eso_18day.txt', sep='\t', header=T)
eso07 <- read.table('eso_7day.txt',  sep='\t', header=T)

scale00 <- read.table('sc_0day.txt',  sep='\t', header=T)
scale12 <- read.table('sc_12day.txt', sep='\t', header=T)
scale18 <- read.table('sc_18day.txt', sep='\t', header=T)
scale07 <- read.table('sc_7day.txt',  sep='\t', header=T)

ear00 <- read.table('ear_0day.txt',  sep='\t', header=T)
ear12 <- read.table('ear_12day.txt', sep='\t', header=T)
ear18 <- read.table('ear_18day.txt', sep='\t', header=T)
ear07 <- read.table('ear_7day.txt',  sep='\t', header=T)

backkm4100 <- read.table('back_km41_0days.txt', sep='\t', header=T)
backkm4105 <- read.table('back_km41_5days.txt', sep='\t', header=T)
backkm4114 <- read.table('back_km41_14days.txt', sep='\t', header=T)

everything <- read.table('all_data.txt', sep='\t', header=T)
fov <- read.table('FoV.txt', sep='\t', header=T)
last <- read.table('last_data.txt', sep='\t', header=T)

#Simulated data (for method validation)
simUnimodal <- read.csv('thousandSimsUnimodal.csv', header=F) #Single population
simUnimodalNoisy <- read.csv('thousandSimsUnimodalWithNoise.csv', header=F) #Single population with noise for initial label taken from eso00
simBimodal <- read.csv('thousandSimsBimodal.csv', header=F) #One 5% subpopulation, 4x slower than the larger
simBimodalNoisy <- read.csv('thousandSimsBimodalWithNoise.csv', header=F) #As above, with noise
simTwoSC <- read.csv('thousandSimsTwoSC.csv', header=F) #One 1/3 subpopulation, 2.6x slower

simFrequent <- read.csv('frequent.csv', header=F) #Demonstration of transient bimodality in division process

#For validating FoV data, we shuffle and shorten original sim data to 200 cells
simUnimodalShort <- head(simUnimodal[sample(nrow(simUnimodal)),],200)
simUnimodalNoisyShort <- head(simUnimodalNoisy[sample(nrow(simUnimodalNoisy)),],200)
simBimodalShort <- head(simBimodal[sample(nrow(simBimodal)),],200)
simBiomodalNoisyShort <- head(simBimodalNoisy[sample(nrow(simBimodalNoisy)),],200)
simTwoSCShort <- head(simTwoSC[sample(nrow(simTwoSC)),],200)

simUnimodalVShort <- head(simUnimodal[sample(nrow(simUnimodal)),],100)
simUnimodalNoisyVShort <- head(simUnimodalNoisy[sample(nrow(simUnimodalNoisy)),],100)
simBimodalVShort <- head(simBimodal[sample(nrow(simBimodal)),],100)
simBiomodalNoisyVShort <- head(simBimodalNoisy[sample(nrow(simBimodalNoisy)),],100)
simTwoSCVShort <- head(simTwoSC[sample(nrow(simTwoSC)),],100)

simUnimodalXShort <- head(simUnimodal[sample(nrow(simUnimodal)),],50)
simUnimodalNoisyXShort <- head(simUnimodalNoisy[sample(nrow(simUnimodalNoisy)),],50)
simBimodalXShort <- head(simBimodal[sample(nrow(simBimodal)),],50)
simBiomodalNoisyXShort <- head(simBimodalNoisy[sample(nrow(simBimodalNoisy)),],50)
simTwoSCXShort <- head(simTwoSC[sample(nrow(simTwoSC)),],50)
