#Previous Growth Rate calculations originally constructed in a separate R script by Nina Singh (11-18-2016)
#Corrected Growth Rate Calculations Natalie Lozano-Huntelman (10-18-2020) OD readings are now log transformed before finding max linear slope
#NLH also made all values <= 0.01 and min ONLY considered positive control slopes > 0.4


setwd("~/Dropbox/5 Drug Project MSSR/Corrected PreProcessing/091916 Run")

library("gdata")
hours <- c(0,4,8,12)
# creates list of all of the files (hour and plate)
expts <- list.files("./MSSR Raw Data/")

# there are 5 hours for each plate
numPlates <- length(expts)/5
#numPlates <- 27

# takes in a plate and subsets the relevant parts, makes them numeric, and normalizes the data
preppedPlate <- function(plate){
  plate <- plate[6:21,1:24]
  plate <- do.call(cbind, lapply(plate, as.numeric))
  avgA <- mean(plate[,1])
  plate <- plate[,2:ncol(plate)]-avgA
  plate[plate<0.01]<- 0.01 #replace all values less than min threshold (0.01) to 0.01
  return(plate)
}

#i=34
for (i in 1:numPlates){
  plateNum = paste("Plate ", i, ".xls", sep = "") # string "Plate i.xls"
  #plateNum = paste("Plate ", i, "a.xls", sep = "") # string "Plate i.xls"
  plateLocations = grep(plateNum, expts) # finds indices containing Plate i in list of plate/hour csvs
  plateNames <- paste("./MSSR Raw Data/",expts[plateLocations], sep = "") # full file names for Plate i for all hours
  
  # finds hour 0 for Plate i and reads in, then preps plate
  hour0 <- read.xls(plateNames[grep("Hour 0",plateNames)], sheet ='Plate', sep = ',', stringsAsFactors = F);
  hour0 <- preppedPlate(hour0)
  
  hour4 <- read.xls(plateNames[grep("Hour 4",plateNames)], sheet ='Plate', sep = ',', stringsAsFactors = F);
  hour4 <- preppedPlate(hour4)
  
  hour8 <- read.xls(plateNames[grep("Hour 8",plateNames)], sheet ='Plate', sep = ',', stringsAsFactors = F);
  hour8 <- preppedPlate(hour8)
  
  hour12 <- read.xls(plateNames[grep("Hour 12",plateNames)], sheet ='Plate', sep = ',', stringsAsFactors = F);
  hour12 <- preppedPlate(hour12)
  
  hour16 <- read.xls(plateNames[grep("Hour 16",plateNames)], sheet ='Plate', sep = ',', stringsAsFactors = F);
  hour16 <- preppedPlate(hour16)
  

  slope40 = (log(hour4)-log(hour0))/4 # calculates slopes by subtracting all data for hour 0 from data from hour 4 then dividing by 4
  #slope40 = matrix(0, 16, 23)
  slope84 = (log(hour8)-log(hour4))/4
  #slope84 = matrix(0, 16, 23)
  slope128 = (log(hour12)-log(hour8))/4
  #slope128 = matrix(0, 16, 23)
  slope1612 = (log(hour16)-log(hour12))/4
  #slope1612 = matrix(0, 16, 23)
  
  slopesPlateI <- c(slope40,slope84,slope128,slope1612) # array of slopes for Plate i
  
  maxSlopePlateI = matrix(nrow = nrow(slope40), ncol = ncol(slope40)) # will store the maxSlope for each well in Plate i
  maxSlopeIndicesPlateI = matrix(nrow = nrow(slope40), ncol = ncol(slope40)) # will store the index of the maxSlope for each well in Plate i
  
  for (j in 1:nrow(slope40)){
    for(k in 1:ncol(slope40)){
      wellSlopes = c(slope40[j,k],slope84[j,k],slope128[j,k],slope1612[j,k]) # finds the slopes for the well in question
      maxSlopePlateI[j,k] = max(wellSlopes) # finds the maximum of the slopes for the well
      maxSlopeIndicesPlateI[j,k] = which.max(wellSlopes) # finds the index of the maximum of the slopes for the well
    }
  }
  
  # normalizes growth rates
  averages <- c(mean(slope40[,ncol(slope40)]),mean(slope84[,ncol(slope84)]),mean(slope128[,ncol(slope128)]),mean(slope1612[,ncol(slope1612)]))
  LastColumn <- maxSlopePlateI[,ncol(maxSlopePlateI)]
  LastColumn <- LastColumn[LastColumn>0.4]
  if (length(LastColumn) == 0){avgLastColumn = 0.6}else{avgLastColumn = mean(LastColumn)}
    
  
  maxSlopePlateI = (maxSlopePlateI[,1:ncol(maxSlopePlateI)-1]/avgLastColumn)*100
  
  # to run through organizing program
  fileI = paste("./GrowthRates/", plateNum, sep = "")
  fileI = paste(substring(fileI, 1, nchar(fileI)-3), "csv", sep = '')
  write.table(maxSlopePlateI, fileI, row.names = FALSE, col.names = FALSE, sep = ',')
  
  # below will be a csv
  fileI = paste("./MaxSlopeIndices/", plateNum, sep = "")
  fileI = paste(substring(fileI, 1, nchar(fileI)-3), "csv", sep = '')
  write.table(maxSlopeIndicesPlateI, fileI, row.names = FALSE, col.names = FALSE, sep = ',')
}
