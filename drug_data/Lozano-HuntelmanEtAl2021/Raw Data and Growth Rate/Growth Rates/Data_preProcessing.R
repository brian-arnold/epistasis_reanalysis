##Confirmed Data organization for 3-4-5 paper... Originally constructed by Nina Singh (12-20-2016). 
##Confirmed and modified for newer versions of R by Natalie Lozano-Huntelman (10-19-2020).

setwd("~/Dropbox/5 Drug Project MSSR/Corrected PreProcessing/110516 Run")


#importing the data
data11 <- read.csv("./GrowthRates/Plate 1.csv", stringsAsFactors = F, header = FALSE)
data12 <- read.csv("./GrowthRates/Plate 2.csv", stringsAsFactors = F, header = FALSE)
data13 <- read.csv("./GrowthRates/Plate 3.csv", stringsAsFactors = F, header = FALSE)
data21 <- read.csv("./GrowthRates/Plate 4.csv", stringsAsFactors = F, header = FALSE)
data22 <- read.csv("./GrowthRates/Plate 5.csv", stringsAsFactors = F, header = FALSE)
data23 <- read.csv("./GrowthRates/Plate 6.csv", stringsAsFactors = F, header = FALSE)
data31 <- read.csv("./GrowthRates/Plate 7.csv", stringsAsFactors = F, header = FALSE)
data32 <- read.csv("./GrowthRates/Plate 8.csv", stringsAsFactors = F, header = FALSE)
data33 <- read.csv("./GrowthRates/Plate 9.csv", stringsAsFactors = F, header = FALSE)


#data14 = as.data.frame(matrix(-111, 16, 22))
#data14 <- read.csv("./GrowthRates/Plate 3a.csv", stringsAsFactors = F, header = FALSE)
#data24= as.data.frame(matrix(-111, 16, 22))
#data24 <- read.csv("./GrowthRates/Plate 4a.csv", stringsAsFactors = F, header = FALSE)
#data34= as.data.frame(matrix(-111, 16, 22))
#data34 <- read.csv("./GrowthRates/Plate 9a.csv", stringsAsFactors = F, header = FALSE)


allData1 <- rbind(data11,data21,data31)
allData2 <- rbind(data12,data22,data32)
allData3 <- rbind(data13,data23,data33)
#allData4 <- rbind(data14,data24,data34)

rowsData <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/rowsToSearch.csv", stringsAsFactors = F, header = TRUE)
colsData <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/colsToSearch.csv", stringsAsFactors = F, header = TRUE)
rowsData <- data.frame(data.matrix(rowsData))
colsData <- data.frame(data.matrix(colsData))
labels <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/labels.csv", stringsAsFactors = F, header = TRUE)

#organizing the data
toView <- NULL
blank <- vector(mode = "character", length = 32)
for(m in 1:243){
  rowNum = m
  elements <- matrix(nrow = 3, ncol = 32);
  for (i in seq(1,31)){
    elements[1,i+1] <- allData1[rowsData[rowNum,i],colsData[rowNum,i]]
    elements[2,i+1] <- allData2[rowsData[rowNum,i],colsData[rowNum,i]]
    elements[3,i+1] <- allData3[rowsData[rowNum,i],colsData[rowNum,i]]
 #   elements[4,i+1] <- allData4[rowsData[rowNum,i],colsData[rowNum,i]]
  }
  labelsForRow <- c(labels[rowNum,31],labels[rowNum,])
  labelsForRow[2:32] <-gsub(pattern = "[123]", replacement = "", labelsForRow[2:32])
  labelsForRow[2:32] <-gsub(pattern = "([ABCDE])([ABCDE])", replacement = "\\1+\\2", labelsForRow[2:32])
  labelsForRow[2:32] <-gsub(pattern = "([ABCDE])([ABCDE])", replacement = "\\1+\\2", labelsForRow[2:32])
  
  #labeling and adding each new set of data to the larger data frame
  toView <- rbind(toView,labelsForRow)
  toView <- rbind(toView,elements)
  toView <- rbind(toView,blank)
}
toView[is.na(toView)] <- ""
#exporting the newly formatted data to a CSV
write.table(toView, file = "./Plates1_9.csv", row.names = FALSE, col.names = FALSE, sep = ",")



#Plates 10-18 2nd Combo of the Run


#importing the data
data11 <- read.csv("./GrowthRates/Plate 10.csv", stringsAsFactors = F, header = FALSE)
data12 <- read.csv("./GrowthRates/Plate 11.csv", stringsAsFactors = F, header = FALSE)
data13 <- read.csv("./GrowthRates/Plate 12.csv", stringsAsFactors = F, header = FALSE)
data21 <- read.csv("./GrowthRates/Plate 13.csv", stringsAsFactors = F, header = FALSE)
data22 <- read.csv("./GrowthRates/Plate 14.csv", stringsAsFactors = F, header = FALSE)
data23 <- read.csv("./GrowthRates/Plate 15.csv", stringsAsFactors = F, header = FALSE)
data31 <- read.csv("./GrowthRates/Plate 16.csv", stringsAsFactors = F, header = FALSE)
data32 <- read.csv("./GrowthRates/Plate 17.csv", stringsAsFactors = F, header = FALSE)
data33 <- read.csv("./GrowthRates/Plate 18.csv", stringsAsFactors = F, header = FALSE)

#data14 = as.data.frame(matrix(-111, 16, 22))
#data14 <- read.csv("./GrowthRates/Plate 11a.csv", stringsAsFactors = F, header = FALSE)
#data24= as.data.frame(matrix(-111, 16, 22))
#data24 <- read.csv("./GrowthRates/Plate 13a.csv", stringsAsFactors = F, header = FALSE)
#data34= as.data.frame(matrix(-111, 16, 22))
#data34 <- read.csv("./GrowthRates/Plate 18a.csv", stringsAsFactors = F, header = FALSE)

allData1 <- rbind(data11,data21,data31)
allData2 <- rbind(data12,data22,data32)
allData3 <- rbind(data13,data23,data33)
#allData4 <- rbind(data14,data24,data34)

rowsData <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/rowsToSearch.csv", stringsAsFactors = F, header = TRUE)
colsData <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/colsToSearch.csv", stringsAsFactors = F, header = TRUE)
rowsData <- data.frame(data.matrix(rowsData))
colsData <- data.frame(data.matrix(colsData))
labels <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/labels.csv", stringsAsFactors = F, header = TRUE)

#organizing the data
toView <- NULL
blank <- vector(mode = "character", length = 32)
for(m in 1:243){
  rowNum = m
  elements <- matrix(nrow = 3, ncol = 32);
  for (i in seq(1,31)){
    elements[1,i+1] <- allData1[rowsData[rowNum,i],colsData[rowNum,i]]
    elements[2,i+1] <- allData2[rowsData[rowNum,i],colsData[rowNum,i]]
    elements[3,i+1] <- allData3[rowsData[rowNum,i],colsData[rowNum,i]]
    #elements[4,i+1] <- allData4[rowsData[rowNum,i],colsData[rowNum,i]]
  }
  labelsForRow <- c(labels[rowNum,31],labels[rowNum,])
  labelsForRow[2:32] <-gsub(pattern = "[123]", replacement = "", labelsForRow[2:32])
  labelsForRow[2:32] <-gsub(pattern = "([ABCDE])([ABCDE])", replacement = "\\1+\\2", labelsForRow[2:32])
  labelsForRow[2:32] <-gsub(pattern = "([ABCDE])([ABCDE])", replacement = "\\1+\\2", labelsForRow[2:32])
  
  #labeling and adding each new set of data to the larger data frame
  toView <- rbind(toView,labelsForRow)
  toView <- rbind(toView,elements)
  toView <- rbind(toView,blank)
}
toView[is.na(toView)] <- ""
#exporting the newly formatted data to a CSV
write.table(toView, file = "./Plates10_18.csv", row.names = FALSE, col.names = FALSE, sep = ",")


#Plates 19-27 3rd Combo of the Run


#importing the data
data11 <- read.csv("./GrowthRates/Plate 19.csv", stringsAsFactors = F, header = FALSE)
data12 <- read.csv("./GrowthRates/Plate 20.csv", stringsAsFactors = F, header = FALSE)
data13 <- read.csv("./GrowthRates/Plate 21.csv", stringsAsFactors = F, header = FALSE)
data21 <- read.csv("./GrowthRates/Plate 22.csv", stringsAsFactors = F, header = FALSE)
data22 <- read.csv("./GrowthRates/Plate 23.csv", stringsAsFactors = F, header = FALSE)
data23 <- read.csv("./GrowthRates/Plate 24.csv", stringsAsFactors = F, header = FALSE)
data31 <- read.csv("./GrowthRates/Plate 25.csv", stringsAsFactors = F, header = FALSE)
data32 <- read.csv("./GrowthRates/Plate 26.csv", stringsAsFactors = F, header = FALSE)
data33 <- read.csv("./GrowthRates/Plate 27.csv", stringsAsFactors = F, header = FALSE)

#data14 = as.data.frame(matrix(-111, 16, 22))
#data14 <- read.csv("./GrowthRates/Plate 20a.csv", stringsAsFactors = F, header = FALSE)
#data24= as.data.frame(matrix(-111, 16, 22))
#data24 <- read.csv("./GrowthRates/Plate 24a.csv", stringsAsFactors = F, header = FALSE)
#data34= as.data.frame(matrix(-111, 16, 22))
#data34 <- read.csv("./GrowthRates/Plate 26a.csv", stringsAsFactors = F, header = FALSE)

allData1 <- rbind(data11,data21,data31)
allData2 <- rbind(data12,data22,data32)
allData3 <- rbind(data13,data23,data33)
#allData4 <- rbind(data14,data24,data34)

rowsData <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/rowsToSearch.csv", stringsAsFactors = F, header = TRUE)
colsData <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/colsToSearch.csv", stringsAsFactors = F, header = TRUE)
rowsData <- data.frame(data.matrix(rowsData))
colsData <- data.frame(data.matrix(colsData))
labels <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/labels.csv", stringsAsFactors = F, header = TRUE)

#organizing the data
toView <- NULL
blank <- vector(mode = "character", length = 32)
for(m in 1:243){
  rowNum = m
  elements <- matrix(nrow = 3, ncol = 32);
  for (i in seq(1,31)){
    elements[1,i+1] <- allData1[rowsData[rowNum,i],colsData[rowNum,i]]
    elements[2,i+1] <- allData2[rowsData[rowNum,i],colsData[rowNum,i]]
    elements[3,i+1] <- allData3[rowsData[rowNum,i],colsData[rowNum,i]]
 #   elements[4,i+1] <- allData4[rowsData[rowNum,i],colsData[rowNum,i]]
  }
  labelsForRow <- c(labels[rowNum,31],labels[rowNum,])
  labelsForRow[2:32] <-gsub(pattern = "[123]", replacement = "", labelsForRow[2:32])
  labelsForRow[2:32] <-gsub(pattern = "([ABCDE])([ABCDE])", replacement = "\\1+\\2", labelsForRow[2:32])
  labelsForRow[2:32] <-gsub(pattern = "([ABCDE])([ABCDE])", replacement = "\\1+\\2", labelsForRow[2:32])
  
  #labeling and adding each new set of data to the larger data frame
  toView <- rbind(toView,labelsForRow)
  toView <- rbind(toView,elements)
  toView <- rbind(toView,blank)
}
toView[is.na(toView)] <- ""
#exporting the newly formatted data to a CSV
write.table(toView, file = "./Plates19_27.csv", row.names = FALSE, col.names = FALSE, sep = ",")


#Plates 28-36 4th Combo of the Run


#importing the data
data11 <- read.csv("./GrowthRates/Plate 28.csv", stringsAsFactors = F, header = FALSE)
data12 <- read.csv("./GrowthRates/Plate 29.csv", stringsAsFactors = F, header = FALSE)
data13 <- read.csv("./GrowthRates/Plate 30.csv", stringsAsFactors = F, header = FALSE)
data21 <- read.csv("./GrowthRates/Plate 31.csv", stringsAsFactors = F, header = FALSE)
data22 <- read.csv("./GrowthRates/Plate 32.csv", stringsAsFactors = F, header = FALSE)
data23 <- read.csv("./GrowthRates/Plate 33.csv", stringsAsFactors = F, header = FALSE)
data31 <- read.csv("./GrowthRates/Plate 34.csv", stringsAsFactors = F, header = FALSE)
data32 <- read.csv("./GrowthRates/Plate 35.csv", stringsAsFactors = F, header = FALSE)
data33 <- read.csv("./GrowthRates/Plate 36.csv", stringsAsFactors = F, header = FALSE)

#data14 = as.data.frame(matrix(-111, 16, 22))
#data14 <- read.csv("./GrowthRates/Plate 34a.csv", stringsAsFactors = F, header = FALSE)
#data24= as.data.frame(matrix(-111, 16, 22))
#data24 <- read.csv("./GrowthRates/Plate 32a.csv", stringsAsFactors = F, header = FALSE)
#data34= as.data.frame(matrix(-111, 16, 22))
#data34 <- read.csv("./GrowthRates/Plate 34a.csv", stringsAsFactors = F, header = FALSE)

allData1 <- rbind(data11,data21,data31)
allData2 <- rbind(data12,data22,data32)
allData3 <- rbind(data13,data23,data33)
#allData4 <- rbind(data14,data24,data34)

rowsData <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/rowsToSearch.csv", stringsAsFactors = F, header = TRUE)
colsData <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/colsToSearch.csv", stringsAsFactors = F, header = TRUE)
rowsData <- data.frame(data.matrix(rowsData))
colsData <- data.frame(data.matrix(colsData))
labels <- read.csv("~/Dropbox/5 Drug Project MSSR/Programs_to generate Elif data/labels.csv", stringsAsFactors = F, header = TRUE)

#organizing the data
toView <- NULL
blank <- vector(mode = "character", length = 32)
for(m in 1:243){
  rowNum = m
  elements <- matrix(nrow = 3, ncol = 32);
  for (i in seq(1,31)){
    elements[1,i+1] <- allData1[rowsData[rowNum,i],colsData[rowNum,i]]
    elements[2,i+1] <- allData2[rowsData[rowNum,i],colsData[rowNum,i]]
    elements[3,i+1] <- allData3[rowsData[rowNum,i],colsData[rowNum,i]]
    #elements[4,i+1] <- allData4[rowsData[rowNum,i],colsData[rowNum,i]]
  }
  labelsForRow <- c(labels[rowNum,31],labels[rowNum,])
  labelsForRow[2:32] <-gsub(pattern = "[123]", replacement = "", labelsForRow[2:32])
  labelsForRow[2:32] <-gsub(pattern = "([ABCDE])([ABCDE])", replacement = "\\1+\\2", labelsForRow[2:32])
  labelsForRow[2:32] <-gsub(pattern = "([ABCDE])([ABCDE])", replacement = "\\1+\\2", labelsForRow[2:32])
  
  #labeling and adding each new set of data to the larger data frame
  toView <- rbind(toView,labelsForRow)
  toView <- rbind(toView,elements)
  toView <- rbind(toView,blank)
}
toView[is.na(toView)] <- ""
#exporting the newly formatted data to a CSV
write.table(toView, file = "./Plates28_36.csv", row.names = FALSE, col.names = FALSE, sep = ",")




