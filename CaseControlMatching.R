# Script for balancing a cohort of samples for a two-plate epigenetics arrays
# Matched by different co-variates, factors are weighted


library(tibble)
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(gdata)
library(beepr)

# INPUT FILES
fileIn <- "/Users/timothystone/Desktop/CaseControlMatches/SPITDataDumpMar2020.csv"

# OUTPUT FILES, CHANGE 
fileOut <- paste0("/Users/timothystone/Desktop/CaseControlMatches/CASECONTROL_",sex, normal, ".xls")

# WEIGHTS
ageWeighting        <- 300
cigaretteWeighting  <- 100
drinkWeighting      <-  50
weightWeighting     <- 300
waistWeighting      <- 300
heartburnWeighting  <- 300
PPIWeighting        <- 300
sexWeighting        <- 500
imcWeighting        <- 500
bmiWeighting        <- 300

# SWITCHES
plotSwitch <- FALSE
disableFilters <- 1

# CONSTANTS
ACCEPTABLE_AGE_DIFFERENCE <- 11

# PSEUDO-CONSTANTS, "normals" can be adjusted to include HG or not, or composite NormHV
sexes <- c("Male", "Female")
normals <- c("Normal", "HV", "NDBE", "hgdim")
#normals <- c("NormHV", "NDBE", "hgdim" )

# INITIALISE MAIN VARIABLES
matchedControl <- c()
globalCount <- 1

distanceMatch <- c("NCigsSmoked", "Nbmi", "Nage", "NPPI", "Nhb", "nGenderVec", "nIMvector")


setwd('/Users/timothystone/Desktop/')
IDandGroups <- read.csv('IDandGroups.csv', stringsAsFactors = F) %>% unique() %>% as_tibble()
colnames(IDandGroups) <- c("Subject.number", "FinalDiagnosis")

samples <- read.csv(filein, stringsAsFactors = F) %>% as_tibble()
samples <- samples %>% mutate(age = ifelse(!is.na(Q_P_age) & Q_P_age > REG_ageConsent, Q_P_age, REG_ageConsent))


#Filter a mismatch

samples <- samples %>% filter(Subject.number != "OES-98-0007")


#This is an imputation
samples <- samples %>% mutate(Q_P_weight_kg = replace(Q_P_weight_kg, Subject.number == "OES-11-0034", 66.65))
samples <- samples %>% mutate(Q_P_weight_kg = replace(Q_P_weight_kg, Subject.number == "OES-07-0050", 85.7))
samples <- samples %>% mutate(Q_P_weight_kg = replace(Q_P_weight_kg, Subject.number == "OES-11-0010", 66.65))
samples <- samples %>% mutate(Q_P_height_cm = replace(Q_P_height_cm, Subject.number == "OES-08-0013", 169))



samples <- samples %>% mutate(age = replace(age,Subject.number == "OES-01-1607", 42))
samples <- samples %>% mutate(age = replace(age,Subject.number == "OES-01-1694", 84))
samples <- samples %>% mutate(age = replace(age,Subject.number == "OES-02-0018", 71))


samples <- samples %>% mutate(Q_P_weight_stones = ifelse (Q_P_weight_stones == "10st 6lbs", 10.42, Q_P_weight_stones))

samples <- samples %>% mutate(Q_P_weight_stones = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_weight_stones))
samples <- samples %>% mutate(Q_P_weight_stones = gsub("\\D+", "", Q_P_weight_stones))

samples <- samples %>% mutate(Q_P_weight_stones = as.numeric(Q_P_weight_stones))
samples <- samples %>% mutate(Q_P_weight_stones =  ifelse(Q_P_weight_stones <6 , NA, Q_P_weight_stones))

samples <- samples %>% mutate(Q_P_weight_pounds = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_weight_pounds))
samples <- samples %>% mutate(Q_P_weight_pounds = gsub("\\D+", "", Q_P_weight_pounds))
samples <- samples %>% mutate(Q_P_weight_stones = as.numeric(Q_P_weight_stones))

samples <- samples %>% mutate(Q_P_weight_kg = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_weight_kg))
samples <- samples %>% mutate(Q_P_weight_kg = gsub("\\D+", "", Q_P_weight_kg))
samples <- samples %>% mutate(Q_P_weight_kg = as.numeric(Q_P_weight_kg))



samples <- samples %>% mutate(Q_P_height_feet = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_height_feet))
samples <- samples %>% mutate(Q_P_height_feet = gsub("\\D+", "", Q_P_height_feet))
samples <- samples %>% mutate(Q_P_height_feet = as.numeric(Q_P_height_feet))
samples <- samples %>% mutate(Q_P_height_inches = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_height_inches))
samples <- samples %>% mutate(Q_P_height_inches = gsub("\\D+", "", Q_P_height_inches))
samples <- samples %>% mutate(Q_P_height_inches = as.numeric(Q_P_height_inches))

samples <- samples %>% mutate(Q_P_height_cm = gsub("(\\d+)(\\.\\d+)", "\\1", Q_P_height_cm))
samples <- samples %>% mutate(Q_P_height_cm = as.numeric(Q_P_height_cm))
samples <- samples %>% mutate(Q_P_height_cm = ifelse(Q_P_height_cm < 100, NA, Q_P_height_cm))


samples <- samples %>% mutate(height = ifelse(grepl("\\d", Q_P_height_feet) & !grepl("\\d", Q_P_height_cm),
                                              round(as.numeric(Q_P_height_feet) * 30.48 + (as.numeric(Q_P_height_inches) * 2.54),2), NA))
samples <- samples %>% mutate(height = ifelse(grepl("\\d", height), height, Q_P_height_cm))



samples <- samples %>% mutate(weight = ifelse(grepl("\\d", Q_P_weight_stones) & !grepl("\\d", Q_P_weight_kg),
                                              round(as.numeric(Q_P_weight_stones) * 6.35029 + (as.numeric(Q_P_weight_pounds) * 0.453592),2), ""))
samples <- samples %>% mutate(weight = ifelse(grepl("\\d", Q_P_weight_kg), Q_P_weight_kg, weight))
samples <- samples %>% mutate(weight = ifelse(weight == "" | weight == 0, NA, weight))
thisweight <- samples %>% select(weight) %>% unlist %>% as.numeric()
samples <- samples %>% mutate(weight = thisweight)
samples <- samples %>% mutate(weight = ifelse(weight > 1000, NA, weight))
samples <- samples %>% mutate(weight = as.numeric(weight))


samples <- samples %>% filter(!is.na(weight))
samples <- samples %>% mutate(bmi = ifelse(grepl("\\d", height) & grepl("\\d", weight), weight/(((height/100))^2), NA))


#waist <- read.csv('SPITWaistSizes.csv', stringsAsFactors = F) %>% as_tibble()
#waist <- mutate(waist, PINF_waist_inches = ifelse(PINF_waist_cm < 50,  PINF_waist_cm, PINF_waist_inches))

#waist <- mutate(waist, final = ifelse(!is.na(PINF_hip_inches), PINF_hip_inches * 2.54 ,PINF_hip_cm))

samples <- left_join(IDandGroups, samples)

ageBoundary <- c(17, 25, 35, 45, 59, 999)
ageRange <- c(8, 10, 10, 15)
averageUnits <- c(5, 15, 25, 35, 45)
alcoholNames <- colnames(samples[,grep("Q_ALC", colnames(samples))])[4:8]
alcoholGreps <- c("teetotal", "light drinker", "moderate drinker", "31", "very heavy drinker")
samples <- samples %>% mutate(CumulAlc = 0)
ages <- select(samples, age) %>% unlist() 
ageBins <- .bincode (ages, ageBoundary)
unitsDrunk <- samples %>% select(CumulAlc) %>% unlist()

for (i in 1:5) {
  drinkYears <- rep(0, length(unitsDrunk))
  drinkYears[which(ageBins == i)] <- ages[which(ageBins == i)] - ageBoundary[i]
  drinkYears[which(ageBins > i)] <- ageRange[i] 
  alcoholSurvey <- select(samples, !!alcoholNames[i]) %>% unlist %>% as.character
  
  for (j in 1:5) {
    whichUnits <- grep(alcoholGreps[j], alcoholSurvey)
    unitsDrunk[whichUnits] <- unitsDrunk[whichUnits] + averageUnits[j] * drinkYears[whichUnits]
  }
}

samples <- mutate(samples, UnitsDrunk = unitsDrunk)

# Smoking create a measure of cumulative cigarettes smoked
ageBoundary <- c(15, 25, 35, 45, 59, 999)
ageRange <- c(10, 10, 10, 15)
averageCigs <- c(5, 15, 30, 50)
smokeNames <- colnames(samples[,grep("Q_SMK", colnames(samples))])[5:9]
smokeGreps <- c("10", "20", "21", "More than")
samples <- samples %>% mutate(CumulSmoke = 0)
ages <- select(samples, age) %>% unlist() 
ageBins <- .bincode (ages, ageBoundary)
cigsSmoked <- samples %>% select(CumulSmoke) %>% unlist()

for (i in 1:5) {
  smokeYears <- rep(0, length(cigsSmoked))
  smokeYears[which(ageBins == i)] <- ages[which(ageBins == i)] - ageBoundary[i]
  smokeYears[which(ageBins > i)] <- ageRange[i] 
  smokeSurvey <- select(samples, !!smokeNames[i]) %>% unlist %>% as.character
  
  for (j in 1:4) {
    whichSmoke <- grep(smokeGreps[j], smokeSurvey)
    cigsSmoked[whichSmoke] <- cigsSmoked[whichSmoke] + averageCigs[j] * smokeYears[whichSmoke]
  }
}

samples <- mutate(samples, CigsSmoked = cigsSmoked)

samples <- mutate(samples, PPIfreq = as.character(Q_HB_ppi_medication_frequency))
samples <- mutate(samples, PPIfreq = gsub("(.+) $", "\\1", PPIfreq, perl=T))
samples <- mutate(samples, PPIfreq = case_when(
  PPIfreq == "Daily" ~ 365,
  PPIfreq == "Few times a month" ~ 24,
  PPIfreq == "Few times a week" ~ 150,
  PPIfreq == "Few times a year" ~ 6,
  PPIfreq == "Never" ~ 0,
  PPIfreq == "" ~ 0))

samples <- mutate(samples, PPIbeg = as.character(Q_HB_ppi_medication_begin))
samples <- mutate(samples, PPIbeg = case_when(
  PPIbeg == "1 to 5 years" ~ 2.5,
  PPIbeg == "10 to 20 years" ~ 15,
  PPIbeg == "5 to 10 years" ~ 7.5,
  PPIbeg == "6 months to 1 years" ~ 0.75,
  PPIbeg == "Less than 6 months" ~ 0.25,
  PPIbeg == "More than 20 years" ~ 30,
  PPIbeg == "" ~ 0))
samples <- mutate(samples, PPIbeg = ifelse(is.na(PPIbeg), 0, PPIbeg))

# Estimate that if someone previous took PPIs they did it for a fifth of their adult life
samples <- mutate(samples, adultLife = age - 18)
samples <- mutate(samples, adultDoseLife = adultLife / 5)
samples <- mutate(samples, PPIprev = as.character(Q_HB_previous_ppi_medication_frequency))


samples <- mutate(samples, PPIprev = case_when(
  PPIprev == "Daily" ~ 365,
  PPIprev == "Few times a month" ~ 24,
  PPIprev == "Few times a week" ~ 150,
  PPIprev == "Few times a year" ~ 6,
  PPIprev == "Never" ~ 0,
  PPIprev == "" ~ 0))
samples <- mutate(samples, ifelse(is.na(PPIprev), 0, PPIprev))
samples <- mutate(samples, PPIprev = ifelse(PPIfreq == 0, PPIprev, 0))
samples <- mutate(samples, PPI = (PPIfreq * PPIbeg) + (PPIprev * adultDoseLife))


samples <- mutate(samples, hb = as.character(Q_HB_heartburn_previous_frequency))
samples <- mutate(samples, hb = case_when(
  hb == "Daily" ~ 365,
  hb == "Few times a month" ~ 24,
  hb == "Few times a week" ~ 150,
  hb == "Few times a year" ~ 6,
  hb == "Never" ~ 0,
  hb == "" ~ 0))


# Count the site

nrow(samples)
ncol(samples)

# Site information
select(samples, Site) %>% table()

# Check Screening numbers are all unique
select(samples, Screening.number) %>% distinct() %>% count()

# Check Subject numbrers
samples %>% select(Subject.number) %>% duplicated() %>% sum

diagnosisRef <- grep("PAT_highest_grade_diagnosis_ever$", colnames(samples), perl=T)
diagnosis <- colnames(samples)[diagnosisRef]
diagnosisNo <- diagnosisRef

groups <- samples[,diagnosisRef] %>% unlist()

gvector <- rep(NA, length(groups))
glevels <- samples[,diagnosisRef] %>% unlist %>% levels
case <- c(glevels[3], glevels[4], glevels[5], glevels[6])
control <- c(glevels[7], glevels[8])
unknown <- c(glevels[1], glevels[2])
other <- c(glevels[9])

gvector[groups %in% case] <- "Case"
gvector[groups %in% control] <- "Control"
gvector[groups %in% unknown] <- "DontKnow"
gvector[groups %in% other] <- "Other"

samples <- mutate(samples, Group = gvector)
everyone <- filter(samples, !is.na(FinalDiagnosis))

cases <- filter(everyone, Group == "Case")
controls <- filter(everyone, Group == "Controls")

everyone <- mutate(everyone, diagnosis = as.character(PAT_highest_grade_diagnosis_ever))
everyone <- mutate(everyone, diagnosis = as.character(diagnosis))



everyone <- mutate(everyone, diagnosis = 
                     case_when(
                       diagnosis == "High Grade Dysplastic Barrett's oesophagus" ~ "High Dys",
                       diagnosis == "Low Grade Dysplastic Barrett's oesophagus" ~ "Low Dys",
                       diagnosis == "Intramucosal Oesophageal Adenocarcinoma" ~ "Intramuc. OAC",
                       diagnosis == "Invasive Oesophageal Adenocarcinoma" ~ "Invasive OAC",
                       diagnosis == "Non dysplastic Barrett's oesophagus" ~ "NDBE",
                       grepl("now", diagnosis) ~ "Don't know",
                       diagnosis== "Normal (no Barrett's)" ~ "Normal",
                       grepl("Other ", diagnosis) ~ "Other")) 


everyone <- mutate(everyone, diagnosis = factor(diagnosis, levels = c("Normal", "NDBE", "Low Dys", "High Dys",
                                                                      "Intramuc. OAC", "Invasive OAC", "DontKnow", "Other")))


#everyone <- filter(everyone, !is.na(diagnosis))

if (!disableFilters) {
  
  #Race Filter
  everyone <- filter(everyone, grepl("White", PINF_ethnicity) | grepl("White", Q_P_ethnicity) )
  
  
  #Low-grade Dysplasia Filter
  #everyone <- filter(everyone, diagnosis != "Low Dys")
  
  
  # Only people with saliva samples
  everyone <- filter(everyone, grepl("aliva", QS_samples_collected))
  
  x <- everyone %>% filter(diagnosis=="Normal") %>% select(age) %>% unlist()
  y <- everyone %>% filter(diagnosis %in% c("Intramuc. OAC", "Invasive OAC")) %>% select(age) %>% unlist()
  
}

everyone <- everyone %>% mutate(Gender = as.character(PINF_gender))

if (plotSwitch == T) {
  jpeg("AllGroupsWEIGHT.jpeg")
  ggplot(data = everyone, aes(x=diagnosis, y = CigsSmoked, col=diagnosis)) +  geom_violin()
  dev.off()
  
  
  jpeg("AllGroupsAGE.jpeg")
  ggplot(data = everyone, aes(x=diagnosis, y = UnitsDrunk, col=diagnosis)) +  geom_boxplot(notch = F)
  dev.off()
}

agevector <- c()

for (i in 20:100) {
  male <- everyone %>% filter(Gender == "Male" & age <= i) %>% nrow() %>% unlist()
  female <- everyone %>% filter(Gender == "Female" & age <= i) %>% nrow() %>% unlist()
  
  agevector <- rbind(agevector, c(i,male, female))
}

colnames(agevector) <- c("Age", "Male", "Female")
agevector <- as.data.frame(agevector)
agevector <- gather(agevector, "Sex","No.of.Younger.People",  Male:Female) 

if (plotSwitch == T) {
  jpeg("AgePlot.jpeg")
  ggplot(data = agevector, aes(x=Age, y=No.of.Younger.People, col=Sex)) + geom_line(size=2)
  dev.off()
  
  jpeg("AllGroupsVIOLIN.jpeg")
  ggplot(data = everyone, aes(x=diagnosis, y = weight, col=diagnosis)) + geom_violin()
  dev.off()
}


cancer <- everyone %>% filter(FinalDiagnosis == "Cancer") %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, age, weight, bmi,
                                                                     CigsSmoked, PPI, hb, QS_samples_collected)
norms <- everyone %>% filter(FinalDiagnosis == "Normal") %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, age, weight, bmi,
                                                                    CigsSmoked, PPI, hb,QS_samples_collected)

ndbe <- everyone %>% filter(FinalDiagnosis == "NDBE") %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, age, weight, bmi,
                                                                 CigsSmoked, PPI, hb, QS_samples_collected)

hgd <- everyone %>% filter(FinalDiagnosis == "HG") %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, age, weight, bmi,
                                                              CigsSmoked, PPI, hb, QS_samples_collected)

hv <- everyone %>% filter(FinalDiagnosis == "HV") %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, age, weight, bmi,
                                                             CigsSmoked, PPI, hb, QS_samples_collected)

im <- everyone %>% filter(FinalDiagnosis == "IM") %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, age, weight, bmi,
                                                             CigsSmoked, PPI, hb, QS_samples_collected)


allControls <- rbind(norms, ndbe, hv, hgd, im)
BigGroup <- rbind(cancer, norms, ndbe)

normalise <- function(x, y=100) { x * y / max(x, na.rm=T) }


cancer <- cancer %>% filter(!is.na(Gender) & (!(is.na(age)))) 
norms <- norms %>% filter(!is.na(Gender) & (!(is.na(age)))) 
ndbe <- ndbe %>% filter(!is.na(Gender) & (!(is.na(age)))) 
hgd <- hgd %>% filter(!is.na(Gender) & (!(is.na(age)))) 
hv <- hv %>% filter(!is.na(Gender) & (!(is.na(age)))) 
im <- im %>% filter(!is.na(Gender) & (!(is.na(age)))) 

hgdim <- rbind(hgd, im)
normHV <- rbind(norms, hv)

everyMatch <- c()

for (bigloop in 1:2) {
  
  if (bigloop == 2) {
    normals <- c("Normal", "NDBE")
  }
  
  
  # loops for each sex / normal
  for (normal in normals) {
    
    if (bigloop == 1) {
      
      if (normal == "Normal")
      {
        all <- rbind(cancer, norms)
      } 
      
      if (normal == "NDBE")
      {
        all <- rbind(cancer, ndbe)
      } 
      if (normal == "hgdim")
      {
        all <- rbind(cancer, hgdim)
      } 
      
      if (normal == "NormHV")
      {
        all <- rbind(cancer, normHV)
      }
      
      if (normal == "HV")
      {
        all <- rbind(cancer, hv)
      } 
    }
    
    if (bigloop == 2) {
      if (normal == "Normal")
      {
        all <- rbind(hgd, norms)
      } 
      
      if (normal == "NDBE")
      {
        all <- rbind(hgd, ndbe)
      } 
      
      
    }
    
    all <- all %>% mutate(Nweight = normalise (weight, weightWeighting))
    all <- all %>% mutate(NCigsSmoked = normalise (CigsSmoked, cigaretteWeighting))
    all <- all %>% mutate(Nbmi = normalise (bmi, bmiWeighting))
    
    #all <- all %>% mutate(NUnitsDrunk = normalise (UnitsDrunk, drinkWeighting))
    all <- all %>% mutate(Nage = normalise (age, ageWeighting))
    all <- all %>% mutate(NPPI = normalise(PPI, PPIWeighting))
    all <- all %>% mutate(Nhb = normalise(hb, heartburnWeighting))
    all <- mutate(all, nGenderVec = normalise(as.numeric(factor(Gender))-1, sexWeighting))
    
    all <- mutate(all, nIMvector = ifelse(FinalDiagnosis == "IM" | FinalDiagnosis == "Cancer", imcWeighting, 0))
    
    
    
    all <- mutate(all, clusterinfo = paste(Subject.number,round(age,0), paste(round(weight,0), "kg", sep=""),
                                           CigsSmoked,FinalDiagnosis,sep="_"))
    all <- mutate(all, weight = round(weight, 1))
    all <- mutate(all, bmi = round(bmi, 1))
    
    
    
    
    #plot(hclust(dist(all[,10:13])), labels= all[,"clusterinfo"] %>% unlist(), col.lab=c("red") )
    
    if (bigloop == 1) {
      cancerName <- "Cancer"
    } else {
      cancerName <- "HG"
    }
    
    nAllCase <- filter(all, FinalDiagnosis == cancerName) %>% nrow %>% unlist()
    nAllControl <- filter(all, FinalDiagnosis != cancerName) %>% nrow %>% unlist()
    
    
    
    allCopy <- all
    
    for (i in 1:nAllCase) {
      
      
      if (normal == "hgdim") {
        ncontrol <- filter(allCopy, grepl("HG|IM", FinalDiagnosis)) %>% nrow %>% unlist()
      }
      
      if (normal == "NormHV") {
        ncontrol <- filter(allCopy, grepl("Normal|HV", FinalDiagnosis)) %>% nrow %>% unlist()
      }
      
      
      if (normal != "hgdim" & normal != "NormHV")
      { 
        ncontrol <- filter(allCopy, grepl(normal, FinalDiagnosis)) %>% nrow %>% unlist()
      }
      ncase <- filter(allCopy, grepl(cancerName, FinalDiagnosis)) %>% nrow %>% unlist()
      
      cat("ncontrol is ", ncontrol,", and ncase is ", ncase, "\n")
      
      if (ncase == 0) { stop() }
      if (ncase >= 1 & ncontrol >= 1) {
        
        
        
        matchStatus <- paste(normal,": MATCHED", sep="")
        distancePairs <- dist(allCopy[, distanceMatch]) %>% as.matrix()
        distancePairsMax <- max(distancePairs)[1][[1]]
        
        # Stop the subjects matching themselves, make the central diagonal huge
        distancePairs[row(distancePairs) == col(distancePairs)]  <- distancePairsMax * 1000
        
        # Stop any Cancer-cancer or Control-control matching, make them huge
        distancePairs[col(distancePairs) <= ncase] <- distancePairsMax * 1000
        distancePairs[row(distancePairs) > ncase] <- distancePairsMax * 1000
        
        # the magic number now is the shortest distance
        magic <- which(distancePairs == min(distancePairs), arr.ind = T)[1,]
        
        archiveMagic <- sort(distancePairs)
        archiveMagic <- archiveMagic[archiveMagic < distancePairsMax * 1000]
        
        closeCases <- lapply(archiveMagic, function(x) { which(distancePairs == x, arr.ind=T)}) %>%
          lapply(., function(x) { x[1,2]}) %>% unlist()
        
        thisCases <- allCopy[closeCases, ] %>% select(Subject.number)
        everyMatch <- rbind(everyMatch, cbind(thisCases, archiveMagic))
        
        magicmatch <- rbind(allCopy[magic[1],1:11], allCopy[magic[2],1:11])
        magicmatch <- cbind(magicmatch)
        
        magicCase <- magicmatch[1,"Subject.number"] %>% unlist() %>% as.character()
        
        
        matchStatus <- normal
        
        
        if (abs(magicmatch[1,"age"] - magicmatch[2,"age"]) > ACCEPTABLE_AGE_DIFFERENCE) {
          matchStatus <- paste(matchStatus, ": Bad Age MATCH", sep="")
        }
        
        if (magicmatch[1,"Gender"] != magicmatch[2,"Gender"]) {
          matchStatus <- paste(matchStatus, ": Gender Mismatch", sep="")
        }
        
        if (globalCount == 1) {
          #write.table(magicmatch, fileout, sep="\t", row.names=F)
          
          matchedControl <- cbind(GlobalMatch = globalCount, MatchStatus =  matchStatus, magicmatch)
          globalCount <- globalCount + 1
          allCopy <- allCopy[-c(magic[1], magic[2]),]
          next()
        }
        
        if(nrow(filter(matchedControl, Subject.number == magicCase)) == 0) {
          
          matchAdd <-  cbind(GlobalMatch = globalCount, MatchStatus =  matchStatus, magicmatch)
          matchedControl <- rbind(matchedControl, matchAdd)
          globalCount <- globalCount + 1
          
        } else {
          thismatch <- filter(matchedControl, Subject.number == magicCase) %>% select(GlobalMatch) %>% unlist() %>% as.integer()
          matchAdd = cbind(GlobalMatch = thismatch, MatchStatus = matchStatus, magicmatch[2,])
          matchedControl <- rbind(matchedControl, matchAdd)
        }
        
        cat("globalCount is ", globalCount, " and I'm removing rows", unlist(magic[1]), " and ", unlist(magic[2]), "\n")
        allCopy <- allCopy[-c(magic[1], magic[2]),]
      }
      
      #if (ncontrol == 0 & ncase != 0) {
      #  
      #  matchStatus <- paste("NO", normal, "MATCH")
      #  magicCase <- allCopy[1,"Subject.number"] 
      #  if(nrow(matchedControl %>% select(Subject.number) %>% intersect(magicCase)) == 0) {
      #    matchAdd <-  cbind(GlobalMatch = globalCount, MatchStatus =  matchStatus, allCopy[1,1:11])
      #    matchAdd <-  matchAdd %>% add_row( GlobalMatch = globalCount, MatchStatus = matchStatus, 
      #                                       Gender = sex, diagnosis = normal, Group = "Control")
      #    matchedControl <- rbind(matchedControl, matchAdd)
      #    globalCount <- globalCount + 1
      #  } else {
      #    
      #    thisGlobal <- matchedControl %>% filter(Subject.number == unlist(magicCase)) %>% 
      #      select(GlobalMatch) %>% as.integer
      #    
      #    matchedControl <- matchedControl %>% add_row(GlobalMatch = thisGlobal, MatchStatus = matchStatus,
      #                                                 Gender = sex, diagnosis = normal, Group = "Control")
      #  }
      #  allCopy <- allCopy[-1,]
      #  
      #  
      #}
      
    }
    
    #colnames(matched) <- colnames(all)[1]
    #plot(hclust(dist(all[,10:13])), labels= all[,"clusterinfo"] %>% unlist(), col.lab=c("red") )
    #dend <- as.dendrogram(hclust(dist(all[,4:7])), labels= all[,"clusterinfo"] %>% unlist())
  }
  
  #matchedControl <- matchedControl %>% mutate(MatchStatus = 
  #                                              ifelse (Group == "Case", FinalDiagnosis, as.character(MatchStatus)))
  
  #matchedControl <- arrange(matchedControl, GlobalMatch, Group) 
  matchedControl <- as.tibble(matchedControl)
  
  matchedControl <- mutate(matchedControl, MatchStatus = as.character(MatchStatus))
  
  #matchedControl <- arrange(matchedControl, GlobalMatch, Group)
  matchedControl <- arrange(matchedControl, GlobalMatch) 
  matchedControl <- mutate(matchedControl, MatchStatus = ifelse(MatchStatus == "Normal" & FinalDiagnosis == cancerName, FinalDiagnosis, MatchStatus))
  
  
  missing <- matchedControl %>% group_by(GlobalMatch) %>% count %>% filter(n !=5) %>% select(GlobalMatch)
  
  if (nrow(missing) != 0) {
    missingNDBE <- left_join(missing, matchedControl) %>% filter(diagnosis == "Normal") %>% select(GlobalMatch)
    missingNormal <- left_join(missing, matchedControl) %>% filter(diagnosis == "NDBE") %>% select(GlobalMatch)
    
    matchedControl <- add_row(matchedControl, GlobalMatch = unlist(missingNDBE), Group = "Control", MatchStatus = "NO NDBE MATCH", diagnosis = "NDBE")
    matchedControl <- add_row(matchedControl, GlobalMatch = unlist(missingNormal), Group = "Control", MatchStatus = "NO NORMAL MATCH", diagnosis = "Normal")
  }
  
  
  used4matching <- matchedControl %>% filter(!is.na(Subject.number)) %>% select(Subject.number) %>% unlist()
  
  everyMatch <- everyMatch[order(everyMatch[,2]),]
  
  if (bigloop == 1) {
    runnersUp <- everyMatch[,1][!everyMatch[,1] %in% used4matching] %>% unique %>% as_tibble()
    colnames(runnersUp) <- "Subject.number"
    
    runnersUp <- left_join(runnersUp, allControls)
    runnersUp <- add_column(runnersUp, GlobalMatch = -999, MatchStatus = "RunnerUp")
    
    hgRunnersUp <- runnersUp %>% filter(FinalDiagnosis == "HG") %>% select(Subject.number)
    normsRunnersUp <- runnersUp %>% filter(FinalDiagnosis == "Normal") %>% select(Subject.number)
    ndbeRunnersUp <- runnersUp %>% filter(FinalDiagnosis == "NDBE") %>% select(Subject.number)
    runnersRef <- rbind(hgRunnersUp, normsRunnersUp, ndbeRunnersUp)
    
    
    norms <- everyone %>% filter(FinalDiagnosis == "Normal") %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, age, weight, bmi,
                                                                        CigsSmoked, PPI, hb, QS_samples_collected)
    norms <- left_join(normsRunnersUp, norms)
    
    ndbe <- everyone %>% filter(FinalDiagnosis == "NDBE") %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, age, weight, bmi,
                                                                     CigsSmoked, PPI, hb, QS_samples_collected)
    ndbe <- left_join(ndbeRunnersUp, ndbe)
    
    hgd <- everyone %>% filter(FinalDiagnosis == "HG") %>% select(Subject.number, FinalDiagnosis, Group, diagnosis, Gender, age, weight, bmi,
                                                                  CigsSmoked, PPI, hb, QS_samples_collected)
    hgd <- left_join(hgRunnersUp, hgd)
  }
  
}


#finalList <- rbind(matchedControl2, runnersUp2)
matchedControl <- mutate(matchedControl, MatchStatus = ifelse(grepl("Normal", MatchStatus) & FinalDiagnosis == "HG", 
                                                              FinalDiagnosis, MatchStatus))



matchedControl <- mutate(matchedControl, MatchStatus = ifelse(MatchStatus == "hgdim", FinalDiagnosis, MatchStatus))
matchedControl <- filter(matchedControl , !is.na(Subject.number))
matchedControl <- mutate(matchedControl, MatchStatus = ifelse(FinalDiagnosis == "Cancer", "Cancer", MatchStatus))



removeGroup <- which(colnames(matchedControl) == "Group")
matchedControl <- matchedControl[,-removeGroup]

norms2match <- anti_join(norms, matchedControl, by="Subject.number") %>% select(Subject.number) %>% unlist()

runners <- filter(runnersUp, Subject.number %in% norms2match)
runners <- runners[,colnames(matchedControl)]
runnerAverage <- matchedControl %>% group_by(Gender) %>% summarise(mean(age, na.rm=T), mean(bmi, na.rm=T), 
                                                                   mean(cigsSmoked, na.rm=T), mean(PPI, na.rm=T), mean(hb, na.rm=T))
colnames(runnerAverage)[2:ncol(runnerAverage)] <- colnames(runners[,7:12][,-2])

for (i in 1:(192 - nrow(matchedControl))){
  
  match1 <- dist(rbind(runnerAverage[1,2:ncol(runnerAverage)],runners[,7:12][,-2])) %>% as.matrix() %>% .[1,-1]
  match1ref <- which(match1 == min(match1))
  run <- runners[match1ref,]
  run[,"GlobalMatch"] <- globalCount
  
  matchedControl <- rbind(matchedControl, runners[match1ref,])
  
  runners <- runners[-match1ref,]
  globalCount <- globalCount +1 
}


write.csv(matchedControl, 'Final_Set_of_Matched_Controls.csv', quote=F, row.names=F)

cat("The Cases have all been Matched")


beep(8)