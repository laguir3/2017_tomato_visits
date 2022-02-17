### DESCRIPTION ####

# Analysis for Aguirre, LA, and Adler, LS. 2022. Interacting Antagonims: Parasite Infection Alters Bombus Impaties (Hymenoptera: Apidae) Responses to herbivory on Tomato Plants. In Review. 

# Code Written by LA Aguirre
# Data Cleaning, Analysis and Visualizations

################## METADATA #
################## counts.csv variables
# bid: Unique bee identifier
# date: Date of observation
# count: Crithidia cell counts per 0.02uL
# ctime: Time of day when Crithidia was counted
# dtime: Time of day when bee was dissected
# wing_raw: Length of marginal cell (right wing) in ocular units at 400x magnification
# 
################## Variables created in script:
# wing_mm: Ocular units converted from wing_raw (13x/20 convertion rate)
# jdate: Julian date from date

################## visits.csv variables
# bid: Unique bee identifier
# start: Plant treatment of first observation
# plant: Plant treatment 
# seconds: Seconds on flower 
# frames: Frames spent on flower (GarageBand application measured recording time in seconds and frames. 25 frames per second)
# 
################## Variables created in script:
# hund: hundredths of second spend on flower from frames
# count: from count in `counts.csv` 
# jdate: Julian date
# # wing_mm: ocular units converted from wing_raw (13x/20 convertion rate)
# status: Infection status (infected/uninfected) from count
# total_visits: number of visits observed for individual bee
# fjdate: Julian date as factor
# fbid: Unique bee identifier as factor

### LOAD #####
# Packages
library(glmmTMB)
library(ggplot2)
library(DHARMa)
library(emmeans)
library(patchwork)
library(tidyverse)
library(car)
library(lmtest) # lrtest 
library(flextable)
library(gdtools) # load fonts
library(grid)

# Set color-blind palette
cb <- c("#000000", 
        "#E69F00", 
        "#56B4E9", 
        "#009E73", 
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7")

# Load crithidia counts data
# setwd("/Users/luis/Documents/Research/2017_tomato/data_files")
counts <- read.csv("tomato_counts.csv", 
                   header = T)

# Cols and rows with data
counts <- counts[1:91, 1:7]

### FORMAT DATA ####
##### Format: Visit Durations ####
# Make new column w/ Julian date
for(i in 1:91){
  if(counts$date[i] == "22-Aug"){
    counts$jdate[i] <- 234
  } else if(counts$date[i] == "23-Aug"){
    counts$jdate[i] <- 235
  } else if(counts$date[i] == "24-Aug"){
    counts$jdate[i] <- 236
  } else if(counts$date[i] == "25-Aug"){
    counts$jdate[i] <- 237
  } else if(counts$date[i] == "26-Aug"){
    counts$jdate[i] <- 238
  }
} 

# Convert ocular units to mm (wing size)
counts$wing_mm <- counts$wing_raw * (13/20)

# Load visitation data
visits <- read.csv("tomato_visits.csv", 
                   header = T)

# Remove the two empty columns that appear in the data frame
visits <- visits[,1:5] 

# Convert frames to hundredths of seconds; output from GarageBand
visits$hund <- (visits$frames * 4) / 100 

# Total time
visits$dur <- visits$seconds + visits$hund 

# Convert the frames that were rounded to twentyfive into 1 second.
for(i in 1:nrow(visits)){
  if(visits$frames[i] == 25){
    visits$frames[i] <- 0
    visits$seconds[i] <- visits$seconds[i] + 1
  }
}

# Check for frames that are equal or greater than 25. 
which(visits$frames >= 25) 
# Two visits with frame values of 30. Their second count is 3, it seems that 
# when they were entered they were incorrectly entered as 30 frames instead of 3
# seconds, and 0 frames
visits$frames[999] <- 3
visits$frames[1102] <- 3 


# Import dates, counts and wing size from counts data frame
for(i in 1:nrow(counts)){
  for(j in 1:nrow(visits)){
    if(counts$bid[i] == visits$bid[j]){
      counts$count[i] -> visits$count[j]
      counts$jdate[i] -> visits$jdate[j]
      counts$wing_mm[i] -> visits$wing_mm[j]
    }
  }
}

# Create a column for infected/uninfected status. Enter 0 as non-infected and 
# 1 as infected
for(i in 1:nrow(visits)){
  if(visits$count[i] >= 1){
    visits$status[i] <- "infected"
  } else if(visits$count[i] == 0){
    visits$status[i] <- "noninfected"
  }
}

# Create new column for number of visits  
# Get individual bee id's
bids <- unique(visits$bid)

# New columns
visits$total_visits <- NA

# Fill column with total visits
for(i in 1:length(bids)){
  total <- nrow(visits[visits$bid == bids[i], ])
  visits$total_visits[visits$bid == bids[i]] <- total
}

# Data types
visits$plant <- as.factor(visits$plant)
visits$status <- as.factor(visits$status)
visits$fjdate <- as.factor(visits$jdate)
visits$fbid <- as.factor(visits$bid)

# Data frame without NA's in count column
visits2 <- na.omit(visits)

# relevel status variabel
# # Set reference levels, non-infected first
visits2$status <- relevel(visits2$status, ref = 2) 
visits$status <- relevel(visits$status, ref = 2)


##### Format: Proportions ####
# Create new data frame
visit_prop <- as.data.frame(matrix(ncol = 6,
                                   nrow = length(unique(visits$bid))))

colnames(visit_prop) <- c("bid", 
                          "visit_prop", 
                          "start_bias", 
                          "count", 
                          "jdate",
                          "total_visits")

# Get bee identifiers
visit_prop$bid <- unique(visits$bid)

# Get counts and date from counts df
for(i in 1:nrow(visit_prop)){
  for(j in 1:nrow(counts)){
    if(visit_prop$bid[i] == counts$bid[j]){
      visit_prop$count[i] <- counts$count[j]
      visit_prop$jdate[i] <- counts$jdate[j]
    }
  }
}

# Get visit proportions and start bias
for(i in 1:nrow(visit_prop)){
  
  # Subset data by bee
  temp <- subset(visits, 
                 bid == visit_prop$bid[i], 
                 select = c(start, 
                            plant, 
                            total_visits))
  
  # Divide number of visits to damage by all visits observed
  visit_prop$visit_prop[i] <- length(temp$plant[temp$plant == 'damage']) / nrow(temp) 
  visit_prop$start_bias[i] <- as.character(temp$start[1])
  
  # Get total visits
  visit_prop$total_visits[i] <- temp$total_visits[1]
}

# Data types
visit_prop$start_bias <- as.factor(visit_prop$start_bias)
visit_prop$jdate <- as.factor(visit_prop$jdate)


##### Format: Sequence Lengths ####
# Create new data frame
visit_seq <- as.data.frame(matrix(ncol = 3, 
                                  nrow = 1))

colnames(visit_seq) <- c("bid", 
                         "plant", 
                         "seq_length")

# Get bee identifiers
bids <- unique(visits$bid)

# For loop to calculate length of consecutive visit sequences 
# to same plant treatment
for(i in 1:length(bids)){
  
  # Subset for each bee
  temp <- subset(visits,
                 bid == bids[i], 
                 select = (1:3))
  
  # temporary data frame 
  temp_values <- as.data.frame(matrix(ncol = 3, 
                                      nrow = 1))
  
  # Consecutive visits counter
  for(j in 1:nrow(temp)){
    
    if(j == 1){ # Initialize
      
      # counter
      temp_seq <- 1 
      
      # place holder for sequences
      temp_values <- c(temp$bid[j], 
                       paste(temp$plant[j]), 
                       temp_seq)
      
    }  
    # If current visit is NOT to same plant treatment as previous visit,
    # begin new count for next plant treatmen
    else if(j >= 2 && 
              j < nrow(temp) && 
              temp$plant[j] != temp$plant[j-1]){
      
      
      visit_seq <- rbind(visit_seq, temp_values)
      
      # reset counter
      temp_seq <- 1
      
      # place holder
      temp_values <- c(temp$bid[j], 
                       paste(temp$plant[j]), 
                       temp_seq)
      
    } 
    # If current visit is to same plant treatment as previous, continue count
    else if(j >= 2 && 
              j < nrow(temp) && 
              temp$plant[j] == temp$plant[j-1]){
      
      # Counter
      temp_seq = temp_seq + 1
      
      # place holder
      temp_values <- c(temp$bid[j], 
                       paste(temp$plant[j]), 
                       temp_seq)
      
    } 
    # If last visit is NOT to same plant as previous, set last sequence to 1
    else if(j == nrow(temp) && 
            temp$plant[j] != temp$plant[j-1]){
      
      visit_seq <- rbind(visit_seq, 
                         temp_values)
      
      # reset count
      temp_seq <- 1
      
      # place holder
      temp_values <- c(temp$bid[j], 
                       paste(temp$plant[j]), 
                       temp_seq)
      
      # bind place holder onto data frame
      visit_seq <- rbind(visit_seq, temp_values)
      
    } 
    # If last visit is to same plant treatment as previous, declare count total
    else if(j == nrow(temp) && 
            temp$plant[j] == temp$plant[j-1]){
      
      # continue counter
      temp_seq <- temp_seq + 1
      
      # place holder
      temp_values <- c(temp$bid[j], 
                       paste(temp$plant[j]), 
                       temp_seq)
      
      # bind place holder onto data frame
      visit_seq <- rbind(visit_seq, temp_values)
      
    }
  }
}

# Remove NA
visit_seq <- visit_seq[-1,]

# New cols for counts and date
visit_seq$count <- NA; visit_seq$jdate <- NA

# Data type
visit_seq[,3] <- as.integer(visit_seq[,3])


# Import dates and counts from counts data frame to visits data frame
for(i in 1:nrow(counts)){
  for(j in 1:nrow(visit_seq)){
    if(counts$bid[i] == visit_seq$bid[j]){
      counts$count[i] -> visit_seq$count[j]
      counts$jdate[i] -> visit_seq$jdate[j]
    }
  }
}

# Data types
visit_seq$jdate <- as.factor(paste((visit_seq$jdate)))
visit_seq$plant <- as.factor(visit_seq$plant)
visit_seq$bid <- as.factor(visit_seq$bid)
visit_seq <- na.omit(visit_seq)

### DATA ANALYSIS ####
##### Pearson Test: Count and Wing Size #### 

# Pearson's correlation test
cor.test(~ count + wing_mm,
         data = visits2) 

# t = -10.522, df = 2361, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   -0.2498377 -0.1727979
# sample estimates:
#   cor 
# -0.2116465 

##### Model 1: Visit Duration ~ Count * Treatment ####
###### Model selection #### 
### Global model
dur_mod1 <- glmmTMB(dur ~ plant + count + fjdate 
                    + wing_mm+
                      plant*count*fjdate +
                      (1|bid),
                    data = visits2, 
                    family = Gamma(link = "log"),
                    contrasts = list(plant = "contr.sum",
                                     fjdate = "contr.sum"))

# Check
summary(dur_mod1)
simulateResiduals(fittedModel = dur_mod1, plot = T)

# Wing size, not significant drop, and use full data set (i.e. visits df)
dur_mod2 <- glmmTMB(dur ~ plant + count + fjdate +
                      plant*count*fjdate +
                      (1|bid),
                    data = visits, 
                    family = Gamma(link = "log"),
                    contrasts = list(plant = "contr.sum", 
                                     fjdate = "contr.sum")) 
# summary(dur_mod2)

# Remove triple interaction
dur_mod3 <- update(dur_mod2, .~. - plant:count:fjdate)

# Compare
AIC(dur_mod2, # best
    dur_mod3) 

# Remove interactions 
dur_mod4.1 <- update(dur_mod3, .~. -plant:count)
dur_mod4.2 <- update(dur_mod3, .~. -plant:fjdate)
dur_mod4.3 <- update(dur_mod3, .~. -count:fjdate)

# Compare
AIC(dur_mod3, 
    dur_mod4.1, 
    dur_mod4.2, 
    dur_mod4.3) # best

# Remove next interactions
dur_mod5.1 <- update(dur_mod4.3, .~. -plant:count)
dur_mod5.2 <- update(dur_mod4.3, .~. -plant:fjdate)

# Compare
AIC(dur_mod4.3, # best
    dur_mod5.1, 
    dur_mod5.2) 

# Model diagnostics
simulateResiduals(fittedModel = dur_mod4.3, 
                  plot = T) # KS-test significant

# Final model 
summary(dur_mod4.3) 
# Family: Gamma  ( log )
# Formula:          
#   dur ~ plant + count + fjdate + (1 | bid) + plant:count + plant:fjdate
# Data: visits
# 
# AIC      BIC   logLik deviance df.resid 
# 10182.3  10263.3  -5077.1  10154.3     2396 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# bid    (Intercept) 0.07661  0.2768  
# Number of obs: 2410, groups:  bid, 90
# 
# Dispersion estimate for Gamma family (sigma^2): 0.319 
# 
# Conditional model:
#                  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)     1.3906712  0.0373464   37.24  < 2e-16 ***
# plant1         -0.0140084  0.0144859   -0.97  0.33352    
# count          -0.0007757  0.0009930   -0.78  0.43475    
# fjdate1         0.1964987  0.0633490    3.10  0.00192 ** 
# fjdate2        -0.1598133  0.0602519   -2.65  0.00799 ** 
# fjdate3         0.0218406  0.0648766    0.34  0.73638    
# fjdate4        -0.0840694  0.0607286   -1.38  0.16625    
# plant1:count    0.0007654  0.0003707    2.06  0.03898 *  
# plant1:fjdate1 -0.0675177  0.0237474   -2.84  0.00447 ** 
# plant1:fjdate2 -0.0180252  0.0234622   -0.77  0.44233    
# plant1:fjdate3  0.0015471  0.0248756    0.06  0.95041    
# plant1:fjdate4  0.0404432  0.0236309    1.71  0.08700 .  


###### Post hoc: Wald test ####
(tab <- Anova(dur_mod4.3, type = "III")) # Save results
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: dur
#                  Chisq Df Pr(>Chisq)    
# (Intercept)  1386.5978  1  < 2.2e-16 ***
# plant           0.9352  1   0.333524    
# count           0.6101  1   0.434748    
# fjdate         15.7065  4   0.003439 ** 
# plant:count     4.2619  1   0.038975 *  
# plant:fjdate   10.8175  4   0.028693 * 

# Double check interaction terms
lrtest(dur_mod4.3, .~. -plant:count) # close
#   #Df  LogLik Df  Chisq Pr(>Chisq)  
# 1  14 -5077.1                       
# 2  13 -5079.3 -1 4.2672    0.03886 *

lrtest(dur_mod4.3, .~. -plant:fjdate) # close
#   #Df  LogLik Df  Chisq Pr(>Chisq)  
# 1  14 -5077.1                       
# 2  10 -5082.5 -4 10.778    0.02917 *


##### Model 2: Visit Proportions ~ Count ##### 
###### Model selection ####
# Global Model
prop_mod1 <- glmmTMB(visit_prop ~ count * start_bias * jdate, 
                     data = visit_prop, 
                     family = gaussian(),
                     contrasts = list(start_bias = "contr.sum",
                                      jdate = "contr.sum"))

# Model diagnostics
prop_mod_check <- simulateResiduals(fittedModel = prop_mod1, 
                                    plot = T)

# Remove triple interaction
prop_mod2 <- update(prop_mod1, .~. -count:start_bias:jdate) 

# Compare
AIC(prop_mod1, 
    prop_mod2) # best lower AIC

# Remove interactions
prop_mod3.1 <- update(prop_mod2, .~. -count:start_bias)
prop_mod3.2 <- update(prop_mod2, .~. -start_bias:jdate)
prop_mod3.3 <- update(prop_mod2, .~. -count:jdate) 

# Compare
AIC(prop_mod2, 
    prop_mod3.1, 
    prop_mod3.2, 
    prop_mod3.3) # best lower AIC

# Remove next interaction
prop_mod4.1 <- update(prop_mod3.3, .~. -count:start_bias) 
prop_mod4.2 <- update(prop_mod3.3, .~. -start_bias:jdate)

# Compare
AIC(prop_mod3.3, 
    prop_mod4.1, # best fewer parameter
    prop_mod4.2)

# Remove last interaction
prop_mod5 <- update(prop_mod4.1, .~. -start_bias:jdate) 

# Compare
AIC(prop_mod4.1, 
    prop_mod5) # best fewer parameters

# Remove single terms 
prop_mod6.1 <- update(prop_mod5, .~. -count)
prop_mod6.2 <- update(prop_mod5, .~. -start_bias)
prop_mod6.3 <- update(prop_mod5, .~. -jdate) 

# Compare
AIC(prop_mod5, 
    prop_mod6.1, 
    prop_mod6.2, 
    prop_mod6.3) # best

# Remove remaining terms
prop_mod7.1 <- update(prop_mod6.3, .~. -count) 
prop_mod7.2 <- update(prop_mod6.3, .~. -start_bias)

# Compare
AIC(prop_mod6.3, 
    prop_mod7.1, # best
    prop_mod7.2)

# Remove last terms
prop_mod8 <- update(prop_mod7.1, .~. -start_bias)

# Compare
AIC(prop_mod7.1, 
    prop_mod8) # best model, fewer parameters

# Final Model 
summary(prop_mod8)

##### Model 3: Sequence Lengths ~ Count * Treatment #####
###### Model Selection ####
# Global Model
seq_mod_norm <- glmmTMB(seq_length ~ plant * count * jdate + 
                          (1|bid),
                        data = visit_seq,
                        family = gaussian(),
                        contrasts = list(plant = "contr.sum",
                                         jdate = "contr.sum"))

seq_mod_pois  <- glmmTMB(seq_length ~ plant * count * jdate + (1|bid), 
                         data = visit_seq,
                         family = poisson(link = "log"),
                         contrasts = list(plant = "contr.sum",
                                          jdate = "contr.sum"))

# NOTE seq_mod_pois fitting spits out warning but troubleshooting vignette 
# suggest this is not a problem

# Compare
AIC(seq_mod_norm, 
    seq_mod_pois)

# Rename
seq_mod1 <- seq_mod_pois

# Remove Random Effect
seq_mod2 <- update(seq_mod1, .~. - (1|bid)); 

# Compare
AIC(seq_mod2, 
    seq_mod1) # best

# Remove triple interaction
seq_mod3 <- update(seq_mod1, .~. -plant:count:jdate) 

# Compare
AIC(seq_mod1, 
    seq_mod3) # Best

# Remove interactions
seq_mod4.1 <- update(seq_mod3, .~. -plant:count)
seq_mod4.2 <- update(seq_mod3, .~. -plant:jdate)
seq_mod4.3 <- update(seq_mod3, .~. -count:jdate) 

# Compare
AIC(seq_mod3, 
    seq_mod4.1, 
    seq_mod4.2, 
    seq_mod4.3) # best

# Remove interactions
seq_mod5.1 <- update(seq_mod4.3, .~. -plant:count) # best
seq_mod5.2 <- update(seq_mod4.3, .~. -plant:jdate)

# Compare
AIC(seq_mod4.3, 
    seq_mod5.1, 
    seq_mod5.2)

# remove last interaction
seq_mod6 <- update(seq_mod5.1, .~. -plant:jdate) 

# Compare
AIC(seq_mod5.1, # best
    seq_mod6)

# remove non-interactive term
seq_mod7 <- update(seq_mod5.1, .~. -count)

# Compare
AIC(seq_mod5.1, 
    seq_mod7) # best

# Final Model 
summary(seq_mod7)

# Wald Test
Anova(seq_mod7, type = 'III')
# Response: seq_length
#                Chisq Df Pr(>Chisq)    
# (Intercept) 614.6033  1    < 2e-16 ***
# plant         0.3495  1    0.55440    
# jdate         1.2181  4    0.87511    
# plant:jdate   9.7571  4    0.04472 *  

##### Alternate analysis with infection status ####

# NOTE: THIS MODEL INCLUDES A DICHOTOMIZED VARIABLE FOR INFECTION STATUS.
# USING DICHOTOMIZED VARIABLES DERIVED FROM CONTINUOUS VARIABLES CAN BE 
# HIGHLY PROBLEMATIC

###### Model selection ####
# Global model
status_mod1 <- glmmTMB(dur ~ plant + status + fjdate + wing_mm +
                         plant*status*fjdate +
                         (1|bid),
                       contrasts = list(plant = "contr.sum",
                                        status = "contr.sum",
                                        fjdate = "contr.sum"),
                       data = visits2, 
                       family = Gamma(link = "log"))
# Check
summary(status_mod1)
simulateResiduals(fittedModel = status_mod1, plot = T)

# Drop wing size, use full data sel
status_mod2 <- glmmTMB(dur ~ plant + status + fjdate 
                       + plant*status*fjdate +
                         (1|bid),
                       contrasts = list(plant = "contr.sum",
                                        status = "contr.sum",
                                        fjdate = "contr.sum"),
                       data = visits, 
                       family = Gamma(link = "log"))
#summary(status_mod2)

# Remove triple interaction
status_mod3 <- update(status_mod2, .~. - plant:status:fjdate)

# Compare
AIC(status_mod2, 
    status_mod3) # best

# Remove interactions 
status_mod4.1 <- update(status_mod3, .~. -plant:status)
status_mod4.2 <- update(status_mod3, .~. -plant:fjdate)
status_mod4.3 <- update(status_mod3, .~. -status:fjdate)

# Compare
AIC(status_mod3, # best
    status_mod4.1, 
    status_mod4.2, 
    status_mod4.3) 

# Model diagnostics
simulateResiduals(fittedModel = status_mod3, 
                  plot = T)

# Final Model
summary(status_mod3)
# Family: Gamma  ( log )
# Formula:          
#   dur ~ plant + status + fjdate + (1 | bid) + plant:status + plant:fjdate +  
#   status:fjdate
# Data: visits
# 
# AIC      BIC   logLik deviance df.resid 
# 10177.4  10281.6  -5070.7  10141.4     2392 
# 
# Random effects:
#   
#   Conditional model:
#   Groups Name        Variance Std.Dev.
# bid    (Intercept) 0.06665  0.2582  
# Number of obs: 2410, groups:  bid, 90
# 
# Dispersion estimate for Gamma family (sigma^2): 0.319 
# 
# Conditional model:
#                  Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      1.405235   0.032430   43.33  < 2e-16 ***
# plant1          -0.001076   0.012735   -0.08  0.93267    
# status1         -0.004474   0.032421   -0.14  0.89024    
# fjdate1          0.281573   0.069461    4.05 5.04e-05 ***
# fjdate2         -0.186678   0.057777   -3.23  0.00123 ** 
# fjdate3         -0.003239   0.062055   -0.05  0.95837    
# fjdate4         -0.091999   0.058267   -1.58  0.11435    
# plant1:status1  -0.032229   0.013008   -2.48  0.01322 *  
# plant1:fjdate1  -0.073284   0.024084   -3.04  0.00234 ** 
# plant1:fjdate2  -0.015574   0.023479   -0.66  0.50712    
# plant1:fjdate3   0.004140   0.024850    0.17  0.86768    
# plant1:fjdate4   0.043673   0.023725    1.84  0.06565 .  
# status1:fjdate1  0.210495   0.069459    3.03  0.00244 ** 
# status1:fjdate2 -0.087782   0.057745   -1.52  0.12847    
# status1:fjdate3 -0.064629   0.062054   -1.04  0.29765    
# status1:fjdate4 -0.082055   0.058280   -1.41  0.15915 

###### Post hoc: Wald test ####
Anova(status_mod3, type = "III")
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: dur
#                   Chisq Df Pr(>Chisq)    
# (Intercept)   1877.5968  1  < 2.2e-16 ***
# plant            0.0071  1  0.9326666    
# status           0.0190  1  0.8902420    
# fjdate          23.1320  4  0.0001192 ***
# plant:status     6.1391  1  0.0132226 *  
# plant:fjdate    12.0183  4  0.0172158 *  
# status:fjdate   11.3210  4  0.0231837 *

### FIGURES AND TABLES ####

##### Figure 1: Visit Duration ~ Count * Treatment ####
###### Model estimates ######
# Get trends
emmod <- emtrends(dur_mod4.3, 
                  pairwise ~ plant, 
                  var = "count", 
                  transform = "response")

# Create df for estimates
mylist <- list(count = seq(min(visits$count),
                           max(visits$count),
                           by = 1), 
               plant = c("control", "damage"))

# Get estimates for figure
est <- emmip(dur_mod4.3,
             plant ~ count,
             at = mylist,
             CIs = TRUE,
             type = "response", 
             style = "numeric",
             plotit = FALSE)

###### Main plot ######
pdata_main <- ggplot(data = est,
                     aes(x = count, 
                         y = yvar,
                         color = plant)) +
  geom_line(size = .75) + 
  geom_point(data = visits, 
             aes(x = count,
                 y = dur,
                 shape = plant,
                 color = plant),
             size = 1,
             alpha = 0.25,
             stroke = .8) +
  theme_classic() +
  xlab(expression(paste(italic("Crithidia"), " Cells/0.02 ", mu,"L"))) + 
  ylab("Visit Duration (s)") +
  scale_shape_manual(name = "Plant Treatment", 
                     values = c(1, 
                                2),
                     label = c("Control",
                               "Herbivory")) +
  scale_color_manual(name = "Plant Treatment",
                     values = c(cb[1],
                                cb[7]),
                     label = c("Control", "Herbivory")) +
  scale_fill_manual(name = "Plant Treatment", 
                    values = c(cb[1],
                               cb[7]),
                    label = c("Control",
                              "Herbivory")) +
  theme(legend.position = "bottom",
        legend.key.height = unit(7, "mm")) +
  guides(shape = guide_legend(override.aes = list(alpha = 1)))
pdate_main

###### Trends inset####
ptrend_ins <- ggplot(data = est, 
                     aes(x = count,
                         y = yvar,
                         color = plant)) + 
  geom_line(size = 1.5) +
  theme_classic() +
  geom_ribbon(aes(ymax = UCL,
                  ymin = LCL,
                  fill = plant),
              alpha = 0.1) +
  xlab(expression(paste(italic("Crithidia"),
                        " Cells/0.02 ", 
                        mu,
                        "L"))) + 
  ylab("Visit Duration (s)") +
  scale_color_manual(name = "Plant Treatment", 
                     values = c(cb[1],
                                cb[7]),
                     label = c("Control",
                               "Herbivory")) +
  scale_fill_manual(name = "Plant Treatment", 
                    values = c(cb[1],
                               cb[7]),
                    label = c("Control",
                              "Herbivory")) +
  theme(legend.position = "none", 
        axis.title = element_text(size = 9)) + 
  scale_x_continuous(breaks = c(0, 100, 200))

ptrend_ins

# Plot with data and trends inset
tom_plot <- pdata_main + inset_element(ptrend_ins,
                                       right = 1,
                                       top = 1,
                                       left = 0.47,
                                       bottom = 0.5)

# Alternative plot, inversed
tom_plot2 <- ptrend_ins + inset_element(pdata_main, 
                                        right = 1,
                                        top = 1,
                                        left = 0.47,
                                        bottom = 0.5)

# Save figure w/ publication qualities
ggsave("Figure_1.pdf",
       tom_plot,
       units = "mm", 
       width = 82,
       height = 125,
       dpi = 1200)

##### Figure S1: Visit Proportions ~ Count ####
figureS1 <- ggplot(visit_prop) + 
  geom_point(aes(y = visit_prop, 
                 x = count)) +
  theme_classic() +
  labs(x = expression(paste(italic("Crithidia "),
                            "Cells/0.02",
                            mu,
                            "l")),
       y = "Proportion of Visits to Flowers on Damaged Plants")

# Save figure
ggsave("Figure_S1.pdf",
       figureS1,
       width = 5, 
       height = 5)

##### Figure S2: Sequence Lengths ~ Count * Treatment ####
# Labels for plot 
labely <- "Number of Consecutive Visits to Flowers\n from the Same Treatment"
stats_ann <- "Plant Treatment: p-value = 0.465 \n Date: p-value = 0.660 \n Treatment x Date:  p-value = 0.045"

# Plot
figureS2 <- ggplot(visits, aes(x = fjdate, 
                               y = dur, 
                               fill = plant), 
                   alpha = .75) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(x = "Julian Date",  # Labels
       y = labely, 
       fill = "Treatment") +
  scale_fill_manual(name = "Plant Treatment", # Set colors
                    values = c(cb[1],
                               cb[7]),
                    label = c("Control",
                              "Herbivory")) + 
  annotate("text", 
           x = 4, 
           y = 30, 
           label = stats_ann, 
           size = 4, 
           color = "black") + 
  theme(legend.position = "bottom")

# Save figure
ggsave("Figure_S2.pdf",
       figureS2,
       width = 5,
       height = 5)

##### Figure S3: Visit Duration ~ Count * Treatment * Date ######
###### Model estimates ####
emmod_date <- emtrends(dur_mod4.3, 
                       pairwise ~ plant | fjdate, 
                       var = "count", 
                       transform = "response",
                       mult.name = "fjdate")

mylist_date <- list(count = seq(min(visits$count),
                                max(visits$count),
                                by = 1), 
                    plant = c("control", 
                              "damage"),
                    fjdate = levels(visits$fjdate))

###### Plot #####
est_date <- emmip(dur_mod4.3,
                  plant ~ count | fjdate,
                  mult.name = "fjdate",
                  at = mylist,
                  CIs = TRUE,
                  type = "response", 
                  style = "numeric",
                  plotit = F) 

figureS3 <- ggplot(est_date,
                   (aes(x = count, 
                        y = yvar,
                        color = plant))) + 
  geom_line(size = 1.5) + 
  facet_wrap(~fjdate) +
  theme_classic() +
  geom_ribbon(aes(ymax = UCL,
                  ymin = LCL,
                  fill = plant),
              alpha = 0.2) +
  xlab(expression(paste(italic("Crithidia"),
                        " Cells/0.02 ", 
                        mu,
                        "L"))) + 
  ylab("Visit Duration (s)") +
  scale_color_manual(name = "Plant Treatment", 
                     values = c(cb[1],
                                cb[7]),
                     label = c("Control",
                               "Herbivory")) +
  scale_fill_manual(name = "Plant Treatment", 
                    values = c(cb[1],
                               cb[7]),
                    label = c("Control",
                              "Herbivory")) +
  theme(legend.position = "bottom", 
        axis.title = element_text(size = 9)) + 
  scale_x_continuous(breaks = c(0, 100, 200))


# Save plot
ggsave("Figure_S3.pdf",
       figureS3,
       units = "mm", 
       width = 125,
       height = 125,
       dpi = 1200)

##### Table 1 #########
###### Format data frame #### 
tab <- as.data.frame(tab)

# Edit data frame
# Round
tab <- round(tab, digits = 3)

# Include row names
tab <- rownames_to_column(tab, "row_names")

# Change row names
tab[,1] <- c("Intercept",
             "Herbivory",
             "Crithidia Count",
             "Julian Date",
             "Herbivory x Crithidia Count",
             "Herbivory x Julian Date")

# Change column names 
names(tab)[1] <- "Fixed Effects"


###### Create flextable ####
myft <- qflextable(tab)

# Header style
myft <- italic(myft, 
               part = "header")
myft <- bold(myft, 
             part = "header")

# Change colnames
chi <- "\u03C7"
myft <- set_header_labels(myft,
                          `Pr(>Chisq)` = "p-value",
                          `Chisq` = paste0("\u03C7","\U00B2"))

# Bold significant
myft <- bold(myft, 
             i = c(1, 4, 5, 6),
             part = "body")

# Make Times New Roman
myft <- font(myft, 
             part = "all",
             fontname = "Times New Roman")

### Save table
save_as_docx("Table 1" = myft, 
             path = "Table1.docx")

