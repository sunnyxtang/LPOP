#############################################

# Compiling LPOP NLP outputs from Winterlight #

#############################################

#This Version: Version 6: 2023 10-26
#New: medium starting set of features

### Plan: 
# Read in data from task-wise data that was already compiled
# By task group, collapse other features
# Remove variables that have no variability
# Remove overly correlated variables (>0.90)
# Run Imputations for missing values
# Run PCA
# Merge with clinical data


######## LOAD PACKAGES ###########

# Data Processing:
library(psych) #psych stats package, includes factor analysis
library(bestNormalize) #yeojohnson and other transformations
library(mice) #multiple imputation
#library(Matching) #propensity score matching
#library(MatchIt) #simpler matching
library(robustbase) #adjusted box plot
library(performance) #check_outliers for mahalanobis distance, influence
#library(NbClust) #Clustering - how many
library(Rtsne) #t-SNE dimensionality reduction
library(missForest)

# General Stats:
library(effsize) #cohen d, etc.
library(lm.beta) #standardized beta
library(DescTools) #Fisher r to z tranformation
library(nlme) #mixed effects
#library(irr) # interrater correlation
library(RVAideMemoire) #Pairwise comparisons and other stats
#library(lme4) #expands on nlme with generalized linear regressions
#library(lavaan) #cFA, mediation analysis
library(moments) #skewness and kurtosis

# Graphing
library(corrplot) # correlation plots

# Standard Packages (all projects)
library(plyr) #data manipulations
library(naniar) #na functions
library(arsenal) #tableby
library(tidyverse) #dplyr, ggplot, stringr, et al


set.seed(123)

load("wll_bytask.R")
load("justvars.R")


######## Acoustic Variables ###########
# Applies to: all tasks

### Identify variables
vars_acoustic <- justvars$feature_name[
      which(justvars$feature.category == "timing" | justvars$feature.category == "acoustic")]
vars_acoustic

# Take out hnr + mfcc
vars_acoustic <- str_subset(vars_acoustic, pattern = "hnr", negate = TRUE)
vars_acoustic <- str_subset(vars_acoustic, pattern = "mfcc_", negate = TRUE)
vars_acoustic <- str_subset(vars_acoustic, pattern = "zcr_", negate = TRUE)

### Initialize empty dataframes by task
table(wll_bytask$task_name)

bytask_PAR <- filter(wll_bytask, task_name == "PAR") %>% select(session_id, all_of(vars_acoustic))

bytask_FLU <- as.data.frame(matrix(nrow = 0, ncol = length(vars_acoustic) + 1))
colnames(bytask_FLU) <- c("session_id", vars_acoustic)

bytask_PIC <- as.data.frame(matrix(nrow = 0, ncol = length(vars_acoustic) + 1))
colnames(bytask_PIC) <- c("session_id", vars_acoustic)

bytask_JOU <- as.data.frame(matrix(nrow = 0, ncol = length(vars_acoustic) + 1))
colnames(bytask_JOU) <- c("session_id", vars_acoustic)

### Loop through each session
lyst.sessionids <- wll_bytask$session_id %>% unique()
n = length(vars_acoustic)

for (i in lyst.sessionids) {
      # Fluency
      temp <- wll_bytask %>% filter(session_id == i & task_name == "FLU") %>% 
            select(session_id, all_of(vars_acoustic)) #per session
      if (nrow(temp) > 0) {
            temp <- colMeans(temp, na.rm = TRUE) #columns means for session
            bytask_FLU[nrow(bytask_FLU) + 1,] <- temp #Add to df
      }
      
      # Picture
      temp <- wll_bytask %>% filter(session_id == i & task_name == "PIC") %>% 
            select(session_id, all_of(vars_acoustic)) #per session
      if (nrow(temp) > 0) {
            temp <- colMeans(temp, na.rm = TRUE) #columns means for session
            bytask_PIC[nrow(bytask_PIC) + 1,] <- temp #Add to df
      }
      
      
      # Journaling
      temp <- wll_bytask %>% filter(session_id == i & task_name == "JOU") %>% 
            select(session_id, all_of(vars_acoustic)) #per session
      if (nrow(temp) > 0) {
            temp <- colMeans(temp, na.rm = TRUE) #columns means for session
            bytask_JOU[nrow(bytask_JOU) + 1,] <- temp #Add to df
      }
      
      
}


### Aggregate features
#Rename
colnames(bytask_PAR) <- c("session_id", paste(vars_acoustic, "PAR", sep = "_"))
colnames(bytask_FLU) <- c("session_id", paste(vars_acoustic, "FLU", sep = "_"))
colnames(bytask_PIC) <- c("session_id", paste(vars_acoustic, "PIC", sep = "_"))
colnames(bytask_JOU) <- c("session_id", paste(vars_acoustic, "JOU", sep = "_"))

#Merge
just_acoustic <- merge(bytask_PAR, bytask_FLU, by="session_id", all = TRUE)
just_acoustic <- merge(just_acoustic, bytask_PIC, by="session_id", all = TRUE)
just_acoustic <- merge(just_acoustic, bytask_JOU, by="session_id", all = TRUE)


head(just_acoustic)



######## Lexical Features ###########
# Applies to fluency, picture, journaling

#Check out categories: lexical, sentiment, 

### Identify variables
vars_lexical <- justvars$feature_name[which(justvars$feature.category=="lexical" | justvars$feature.category=="sentiment")]

vars_add2open <- vars_lexical[c(3:8, 11, 13:17, 19:27, 46:47, 51:104, 106:113)]
vars_lexical <- vars_lexical[which(!vars_lexical %in% vars_add2open)]

vars_lexical

### Initialize empty dataframes by task
bytask_FLU <- as.data.frame(matrix(nrow = 0, ncol = length(vars_lexical) + 1))
colnames(bytask_FLU) <- c("session_id", vars_lexical)

bytask_PIC <- as.data.frame(matrix(nrow = 0, ncol = length(vars_lexical) + 1))
colnames(bytask_PIC) <- c("session_id", vars_lexical)

bytask_JOU <- as.data.frame(matrix(nrow = 0, ncol = length(vars_lexical) + 1))
colnames(bytask_JOU) <- c("session_id", vars_lexical)

### Loop through each session
n = length(vars_lexical)

for (i in lyst.sessionids) {
      # Fluency
      temp <- wll_bytask %>% filter(session_id == i & task_name == "FLU") %>% 
            select(session_id, all_of(vars_lexical)) #per session
      if (nrow(temp) > 0) {
            temp <- colMeans(temp, na.rm = TRUE) #columns means for session
            bytask_FLU[nrow(bytask_FLU) + 1,] <- temp #Add to df
      }
      
      # Picture
      temp <- wll_bytask %>% filter(session_id == i & task_name == "PIC") %>% 
            select(session_id, all_of(vars_lexical)) #per session
      if (nrow(temp) > 0) {
            temp <- colMeans(temp, na.rm = TRUE) #columns means for session
            bytask_PIC[nrow(bytask_PIC) + 1,] <- temp #Add to df
      }
      
      
      # Journaling
      temp <- wll_bytask %>% filter(session_id == i & task_name == "JOU") %>% 
            select(session_id, all_of(vars_lexical)) #per session
      if (nrow(temp) > 0) {
            temp <- colMeans(temp, na.rm = TRUE) #columns means for session
            bytask_JOU[nrow(bytask_JOU) + 1,] <- temp #Add to df
      }
      
      
}


### Aggregate features
#Rename
colnames(bytask_FLU) <- c("session_id", paste(vars_lexical, "FLU", sep = "_"))
colnames(bytask_PIC) <- c("session_id", paste(vars_lexical, "PIC", sep = "_"))
colnames(bytask_JOU) <- c("session_id", paste(vars_lexical, "JOU", sep = "_"))

#Merge
just_lexical <- bytask_FLU
just_lexical <- merge(just_lexical, bytask_PIC, by="session_id", all = TRUE)
just_lexical <- merge(just_lexical, bytask_JOU, by="session_id", all = TRUE)

head(just_lexical)



######## Discourse-level Features ###########
# Applies to picture, journaling
table(justvars$feature.category)
#Check out categories: discourse, gloabl + local coherence, morphological, syntactic

### Identify variables
vars_discourse <- justvars$feature_name[
      which(justvars$feature.category %in% c("discourse", "local coherence",
                                             "syntactic"))
]

vars_discourse <- c(vars_discourse, vars_add2open)

vars_discourse <- str_subset(vars_discourse, pattern = "NOUN_", negate = TRUE)
vars_discourse <- str_subset(vars_discourse, pattern = "VERB_", negate = TRUE)

vars_discourse


# Take out tags, constructions
vars_discourse <- str_subset(vars_discourse, pattern = "tag_", negate = TRUE)


### Initialize empty dataframes by task
bytask_PIC <- as.data.frame(matrix(nrow = 0, ncol = length(vars_discourse) + 1))
colnames(bytask_PIC) <- c("session_id", vars_discourse)

bytask_JOU <- as.data.frame(matrix(nrow = 0, ncol = length(vars_discourse) + 1))
colnames(bytask_JOU) <- c("session_id", vars_discourse)

### Loop through each session
n = length(vars_discourse)

for (i in lyst.sessionids) {
      # Picture
      temp <- wll_bytask %>% filter(session_id == i & task_name == "PIC") %>% 
            select(session_id, all_of(vars_discourse)) #per session
      if (nrow(temp) > 0) {
            temp <- colMeans(temp, na.rm = TRUE) #columns means for session
            bytask_PIC[nrow(bytask_PIC) + 1,] <- temp #Add to df
      }
      
      
      # Journaling
      temp <- wll_bytask %>% filter(session_id == i & task_name == "JOU") %>% 
            select(session_id, all_of(vars_discourse)) #per session
      if (nrow(temp) > 0) {
            temp <- colMeans(temp, na.rm = TRUE) #columns means for session
            bytask_JOU[nrow(bytask_JOU) + 1,] <- temp #Add to df
      }
      
      
}


### Aggregate features
#Rename
colnames(bytask_PIC) <- c("session_id", paste(vars_discourse, "PIC", sep = "_"))
colnames(bytask_JOU) <- c("session_id", paste(vars_discourse, "JOU", sep = "_"))

#Merge
just_discourse <- bytask_PIC
just_discourse <- merge(just_discourse, bytask_JOU, by="session_id", all = TRUE)

head(just_discourse)



######## Fluency task only ###########

### Identify variables
vars_fluency <- c("phonemic_fluency_unique_correct_words_with_pos_filtering", "semantic_fluency_unique_correct_words")
temp <- wll_bytask[which(wll_bytask$task_name=="FLU"), c("task_name", "phonemic_fluency_unique_correct_words_with_pos_filtering", "semantic_fluency_unique_correct_words")]
head(temp) # only defined for 1 task at a time

### Get separate fluency features and merge
just_fluency <- select(wll_bytask, session_id, phonemic_fluency_unique_correct_words_with_pos_filtering) %>%
      filter(!is.na(phonemic_fluency_unique_correct_words_with_pos_filtering))

temp <- select(wll_bytask, session_id, semantic_fluency_unique_correct_words) %>%
      filter(!is.na(semantic_fluency_unique_correct_words))

just_fluency <- merge(just_fluency, temp, by="session_id", all = TRUE)
head(just_fluency)

#Rename
n = ncol(just_fluency)
colnames(just_fluency)[2:n] <- paste(colnames(just_fluency)[2:n], "FLU", sep = "_")

head(just_fluency)



######## Picture task only ###########
table(justvars$feature.category)

### Identify variables
vars_picture <- justvars$feature_name[which(justvars$feature.category=="information content" | 
                                                  justvars$feature.category=="global coherence")]
vars_picture

### Initialize empty dataframe
just_picture <- as.data.frame(matrix(nrow = 0, ncol = length(vars_picture) + 1))
colnames(just_picture) <- c("session_id", vars_picture)

### Loop through each session
n = length(vars_picture)

for (i in lyst.sessionids) {
      temp <- wll_bytask %>% filter(session_id == i & task_name == "PIC") %>% 
            select(session_id, all_of(vars_picture)) #per session
      if (nrow(temp) > 0) {
            temp <- colMeans(temp, na.rm = TRUE) #columns means for session
            just_picture[nrow(just_picture) + 1,] <- temp #Add to df
      }
}



#Rename
n = ncol(just_picture)
colnames(just_picture)[2:n] <- paste(colnames(just_picture)[2:n], "PIC", sep = "_")

head(just_picture)



######## Merge collapsed variables ###########

# Merging #760 raw variables - medium set
wll_collapsed <- just_acoustic
wll_collapsed <- merge(wll_collapsed, just_lexical, by="session_id")
wll_collapsed <- merge(wll_collapsed, just_discourse, by="session_id")
wll_collapsed <- merge(wll_collapsed, just_fluency, by="session_id")
wll_collapsed <- merge(wll_collapsed, just_picture, by="session_id")

colnames(wll_collapsed)

vars_allcollapsed <- names(wll_collapsed[2:dim(wll_collapsed)[2]])
write.csv(vars_allcollapsed, file="lyst.allcollapsedvars_medium.csv")



######## Getting rid of NA/0SD, etc. ###########
wll_raw <- wll_collapsed

### Replace Nan + Inf values with NA
for (feature in vars_allcollapsed) {
      wll_raw[[feature]][which(is.nan(wll_raw[[feature]]))] <- NA
      wll_raw[[feature]][which(is.infinite(wll_raw[[feature]]))] <- NA
}


###Checking standard deviations - all non-zero? removed 18 variables
temp <- lapply(wll_raw, sd, na.rm=TRUE)
wll_raw <- wll_raw[,-which(temp==0 | is.na(temp))]
761-743

# Recheck 
temp <- lapply(wll_raw, sd, na.rm=TRUE)
temp[1:10]

any(temp==0)
any(is.na(temp))


### Remove variables where missing > 20% - none
temp <- is.na(wll_raw) %>% colSums()
summary(temp)
0.2*nrow(wll_raw) # threshold


### Any rows missing all columns? Range 0-569, Med=0, mean =5.87, 3rd Q = 1
temp <- is.na(wll_raw) %>% rowSums()
hist(temp)

which(temp>25)
wll_raw$session_id[which(temp>25)] #REMOVE 2 outliers with a lot of missing features

lplong[which(lplong$wll_sessionid %in% c("18108", "20138")),] # 2 session 1's - both have no followup. Can exclude.

wll_raw <- wll_raw[which(!wll_raw$session_id %in% c("18108", "20138")),]



######## Remove variables with low range in values ###########
# Rationale - trends could be driven by bad data

# Calculate kurtosis
just_kurtosis <- lapply(wll_raw, kurtosis, na.rm=TRUE) %>% unlist()
hist(just_kurtosis)

# Evaluate quantiles
quantile(just_kurtosis, probs = c(.5, .75, 0.9, 0.95, 0.99))

just_kurtosis[which(just_kurtosis>19.63)]

# Cutoff 90% - cut 75 for kurtosis
n = quantile(just_kurtosis, probs = 0.90)
vars_kurtotic <- names(just_kurtosis[which(just_kurtosis>n)])

wll_raw <- select(wll_raw, -any_of(vars_kurtotic))
length(vars_kurtotic)

### Saving the variables at this stage
vars_interim <- colnames(wll_raw)[2:ncol(wll_raw)]



######## Remove covaried variables: short version ###########

# Checking intra-correlations among variables - absolute values
# Strategy: remove those with correlations >0.80 in reverse order of covariation with other variables
# Keep the most overall unique variables

# Reset dataframe
wll_raw <- select(wll_collapsed, session_id, all_of(vars_interim)) %>% 
      filter(!session_id %in% c("18108", "20138"))


### Calculate total correlatedness
# Get matrix
cor.wll_raw <- cor(as.matrix(wll_raw[2:ncol(wll_raw)]), method = "spearman", use = "complete.obs")
cor.wll_raw <- as.data.frame(cor.wll_raw)
cor.wll_raw[1:10, 1:10]

# Abs values
cor.abs.raw <- lapply(cor.wll_raw, abs) %>% as.data.frame()
cor.abs.raw[1:10, 1:10]
rownames(cor.abs.raw) <- colnames(cor.abs.raw)

# Total
just_totalcor <- rowSums(cor.abs.raw, na.rm = TRUE)
head(just_totalcor)

# Order by absolute value of correlation
cor.abs.raw$totalcor <- just_totalcor
cor.abs.raw <- arrange(cor.abs.raw, totalcor)
cor.abs.raw[658:667, 658:668]

# Get variables in order
vars_ordered <- rownames(cor.abs.raw)
vars_ordered <- rev(vars_ordered) # High to Low total correlation
vars_ordered[1:10] # Check!
vars_ordered[658:667]

# Remove totalcor
cor.abs.raw$totalcor <- NULL


### Remove sequentially by total correlation if abs value of correlation is >0.85

vars_pca <- NULL #Initialize the list
temp <- cor.abs.raw #dummy df

for (feature in vars_ordered) {
      temp[feature, feature] <- NA
      
      if (max(temp[[feature]], na.rm=TRUE) > 0.85) {
            temp[[feature]] <- NULL #Get highly correlated features out of the way
      }
      
      else {
            # If not, then add it to the list of variables we are keeping
            vars_pca <- append(vars_pca, feature)
      }
}

vars_pca


######## Checking variable types ###########

# selecting just the uncorrelated ones
wll_pca <- select(wll_raw, session_id, all_of(vars_pca))
wll_pca[1:10, 1:10]


### Checking correlations

# Correlation Matrix - just the variables
dat.pca <- select(wll_pca, all_of(vars_pca))

# Matrix
cor.pca <- cor(dat.pca, method = "spearman", use = "complete.obs")

# Check covariances - all <0.85? - why are there still 4 variables >0.85?
hist(cor.pca)
temp <- round(cor.pca, digits = 2)
table(temp)

for (i in 1:359) {
      for (j in 1:359) {
            if ((i != j) & (temp[i,j] > 0.85)) {
                  print(rownames(temp)[i])
                  print(colnames(temp)[j])
                  print(temp[i,j])
            }
      }
}

# Manually removing
vars_pca <- str_subset(vars_pca, pattern = "Lu_CP.T_JOU", negate = TRUE)
vars_pca <- str_subset(vars_pca, pattern = "NOUN_age_of_acquisition_PIC", negate = TRUE)

# Recalculate
dat.pca <- select(wll_pca, all_of(vars_pca))
cor.pca <- cor(dat.pca, method = "spearman", use = "complete.obs")

# Scale
wll_scraw <- lapply(dat.pca, scale) %>% as.data.frame()
wll_scraw[1:10, 1:10]
wll_scraw$session_id <- wll_raw$session_id

### What kinds of variables are these? Total = 357

# How many derived from each type of task?
select(dat.pca, ends_with("PAR")) %>% ncol() #7 paragraph
select(dat.pca, ends_with("FLU")) %>% ncol() #37 fluency
select(dat.pca, ends_with("PIC")) %>% ncol() #160 picture description
select(dat.pca, ends_with("JOU")) %>% ncol() #153 journaling

# Remove parts at end of variable names to match with dictionary
vars_checktype <- str_replace_all(vars_pca, pattern = "_PAR", replacement = "")
vars_checktype <- str_replace_all(vars_checktype, pattern = "_FLU", replacement = "")
vars_checktype <- str_replace_all(vars_checktype, pattern = "_PIC", replacement = "")
vars_checktype <- str_replace_all(vars_checktype, pattern = "_JOU", replacement = "")

vars_checktype <- as.data.frame(vars_checktype)
colnames(vars_checktype) <- "feature_name"
vars_checktype <- merge(vars_checktype, justvars, all.x=TRUE)

table(vars_checktype$feature.category)
#acoustic 13, discourse 18, global coherence 2, information content 6, lexical 79, local coherence 5,
#morphological 0, picture aggregate 0, sentiment 21, syntactic 177, timing 31, task 2

#From my categorization
vars_checktype$feature_name[which(vars_checktype$feature_name %in% vars_acoustic)] %>% 
      length() #44 acoustic and timing
vars_checktype$feature_name[which(vars_checktype$feature_name %in% vars_lexical)] %>% 
      length() #58 lexical and sentiment
vars_checktype$feature_name[which(vars_checktype$feature_name %in% vars_discourse)] %>% 
      length() #215 discourse organization
vars_checktype$feature_name[which(vars_checktype$feature_name %in% vars_fluency)] %>%
      length() #2 - both fluency task scores
vars_checktype$feature_name[which(vars_checktype$feature_name %in% vars_picture)] %>% 
      length() #8 picture des aggregates

### Save variables going into PCA
vars_checktype$feature_fullname <- vars_pca
write.csv(vars_checktype, file = "variables_for_pca_medium.csv")


######## PCA ###########

### Kaiser-Meyer-Olkin Measure of Sampling Adequacy - "MSA" Measure of sampling adequacy = 0.5. Meh...
cor.pca.kmo <- KMO(cor.pca)
cor.pca.kmo$MSA

### Bartlett Test of Sphericity: - Want <0.05 to indicate significantly different from 0.
cortest.bartlett(dat.pca)

### Determinant  To check for multicolinearity, want >0.00001
det(cor.pca)

#1 components
set.seed(123)
model_1 <- principal(dat.pca, nfactors = 1, rotate = "promax")
model_1$loadings
summary(model_1)

# % Variance explained - 6.7%
model_1$values[1]/length(model_1$values)


######## Imputation ###########
### Which is missing?
# Which variables? - only 12 of them! max missing 6 features
temp <- is.na(dat.pca) %>% colSums() %>% as.vector()
sum(temp>0)
summary(temp)


### Split off vars with missing + 100 random variables

# New df
wll_imp <- wll_scraw
wll_imp$session_id <- NULL

# which variables have missing values?
temp <- which(colSums(is.na(wll_imp)) > 0)
vars_withmissing <- colnames(wll_imp)[temp]
vars_withmissing

# Check:
wll_imp$VERB_imageability_JOU


### Do imputation with RF impute

# impute
wll_imp <- missForest(wll_imp)$ximp

# Re-scale
wll_imp <- scale(wll_imp) %>% as.data.frame()
head(wll_imp)

# Any NAs? - NOPE!
temp <- is.na(wll_imp) %>% colSums() %>% as.vector()
sum(temp>0)

# Check:
wll_imp$VERB_imageability_JOU



######## Component Scores ########### 

### 1 Component model
justmodel1scores <- predict(model_1, wll_imp)
head(justmodel1scores)

justmodel1scores <- scale(justmodel1scores) %>% as.vector()
hist(justmodel1scores)

# Which were top-loaded items?

loadings.model1 <- model_1$loadings[,1] %>% sort(decreasing=TRUE)
loadings.model1[c(1:20, 337:357)]



######## Clean up + save outputs ###########

### Calculated features
out_features.medium <- wll_raw$session_id %>% as.data.frame()
colnames(out_features.medium) <- "session_id"

out_features.medium <- base::cbind(out_features.medium, justmodel1scores)
head(out_features.medium)
colnames(out_features.medium)[2] <- "C.single.medium"
      
save(out_features.medium, file = "output_features_medium.R")

### Imputed raw data

out_fraw.medium <- wll_raw$session_id %>% as.data.frame()
colnames(out_fraw.medium) <- "session_id"

out_fraw.medium <- base::cbind(out_fraw.medium, wll_imp)
head(out_fraw.medium)

save(out_fraw.medium, file = "output_features_raw_medium.R")

### Variables

out_varsmedium <- vars_pca
save(out_varsmedium, file = "output_variablenames_medium.R")

### Loadings
save(loadings.model1, file="loadings.model1.R")
write.csv(loadings.model1, file="loadings.model1.csv")
