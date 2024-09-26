##load packages
require(dplyr)
require(magrittr)
require(ggplot2)
require(ggpubr)
require(survival)
require(survminer)
require(tidyr)
require(gtsummary)
require(data.table)
require(janitor)

##set directory
directory <- "C://Users/U1103123/OneDrive - IQVIA/"

##import all files with pattern
files <- list.files(directory, pattern = "^all_centers_cohort_.*\\.csv$", full.names = TRUE)
alldat <- lapply(files, read.csv)
names(alldat) <- gsub("^all_centers_cohort_|\\.csv$", "", files)
all_data <- rbindlist(alldat,  fill=FALSE, idcol=TRUE)
list2env(all_data, envir = .GlobalEnv)

##add categorical cohort_ids
all_data <- all_data %>% mutate(cohort_id=case_when(grepl("_1029", .id) ~ "mNSCLC all",
                                       grepl("_1030", .id) ~ "Brain mets only",
                                       grepl("_1031", .id) ~ "Liver mets only",
                                       grepl("_1032", .id) ~ "Adrenal gland mets only",
                                       grepl("_1033", .id) ~ "Bone mets only",
                                       grepl("_1034", .id) ~ "Lung mets only",
                                       grepl("_1035", .id) ~ "Other",
                                       grepl("_1036", .id) ~ "Multiple sites with brain mets",
                                       grepl("_1037", .id) ~ "Multiple sites without brain mets",
                                       grepl("_1038", .id) ~ "Mets unknown site",
                                       grepl("_1060", .id) ~ "ECOG 0",
                                       grepl("_1061", .id) ~ "ECOG 1",
                                       grepl("_1062", .id) ~ "ECOG 2",
                                       grepl("_1063", .id) ~ "ECOG 3",
                                       grepl("_1064", .id) ~ "ECOG 4",
                                       grepl("_1080", .id) ~ "Metastatic at presentation",
                                       grepl("_1081", .id) ~ "Non-metastatic at presentation",
                                       TRUE ~ "Invalid_cohort")) %>%
  group_by(cohort_id) %>%
  ungroup()

##reshape data to a format that can be read by survival packages by first uncounting deaths
survival_long <- all_data %>%
  select(TIME_AT_RISK, at_risk, observed, cohort_id) %>%
  mutate(observed2=case_when(observed>=1 ~ 1,
                             TRUE ~ 0)) %>%
  uncount(observed) %>%
  mutate(censored2=0)
  
##second dataset by uncounting censors
survival_long2 <- all_data %>%
  select(TIME_AT_RISK, at_risk, censored, cohort_id) %>%
  mutate(censored2=case_when(censored>=1 ~ 1,
                             TRUE ~ 0)) %>%
  uncount(censored) %>%
  mutate(observed2=0)

##bind censor events and observed events tables together for long format
surv_table <- rbind(survival_long, survival_long2) %>%
  mutate(vital_status=case_when(observed2==1 ~ 1,
                                TRUE ~ 0))



#########OUTPUTS###########

##############Produce overall survival curve - fit model########################

Overall_survival <- survfit(Surv(TIME_AT_RISK*12/365.25, vital_status) ~ 1, subset(surv_table, cohort_id=="mNSCLC all"))

##apply censoring for counts less than 5
##set minimum number at risk at which point KM and risk table is censored
atrisk <- 5

# subset each stratum separately and set survival, upper limit, lower limit and number at risk to missing where <5
maxcutofftime = 0 # for plotting
  cutofftime <- min(Overall_survival$time[Overall_survival$n.risk < atrisk])
  maxcutofftime = max(maxcutofftime, cutofftime)
  cutoffs <- which(Overall_survival$n.risk < atrisk)
  Overall_survival$lower[cutoffs] <- NA
  Overall_survival$upper[cutoffs] <- NA
  Overall_survival$surv[cutoffs] <- NA
  Overall_survival$n.risk[cutoffs] <- NA


# Draw overall curve
Overall_KM <- survminer::ggsurvplot(Overall_survival, risk.table = TRUE, conf.int = TRUE, xlim = c(0,maxcutofftime),
                                    break.x.by = 6, axes.offset = F, legend.lab="mNSCLC all")


##extract risk table object to replace missings with asterisks in risk table to show suppressed counts
tab <- Overall_KM$table
tab$data[is.na(tab$data)] <- "*"

tab$layers = NULL # clear labels
tab <- tab + 
  geom_text(aes(x = time, y = rev(strata), label = llabels))

##reassigns object to Mets_KM
Overall_KM$table <- tab

##calls Mets_KM
print(Overall_KM)  

##survival probability

numbers_overall <- survival:::survmean(Overall_survival, rmean="none")$matrix %>%
  t() %>%
  as.data.frame() %>%
  mutate(`Number censored`=n.start-events) %>%
  rename(`Number of patients`=records) %>%
  rename(`Number of deaths` = events) %>%
  select(`Number of patients`, `Number of deaths`, `Number censored`) %>%
  mutate(across(everything(),  ~ifelse(. < 5, NA, .)))


quantiles_survival_overall <- quantile(Overall_survival) %>%
 
  as.data.frame() %>%
  t() %>%
  as.data.frame() %>%
  filter(row_number() == 1) %>%
  select(`25`, `50`, `75`) %>%
  mutate(across(everything(), round, 2)) %>%
  rename(`Median survival`=`50`)


joined_survival_overall <- cbind(numbers_overall, quantiles_survival_overall) %>%
  mutate(`Q1-Q3`=paste0(`25`, " - ", `75`)) %>%
  
  
  select(`Number of patients`, `Number of deaths`, `Number censored`, `Median survival`, `Q1-Q3`) %>%
  t() %>%
  as.data.frame() %>%
  rename(`mNSCLC all`=quantile) %>%
  mutate(Variable=rownames(.)) 


survival_times_overall <- Overall_survival %>%
  tbl_survfit(times = c(6, 12, 18, 24), label = "mNSCLC all", label_header = "Survival Probability at {time} Month", estimate_fun = ~gtsummary::style_sigfig(.x, digits = 3, scale = 100) ) %>%
  
  as.data.frame() %>%
  
  rename(Variable=`**Characteristic**`) %>%
  
  t() %>%
  as.data.frame() %>%
  mutate(Variable=rownames(.)) %>%
  row_to_names(row_number = 1)
  

patients_at_time_points <- summary(Overall_survival,times= c(6, 12, 18, 24))

cols <- lapply(c(2, 3:5) , function(x) patients_at_time_points[x])
# Combine the columns into a data frame
output_overall <- do.call(data.frame, cols) %>%
  
  
  
  mutate(n=n.risk+n.event+n.censor) %>%
  select(time, n) %>%
  rename(`mNSCLC all`=n) %>%
  rename(Variable=time)

overall_survival_table_output <- rbind(joined_survival_overall, survival_times_overall, output) %>%
  mutate(Variable=case_when(is.na(Variable) ~ rownames(.),
                            Variable=="6" ~ "Number of Patients at 6 Month",
                            Variable=="12" ~ "Number of Patients at 12 Month",
                            Variable=="18" ~ "Number of Patients at 18 Month",
                            Variable=="24" ~ "Number of Patients at 24 Month",
                            TRUE ~ Variable)) %>%
  replace(is.na(.), "*") %>%
  select(Variable, 1) 
rownames(overall_survival_table_output) <- NULL

Overall_survival_table_output <- overall_survival_table_output[c(1:5,10,6,11,7,12,8,13,9), ] 

########Produce site of mets curve##############

##subset data
Mets_survival <- survfit(Surv((TIME_AT_RISK*12/365.25), vital_status) ~ cohort_id, subset(surv_table, cohort_id %in% c(
                                                                                                           "Brain mets only",
                                                                                                           "Liver mets only",
                                                                                                           "Adrenal gland mets only",
                                                                                                           "Bone mets only",
                                                                                                           "Lung mets only",
                                                                                                           "Other",
                                                                                                           "Multiple sites with brain mets",
                                                                                                           "Multiple sites without brain mets",
                                                                                                           "Mets unknown site")))


# subset each stratum separately and set survival, upper limit, lower limit and number at risk to missing where <5
maxcutofftime = 0 # for plotting
strata <- rep(names(Mets_survival$strata), Mets_survival$strata)
for (i in names(Mets_survival$strata)){
  cutofftime <- min(Mets_survival$time[Mets_survival$n.risk < atrisk & strata == i])
  maxcutofftime = max(maxcutofftime, cutofftime)
  cutoffs <- which(Mets_survival$n.risk < atrisk & strata == i)
  Mets_survival$lower[cutoffs] <- NA
  Mets_survival$upper[cutoffs] <- NA
  Mets_survival$surv[cutoffs] <- NA
  Mets_survival$n.risk[cutoffs] <- NA
}
                                                                                                           
##make KM curve
Mets_KM <- ggsurvplot(Mets_survival, risk.table = TRUE, conf.int=TRUE, risk.table.height = 0.5, xlim = c(0,maxcutofftime),
                      break.x.by = 6, axes.offset = F,
                      legend.lab= c("Adrenal gland mets only",
                                                                      "Bone mets only",
                                                                      "Brain mets only",
                                                                      "Liver mets only",
                                                                      "Lung mets only",
                                                                      "Mets unknown site",
                                                                      "Multiple sites with brain mets",
                                                                      "Multiple sites without brain mets",
                                                                      "Other"))

##extract risk table object to replace missings with asterisks in risk table to show suppressed counts
tab <- Mets_KM$table
tab$data[is.na(tab$data)] <- "*"

tab$layers = NULL # clear labels
tab <- tab + 
  geom_text(aes(x = time, y = rev(strata), label = llabels))

##reassigns object to Mets_KM
Mets_KM$table <- tab

##calls Mets_KM
print(Mets_KM)  

##pairwise comparisons mets
p_mets <- pairwise_survdiff(Surv(TIME_AT_RISK, vital_status) ~ cohort_id, subset(surv_table, cohort_id %in% c("Brain mets only",
                                                                                                         "Liver mets only",
                                                                                                         "Adrenal gland mets only",
                                                                                                         "Bone mets only",
                                                                                                         "Lung mets only",
                                                                                                         "Other",
                                                                                                         "Multiple sites with brain mets",
                                                                                                         "Multiple sites without brain mets",
                                                                                                         "Mets unknown site")))
logrank_mets <-p$p.value %>%
  as.data.frame() %>%
  mutate(across(everything(), round, 3))



##survival probability

numbers_mets <- survival:::survmean(Mets_survival, rmean="none")$matrix %>%
  as.data.frame() %>%
  mutate(`Number censored`=n.start-events) %>%
  rename(`Number of patients`=records) %>%
  rename(`Number of deaths` = events) %>%
  select(`Number of patients`, `Number of deaths`, `Number censored`) %>%
  mutate(across(everything(),  ~ifelse(. < 5, NA, .)))
  

quantiles_survival_mets <- quantile(Mets_survival) %>%
  as.data.frame() %>%
  select(quantile.25, quantile.50, quantile.75) %>%
  mutate(across(everything(), round, 2)) %>%
  rename(`Median survival`=quantile.50)


joined_survival_mets <- cbind(numbers_mets, quantiles_survival_mets) %>%
  mutate(`Q1-Q3`=paste0(quantile.25, " - ", quantile.75)) %>%
  mutate(cohort_id=rownames(.)) %>%
  mutate(cohort_id=gsub(".*\\=", "", cohort_id) ) %>%
  
  select(cohort_id, `Number of patients`, `Number of deaths`, `Number censored`, `Median survival`, `Q1-Q3`) %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate(Variable=NA)
  

survival_times <- Mets_survival %>%
  tbl_survfit(times = c(6, 12, 18, 24), label = "Variable", label_header = "Survival Probability at {time} Month", estimate_fun = ~gtsummary::style_sigfig(.x, digits = 3, scale = 100) ) %>%
  
  as.data.frame() %>%
  t() %>%
  
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  
  mutate(Variable=rownames(.)) 

patients_at_time_points <- summary(Mets_survival,times= c(6, 12, 18, 24))

cols <- lapply(c(2, 3:5, 10) , function(x) patients_at_time_points[x])
# Combine the columns into a data frame
output <- do.call(data.frame, cols) %>%


mutate(cohort_id=gsub(".*\\=", "", strata) ) %>%
mutate(n=n.risk+n.event+n.censor) %>%
select(cohort_id, time, n) %>%
  pivot_wider(names_from = cohort_id, values_from = n) %>%
  rename(Variable=time)

Mets_survival_table_output <- rbind(joined_survival_mets, survival_times, output) %>%
  mutate(Variable=case_when(is.na(Variable) ~ rownames(.),
                            Variable=="6" ~ "Number of Patients at 6 Month",
                            Variable=="12" ~ "Number of Patients at 12 Month",
                            Variable=="18" ~ "Number of Patients at 18 Month",
                            Variable=="24" ~ "Number of Patients at 24 Month",
                            TRUE ~ Variable)) %>%
  replace(is.na(.), "*") %>%
  select(Variable, 1:10) 
rownames(Mets_survival_table_output) <- NULL

Mets_survival_table_output <- Mets_survival_table_output[c(1:5,10,6,11,7,12,8,13,9), ]



#######Metastatic at presentation outputs##################

##subset data
Metsyesno_survival <- survfit(Surv(TIME_AT_RISK*12/365.25, vital_status) ~ cohort_id, subset(surv_table, cohort_id %in% c("Metastatic at presentation",
                                                                                                            "Non-metastatic at presentation")))

# subset each stratum separately and set survival, upper limit, lower limit and number at risk to missing where <5
maxcutofftime = 0 # for plotting
strata <- rep(names(Metsyesno_survival$strata), Metsyesno_survival$strata)
for (i in names(Metsyesno_survival$strata)){
  cutofftime <- min(Metsyesno_survival$time[Metsyesno_survival$n.risk < atrisk & strata == i])
  maxcutofftime = max(maxcutofftime, cutofftime)
  cutoffs <- which(Metsyesno_survival$n.risk < atrisk & strata == i)
  Metsyesno_survival$lower[cutoffs] <- NA
  Metsyesno_survival$upper[cutoffs] <- NA
  Metsyesno_survival$surv[cutoffs] <- NA
  Metsyesno_survival$n.risk[cutoffs] <- NA
}


##make KM curve
Metsyesno_KM <- ggsurvplot(Mets_survival, risk.table = TRUE, conf.int=TRUE, xlim = c(0, maxcutofftime),
                           break.x.by = 6, axes.offset = F,
                           legend.lab= c("Metastatic at presentation",
                                         "Non-metastatic at presentation"))


##extract risk table object to replace missings with asterisks in risk table to show suppressed counts
tab <- Metsyesno_KM$table
tab$data[is.na(tab$data)] <- "*"

tab$layers = NULL # clear labels
tab <- tab + 
  geom_text(aes(x = time, y = rev(strata), label = llabels))

##reassigns object to Mets_KM
Metsyesno_KM$table <- tab

##calls Mets_KM
print(Metsyesno_KM)  

##pairwise comparisons mets
p_metsyesno <- pairwise_survdiff(Surv(TIME_AT_RISK, vital_status) ~ cohort_id, subset(surv_table, cohort_id %in% c("Metastatic at presentation",
                                                                                                         "Non-metastatic at presentation")))
logrank_metsyesno <-p$p.value %>%
  as.data.frame() %>%
  mutate(across(everything(), round, 3))

##survival probability

numbers_metsyesno <- survival:::survmean(Metsyesno_survival, rmean="none")$matrix %>%
  as.data.frame() %>%
  mutate(`Number censored`=n.start-events) %>%
  rename(`Number of patients`=records) %>%
  rename(`Number of deaths` = events) %>%
  select(`Number of patients`, `Number of deaths`, `Number censored`) %>%
  mutate(across(everything(),  ~ifelse(. < 5, NA, .)))


quantiles_survival_metsyesno <- quantile(Metsyesno_survival) %>%
  as.data.frame() %>%
  select(quantile.25, quantile.50, quantile.75) %>%
  mutate(across(everything(), round, 2)) %>%
  rename(`Median survival`=quantile.50)


joined_survival_metsyesno <- cbind(numbers_metsyesno, quantiles_survival_metsyesno) %>%
  mutate(`Q1-Q3`=paste0(quantile.25, " - ", quantile.75)) %>%
  mutate(cohort_id=rownames(.)) %>%
  mutate(cohort_id=gsub(".*\\=", "", cohort_id) ) %>%
  
  select(cohort_id, `Number of patients`, `Number of deaths`, `Number censored`, `Median survival`, `Q1-Q3`) %>%
  t() %>%
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  mutate(Variable=NA)


survival_times_metsyesno <- Metsyesno_survival %>%
  tbl_survfit(times = c(6, 12, 18, 24), label = "Variable", label_header = "Survival Probability at {time} Month", estimate_fun = ~gtsummary::style_sigfig(.x, digits = 3, scale = 100) ) %>%
  
  as.data.frame() %>%
  t() %>%
  
  row_to_names(row_number = 1) %>%
  as.data.frame() %>%
  
  mutate(Variable=rownames(.)) 

patients_at_time_points <- summary(Metsyesno_survival,times= c(6, 12, 18, 24))

cols <- lapply(c(2, 3:5, 10) , function(x) patients_at_time_points[x])
# Combine the columns into a data frame
output <- do.call(data.frame, cols) %>%
  
  
  mutate(cohort_id=gsub(".*\\=", "", strata) ) %>%
  mutate(n=n.risk+n.event+n.censor) %>%
  select(cohort_id, time, n) %>%
  pivot_wider(names_from = cohort_id, values_from = n) %>%
  rename(Variable=time)

Metsyesno_survival_table_output <- rbind(joined_survival_metsyesno, survival_times_metsyesno, output) %>%
  mutate(Variable=case_when(is.na(Variable) ~ rownames(.),
                            Variable=="6" ~ "Number of Patients at 6 Month",
                            Variable=="12" ~ "Number of Patients at 12 Month",
                            Variable=="18" ~ "Number of Patients at 18 Month",
                            Variable=="24" ~ "Number of Patients at 24 Month",
                            TRUE ~ Variable)) %>%
  replace(is.na(.), "*") %>%
  select(Variable, 1:10) 
rownames(Metsyesno_survival_table_output) <- NULL

Metsyesno_survival_table_output <- Mets_survival_table_output[c(1:5,10,6,11,7,12,8,13,9), ]










