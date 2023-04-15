library(survival)
library(survminer)

##############################################################
# load data
setwd('S://0417_Antidepressants - Association of long-term Antidepressants use in youth with t//Whitehall data')
s1quest = read.csv('s1quest_0417_S.csv',header=TRUE)
s2quest = read.csv('s2quest_0417_S.csv',header=TRUE)
s3quest = read.csv('s3quest_0417_S.csv',header=TRUE)
s4quest = read.csv('s4quest_0417_S.csv',header=TRUE)
s5quest = read.csv('s5quest_0417_S.csv',header=TRUE)
s6quest = read.csv('s6quest_0417_S.csv',header=TRUE)
s7quest = read.csv('s7quest_0417_S.csv',header=TRUE)
s8quest = read.csv('s8quest_0417_S.csv',header=TRUE)
s9quest = read.csv('s9quest_0417_S.csv',header=TRUE)
s11quest = read.csv('s11quest_0417_S.csv',header=TRUE)
s12quest = read.csv('s12quest_0417_S.csv',header=TRUE)
s1screen = read.csv('s1screen_0417_S.csv',header=TRUE)
s3screen = read.csv('s3screen_0417_S.csv',header=TRUE)
s5screen = read.csv('s5screen_0417_S.csv',header=TRUE)
s7screen = read.csv('s7screen_0417_S.csv',header=TRUE)
s9screen = read.csv('s9screen_0417_S.csv',header=TRUE)
s11screen = read.csv('s11screen_0417_S.csv',header=TRUE)
s12screen = read.csv('s12screen_0417_S.csv',header=TRUE)

# reform data
s1 = merge(s1quest,s1screen)
s2 = s2quest
s3 = merge(s3quest,s3screen)
s4 = s4quest
s5 = merge(s5quest,s5screen)
s6 = s6quest
s7 = merge(s7quest,s7screen)
s8 = s8quest
s9 = merge(s9quest,s9screen)
s11 = merge(s11quest,s11screen)
s12 = merge(s12quest,s12screen)

data = merge(s1,s2)
data = merge(data,s3)
data = merge(data,s4)
data = merge(data,s5)
data = merge(data,s6)
data = merge(data,s7)
data = merge(data,s8)
data = merge(data,s9)
data = merge(data,s11)
data = merge(data,s12)

#write.csv(data,'reformed data.csv')


#############################################################
###survival analysis

### 0. group patient by times of CNSDRG
data_CNS = data[c(grep('CNSDRG',colnames(data)))]
for (i in seq(10308)) {
  data_CNS$sum[i] = sum(data_CNS[i,1:6],na.rm=TRUE)
}
hist(data_CNS$sum)
sum(data_CNS$sum>1)
CNS_group = data_CNS$sum>1

### 1. select coveriates
sex = data$SEX
age = data$AGE_Q_5
marriage = data$STATUS
educ = data$EDUC
ethic = data$ETHN_DS
lsid = data$LILLC04
coveriate = cbind(age,sex,marriage,educ,ethic,lsid)
summary(coveriate)

### 2. breathing related illness---------------------------

# 2.1 extract needed variables
data_athma = data[c(grep('LSI_ASTHMA',colnames(data)))]
data_CHD = data[c(grep('LSI_CHD',colnames(data)))]
data_COPD = data[c(grep('LSI_COPD',colnames(data)))]
breath_factor = data$LILLC09
breath_data = cbind(breath_factor,data_athma,data_CHD,data_COPD,data_CNS$sum)

# 2.2 athma survival analysis

# reform data and add coveriates
data_athma = cbind(breath_factor,CNS_group,data_athma,sex,age,marriage,ethic,educ,lsid)
data_athma = subset(data_athma,data_athma$breath_factor==0)
times = c(1,2,3,4,5,6)
data_athma$time = 6
data_athma$status = 0
for (i in seq(7231)){
  for (j in 3:8){
    if (is.na(data_athma[i,j])==FALSE){
      data_athma$time[i] = times[j-2]
      data_athma$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_athma)
ggsurvplot(surv_model,data = data_athma,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_athma)
summary(cox_model)
ggforest(cox_model,data=data_athma)


# 2.3 CHD survival analysis

# reform data and add coveriates
data_CHD = cbind(breath_factor,CNS_group,data_CHD,sex,age,marriage,ethic,educ,lsid)
data_CHD = subset(data_CHD,data_CHD$breath_factor==0)
times = c(1,2,3,4,5,6)
data_CHD$time = 6
data_CHD$status = 0
for (i in seq(7231)){
  for (j in 3:8){
    if (is.na(data_CHD[i,j])==FALSE){
      data_CHD$time[i] = times[j-2]
      data_CHD$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_CHD)
ggsurvplot(surv_model,data = data_CHD,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_CHD)
summary(cox_model)
ggforest(cox_model,data=data_CHD)



# 2.3 COPD survival analysis

# reform data and add coveriates
data_COPD = cbind(breath_factor,CNS_group,data_COPD,sex,age,marriage,ethic,educ,lsid)
data_COPD = subset(data_COPD,data_COPD$breath_factor==0)
times = c(1,2,3,4,5,6)
data_COPD$time = 6
data_COPD$status = 0
for (i in seq(7231)){
  for (j in 3:8){
    if (is.na(data_COPD[i,j])==FALSE){
      data_COPD$time[i] = times[j-2]
      data_COPD$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_COPD)
ggsurvplot(surv_model,data = data_COPD,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_COPD)
summary(cox_model)
ggforest(cox_model,data=data_COPD)

### 3. diabetes related illness---------------------------

# 3.1 extract needed variables
diabetes_factor = data$DIABETES

# 3.2 diabetes survival analysis

# reform data and add coveriates
data_diabetes = data[c(grep('LSI_DIAB',colnames(data)))]
data_diabetes = cbind(diabetes_factor,CNS_group,data_diabetes,sex,age,marriage,ethic,educ,lsid)
data_diabetes = subset(data_diabetes,data_diabetes$diabetes_factor==2)
times = c(1,2,3,4,5,6)
data_diabetes$time = 6
data_diabetes$status = 0
for (i in seq(10138)){
  for (j in 3:8){
    if (is.na(data_diabetes[i,j])==FALSE){
      data_diabetes$time[i] = times[j-2]
      data_diabetes$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_diabetes)
ggsurvplot(surv_model,data = data_diabetes,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_diabetes)
summary(cox_model)
ggforest(cox_model,data=data_diabetes)

### 4.dementia survival analysis

# reform data and add coveriates
data_dementia = data[c(grep('LSI_DEM',colnames(data)))]
data_dementia = cbind(CNS_group,data_dementia,sex,age,marriage,ethic,educ,lsid)
times = c(1,2,3,4,5,6)
data_dementia$time = 6
data_dementia$status = 0
for (i in seq(10308)){
  for (j in 2:7){
    if (is.na(data_dementia[i,j])==FALSE){
      data_dementia$time[i] = times[j-1]
      data_dementia$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_dementia)
ggsurvplot(surv_model,data = data_dementia,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_dementia)
summary(cox_model)
ggforest(cox_model,data=data_dementia)


### 5. parlkinson survival analysis

# reform data and add coveriates
data_parlkinson = data[c(grep('LSI_PARKSM',colnames(data)))]
data_parlkinson = cbind(CNS_group,data_parlkinson,sex,age,marriage,ethic,educ,lsid)
times = c(1,2,3,4,5,6)
data_parlkinson$time = 6
data_parlkinson$status = 0
for (i in seq(10308)){
  for (j in 2:7){
    if (is.na(data_parlkinson[i,j])==FALSE){
      data_parlkinson$time[i] = times[j-1]
      data_parlkinson$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_parlkinson)
ggsurvplot(surv_model,data = data_parlkinson,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_parlkinson)
summary(cox_model)
ggforest(cox_model,data=data_parlkinson)


### 6. stroke survival analysis

# reform data and add coveriates
data_stroke = data[c(grep('LSI_STROKE',colnames(data)))]
data_stroke = cbind(CNS_group,data_stroke,sex,age,marriage,ethic,educ,lsid)
times = c(1,2,3,4,5,6)
data_stroke$time = 6
data_stroke$status = 0
for (i in seq(10308)){
  for (j in 2:7){
    if (is.na(data_stroke[i,j])==FALSE){
      data_stroke$time[i] = times[j-1]
      data_stroke$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_stroke)
ggsurvplot(surv_model,data = data_stroke,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_stroke)
summary(cox_model)
ggforest(cox_model,data=data_stroke)

### 7. ms related illness---------------------------

# 7.1 extract needed variables
ms_factor = data$LILLC12

# 7.2 art survival analysis

# reform data and add coveriates
data_art = data[c(grep('OST_ART',colnames(data)))]
data_art = cbind(ms_factor,CNS_group,data_art,sex,age,marriage,ethic,educ,lsid)
data_art = subset(data_art,data_art$ms_factor==0)
times = c(1,2,3,4,5,6)
data_art$time = 6
data_art$status = 0
for (i in seq(6758)){
  for (j in 3:6){
    if (is.na(data_art[i,j])==FALSE){
      data_art$time[i] = times[j-2]
      data_art$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_art)
ggsurvplot(surv_model,data = data_art,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_art)
summary(cox_model)
ggforest(cox_model,data=data_art)

# 7.3 ayr survival analysis

# reform data and add coveriates
data_ayr = data[c(grep('OST_AYR',colnames(data)))]
data_ayr = cbind(ms_factor,CNS_group,data_ayr,sex,age,marriage,ethic,educ,lsid)
data_ayr = subset(data_ayr,data_ayr$ms_factor==0)
times = c(1,2,3,4,5,6)
data_ayr$time = 4
data_ayr$status = 0
for (i in seq(6758)){
  for (j in 3:6){
    if (is.na(data_ayr[i,j])==FALSE){
      data_ayr$time[i] = times[j-2]
      data_ayr$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_ayr)
ggsurvplot(surv_model,data = data_ayr,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_ayr)
summary(cox_model)
ggforest(cox_model,data=data_ayr)

# 7.4 por survival analysis

# reform data and add coveriates
data_por = data[c(grep('OST_POR',colnames(data)))]
data_por = cbind(ms_factor,CNS_group,data_por,sex,age,marriage,ethic,educ,lsid)
data_por = subset(data_por,data_por$ms_factor==0)
times = c(1,2,3,4,5,6)
data_por$time = 4
data_por$status = 0
for (i in seq(6758)){
  for (j in 3:6){
    if (is.na(data_por[i,j])==FALSE){
      data_por$time[i] = times[j-2]
      data_por$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_por)
ggsurvplot(surv_model,data = data_por,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_por)
summary(cox_model)
ggforest(cox_model,data=data_por)

# 7.5 pyr survival analysis

# reform data and add coveriates
data_pyr = data[c(grep('OST_PYR',colnames(data)))]
data_pyr = cbind(ms_factor,CNS_group,data_pyr,sex,age,marriage,ethic,educ,lsid)
data_pyr = subset(data_pyr,data_pyr$ms_factor==0)
times = c(1,2,3,4,5,6)
data_pyr$time = 4
data_pyr$status = 0
for (i in seq(6758)){
  for (j in 3:6){
    if (is.na(data_pyr[i,j])==FALSE){
      data_pyr$time[i] = times[j-2]
      data_pyr$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_pyr)
ggsurvplot(surv_model,data = data_pyr,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_pyr)
summary(cox_model)
ggforest(cox_model,data=data_pyr)


# 7.6 rart survival analysis

# reform data and add coveriates
data_rart = data[c(grep('RHE_ART',colnames(data)))]
data_rart = cbind(ms_factor,CNS_group,data_rart,sex,age,marriage,ethic,educ,lsid)
data_rart = subset(data_rart,data_rart$ms_factor==0)
times = c(1,2,3,4,5,6)
data_rart$time = 4
data_rart$status = 0
for (i in seq(6758)){
  for (j in 3:6){
    if (is.na(data_rart[i,j])==FALSE){
      data_rart$time[i] = times[j-2]
      data_rart$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_rart)
ggsurvplot(surv_model,data = data_rart,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_rart)
summary(cox_model)
ggforest(cox_model,data=data_rart)


# 7.7 rayr survival analysis

# reform data and add coveriates
data_rayr = data[c(grep('RHE_AYR',colnames(data)))]
data_rayr = cbind(ms_factor,CNS_group,data_rayr,sex,age,marriage,ethic,educ,lsid)
data_rayr = subset(data_rayr,data_rayr$ms_factor==0)
times = c(1,2,3,4,5,6)
data_rayr$time = 4
data_rayr$status = 0
for (i in seq(6758)){
  for (j in 3:6){
    if (is.na(data_rayr[i,j])==FALSE){
      data_rayr$time[i] = times[j-2]
      data_rayr$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_rayr)
ggsurvplot(surv_model,data = data_rayr,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_rayr)
summary(cox_model)
ggforest(cox_model,data=data_rayr)


# 7.8 frbever survival analysis

# reform data and add coveriates
data_frbever = data[c(grep('FRBEVER',colnames(data)))]
data_frbever = cbind(ms_factor,CNS_group,data_frbever,sex,age,marriage,ethic,educ,lsid)
data_frbever= subset(data_frbever,data_frbever$ms_factor==0)
times = c(1,2,3,4,5,6)
data_frbever$time = 3
data_frbever$status = 0
for (i in seq(6758)){
  for (j in 3:5){
    if (is.na(data_frbever[i,j])==FALSE){
      data_frbever$time[i] = times[j-2]
      data_frbever$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_frbever)
ggsurvplot(surv_model,data = data_frbever,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_frbever)
summary(cox_model)
ggforest(cox_model,data=data_frbever)


### 8. gout related illness---------------------------

# 8.1 gout survival analysis

# reform data and add coveriates
data_gout = data[c(grep('GOUT',colnames(data)))]
data_gout = data_gout[-c(grep('GOUT_YR',colnames(data_gout)))]
data_gout = cbind(CNS_group,data_gout,sex,age,marriage,ethic,educ,lsid)
times = c(1,2,3,4,5,6)
data_gout$time = 4
data_gout$status = 0
for (i in seq(10308)){
  for (j in 2:5){
    if (is.na(data_gout[i,j])==FALSE){
      data_gout$time[i] = times[j-1]
      data_gout$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_gout)
ggsurvplot(surv_model,data = data_gout,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_gout)
summary(cox_model)
ggforest(cox_model,data=data_gout)


# 8.2 gout_yr survival analysis

# reform data and add coveriates
data_gout_yr = data[c(grep('GOUT_YR',colnames(data)))]
data_gout_yr = cbind(CNS_group,data_gout_yr,sex,age,marriage,ethic,educ,lsid)
times = c(1,2,3,4,5,6)
data_gout_yr$time = 4
data_gout_yr$status = 0
for (i in seq(10308)){
  for (j in 2:5){
    if (is.na(data_gout_yr[i,j])==FALSE){
      data_gout_yr$time[i] = times[j-1]
      data_gout_yr$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_gout_yr)
ggsurvplot(surv_model,data = data_gout_yr,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_gout_yr)
summary(cox_model)
ggforest(cox_model,data=data_gout_yr)

# 8.3 gpvisyr survival analysis

# reform data and add coveriates
data_gpvisyr = data[c(grep('GPVISYR',colnames(data)))]
data_gpvisyr = cbind(CNS_group,data_gpvisyr,sex,age,marriage,ethic,educ,lsid)
times = c(1,2,3,4,5,6)
data_gpvisyr$time = 4
data_gpvisyr$status = 0
for (i in seq(10308)){
  for (j in 2:5){
    if (is.na(data_gpvisyr[i,j])==FALSE){
      data_gpvisyr$time[i] = times[j-1]
      data_gpvisyr$status[i] = 1
    }
  }
}

# plot KM curve
surv_model = survfit(Surv(time,status)~CNS_group,data=data_gpvisyr)
ggsurvplot(surv_model,data = data_gpvisyr,
           conf.int = TRUE,
           pval = TRUE,
           fun = 'pct',
           linetype = 'strata',
           size = 1,
           legend.title = 'CNSDRG groups',
           legend = c(0.3,0.2),
           legend.labs=c('use CNS less than twice','use CNS more than once'),
           surv.median.line='hv')

# cox model
cox_model = coxph(Surv(time,status)~CNS_group+age+sex+marriage+ethic+educ+lsid,data=data_gpvisyr)
summary(cox_model)
ggforest(cox_model,data=data_gpvisyr)
