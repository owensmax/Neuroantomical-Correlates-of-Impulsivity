#BUILD DATASET FOR NEUROANATOMICAL CORRELATES OF IMPULSIVE TRAITS IN CHILDREN AGED 9-10

#add upps data and filter for relevant timepoints and variables
UPPS <-read.csv("/home/max/Documents/ABCD_questionnaire_data/Y_UPPS_P.csv")
UPPS<-UPPS[ which(UPPS$redcap_event_name=='baseline_year_1_arm_1'),]
upps<-c('id_redcap','upps_y_ss_negative_urgency','upps_y_ss_lack_of_planning','upps_y_ss_sensation_seeking','upps_y_ss_positive_urgency','upps_y_ss_lack_of_perseverance')
UPPS<-UPPS[upps]
mydata<-UPPS

#add demographic data and filter for relevant timepoints and variables
demo<-read.csv("/home/max/Documents/ABCD_questionnaire_data/P_Demographics.csv")
demo<-demo[ which(demo$redcap_event_name=='baseline_year_1_arm_1'),]
mydata<-merge(demo[c(1,9)],mydata,by="id_redcap")

#add scanner data
scanner_type<-read.delim("/home/max/Documents/ABCDFixRelease2p0p1/abcd_mri01.txt",na.strings=c(""))
names(scanner_type)[4]<-"id_redcap"
scanner_type<-scanner_type[-c(1),]
scanner_type<-scanner_type[c(4,13)]
mydata<-merge(scanner_type,mydata,by="id_redcap")

#add family structure data
ACS<-read.csv("/home/max/Documents/ABCD_questionnaire_data/ACS.csv")
mydata<-merge(ACS[c(1,40)],mydata, by="id_redcap" )

mydata<-mydata[complete.cases(mydata), ]

#sMRI
#add freesurfer data  and filter for relevant timepoints and variables
smri<-read.delim("/home/max/Documents/ABCDFixRelease2p0p1/ABCDFixRelease2p0p1_renamed_files/abcd_smrip101_renamed.txt",na.strings=c(""))
#write.csv(br_mydata,"br_mydata.csv",row.names = FALSE)
smri2<-read.delim("/home/max/Documents/ABCDFixRelease2p0p1/ABCDFixRelease2p0p1_renamed_files/abcd_smrip201_renamed.txt",na.strings=c(""))
smri<-smri[c(5,11:78,225:292)]
smri2<-smri2[c(4,373,330:361)]
names(smri)[1]<-"id_redcap"
names(smri2)[1]<-"id_redcap"
smri<-merge(smri2,smri,by="id_redcap")

smri<-smri[, -grep("csf", names(smri)) ]
smri<-smri[, -grep("ventricle", names(smri)) ]
smri<-smri[, -grep("lat.vent", names(smri)) ]
smri<-smri[, -grep("white.matter", names(smri)) ]
smri<-smri[, -grep("lesion", names(smri)) ]

#merge behavioral and smri data
mydata_smri<-merge(mydata,smri,by="id_redcap")

#add quality control data
smri_qc <- read.delim("/home/max/Documents/ABCDFixRelease2p0p1/freesqc01.txt",na.strings=c(""))
cols <- smri_qc[c(4,15)]
names(cols)[1]<-"id_redcap"
cols<-cols[-c(1),]
mydata_smri <- merge(mydata_smri,cols, by="id_redcap")
mydata_smri <- mydata_smri[mydata_smri$fsqc_qc == 1, ]
mydata_smri<-mydata_smri[c(1:165)]
mydata_smri<-mydata_smri[complete.cases(mydata_smri), ]

#add dti data
dti<-read.delim("/home/max/Documents/ABCDFixRelease2p0p1/abcd_dti_p101.txt",na.strings=c(""))
dti<-dti[c(4,14:97)]
dti<-dti[-c(1),]
names(dti)[1]<-"id_redcap"

dtiqc<-read.delim("/home/max/Documents/ABCD2p0NDA/mriqcrp102.txt")
names(dtiqc)[4]<-"id_redcap"
dtiqc<-dtiqc[-c(1),]
dtiqc<-dtiqc[c(4,367)]
dtiqc[, 2] <- as.numeric(as.character( dtiqc[, 2] ))

mydata_dti<-Reduce(function(x,y) merge(x,y,by="id_redcap"), list(mydata,dti,dtiqc))
mydata_dti <- mydata_dti[mydata_dti$iqc_dmri_good_ser > 0, ]
mydata_dti<-mydata_dti[, -grep("iqc_dmri_good_ser", names(mydata_dti)) ]
mydata_dti<-mydata_dti[, -grep("_all", names(mydata_dti)) ]
mydata_dti<-mydata_dti[complete.cases(mydata_dti), ]

#winsorize data
win_mydata_smri<-mydata_smri
library('DescTools')
for (i in 10:length(win_mydata_smri))
{
  win_mydata_smri[c(i)]<-Winsorize(win_mydata_smri[,i], probs = c(0.05, 0.95),type=7)
}

win_mydata_dti<-mydata_dti
for (i in 10:length(win_mydata_dti))
{
  win_mydata_dti[c(i)]<-Winsorize(win_mydata_dti[,i], probs = c(0.05, 0.95),type=7)
}

#merge ICV into dti data and tidy up
win_mydata_dti = merge(win_mydata_smri[c(1,10)],win_mydata_dti,by='id_redcap')
win_mydata_dti = win_mydata_dti[c(1,3:10,2,11:length(win_mydata_dti))]

#change labels as needed for GAMM4 script
names(win_mydata_smri)[1]<-"src_subject_id"
names(win_mydata_smri)[4]<-'sex'
names(win_mydata_dti)[1]<-"src_subject_id"
names(win_mydata_dti)[4]<-'sex'

#make data into factors
win_mydata_smri[4]<-as.factor(win_mydata_smri$sex)
win_mydata_smri[2]<-as.factor(win_mydata_smri$rel_family_id)
win_mydata_dti[4]<-as.factor(win_mydata_dti$sex)
win_mydata_dti[2]<-as.factor(win_mydata_dti$rel_family_id)

#remove missing data
win_mydata_smri<-win_mydata_smri[complete.cases(win_mydata_smri),]
win_mydata_dti<-win_mydata_dti[complete.cases(win_mydata_dti),]

#write datasets to csvs for mixed effect and elastic net analyses
write.csv(win_mydata_smri,"/home/max/Documents/UPPS_P/fin_UPPS_smri.csv",row.names = FALSE)
write.csv(win_mydata_dti,"/home/max/Documents/UPPS_P/fin_UPPS_dti.csv",row.names = FALSE)
