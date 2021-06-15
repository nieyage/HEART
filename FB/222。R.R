
meta<-matrix(data=NA,nrow=17320,ncol=3);
colnames(meta)=c("group","cellbarcode","subtypename");
meta[,1] = as.character(hisall$subtypename)
meta[,2] = as.character(hisall$sample)
meta[,3] = "other"
for ( i in 1:nrow(meta))
{
  if (meta[i,1]=="EC" & meta[i,2]=="P1MI")
    meta[i,3] = "EC_P1MI"
  if (meta[i,1]=="EC" & meta[i,2]=="P8MI")
    meta[i,3] = "EC_P8MI"
  if (meta[i,1]=="EC" & meta[i,2]=="P1Sham")
    meta[i,3] = "EC_P1Sham"
  if (meta[i,1]=="EC" & meta[i,2]=="P8Sham")
    meta[i,3] = "EC_P8Sham"
  if (meta[i,1]=="FB" & meta[i,2]=="P1MI")
    meta[i,3] = "FB_P1MI"
  if (meta[i,1]=="FB" & meta[i,2]=="P8MI")
    meta[i,3] = "FB_P8MI"
  if (meta[i,1]=="FB" & meta[i,2]=="P1Sham")
    meta[i,3] = "FB_P1Sham"
  if (meta[i,1]=="FB" & meta[i,2]=="P8Sham")
    meta[i,3] = "FB_P8Sham"
  if (meta[i,1]=="SMC" & meta[i,2]=="P1MI")
    meta[i,3] = "SMC_P1MI"
  if (meta[i,1]=="SMC" & meta[i,2]=="P8MI")
    meta[i,3] = "SMC_P8MI"
  if (meta[i,1]=="SMC" & meta[i,2]=="P1Sham")
    meta[i,3] = "SMC_P1Sham"
  if (meta[i,1]=="SMC" & meta[i,2]=="P8Sham")
    meta[i,3] = "SMC_P8Sham"
  if (meta[i,1]=="IC" & meta[i,2]=="P1MI")
    meta[i,3] = "IC_P1MI"
  if (meta[i,1]=="IC" & meta[i,2]=="P8MI")
    meta[i,3] = "IC_P8MI"
  if (meta[i,1]=="IC" & meta[i,2]=="P1Sham")
    meta[i,3] = "IC_P1Sham"
  if (meta[i,1]=="IC" & meta[i,2]=="P8Sham")
    meta[i,3] = "IC_P8Sham"
  if (meta[i,1]=="EPI" & meta[i,2]=="P1MI")
    meta[i,3] = "EPI_P1MI"
  if (meta[i,1]=="EPI" & meta[i,2]=="P8MI")
    meta[i,3] = "EPI_P8MI"
  if (meta[i,1]=="EPI" & meta[i,2]=="P1Sham")
    meta[i,3] = "EPI_P1Sham"
  if (meta[i,1]=="EPI" & meta[i,2]=="P8Sham")
    meta[i,3] = "EPI_P8Sham"
  if (meta[i,1]=="CM" & meta[i,2]=="P1MI")
    meta[i,3] = "CM_P1MI"
  if (meta[i,1]=="CM" & meta[i,2]=="P8MI")
    meta[i,3] = "CM_P8MI"
  if (meta[i,1]=="CM" & meta[i,2]=="P1Sham")
    meta[i,3] = "CM_P1Sham"
  if (meta[i,1]=="CM" & meta[i,2]=="P8Sham")
    meta[i,3] = "CM_P8Sham"
  if (meta[i,1]=="Pericyte" & meta[i,2]=="P1MI")
    meta[i,3] = "Pericyte_P1MI"
  if (meta[i,1]=="Pericyte" & meta[i,2]=="P8MI")
    meta[i,3] = "Pericyte_P8MI"
  if (meta[i,1]=="Pericyte" & meta[i,2]=="P1Sham")
    meta[i,3] = "Pericyte_P1Sham"
  if (meta[i,1]=="Pericyte" & meta[i,2]=="P8Sham")
    meta[i,3] = "Pericyte_P8Sham"

}