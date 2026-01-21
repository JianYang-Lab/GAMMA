args <- commandArgs(trailingOnly = TRUE)
trait_name=args[1]
OUTPUT=args[2]
QTL_name_mapping_file=args[3]
gsmr_pQTL_file=args[4]
pQTL_epi_file=args[5]
gamma_gene_file=args[6]
R_functions_file=args[7]
L2G_column_file=args[8]


source(R_functions_file)
source(L2G_column_file)
suppressMessages({
library(data.table);
library(stringr);
library(dplyr)
})

### ******************************************************************************************** 
# genecode summary
## ********************************************************************************************
# 0. loading gene annotation -----------------------------------------------------
result=fread(gamma_gene_file)
result$GWAS_LOCUS=NA;result$Lead_SNP=NA;result$Lead_SNP_BP=NA;

# 1. COJO locus -----------------------------------------------------
COJO_locus=fread(paste0(OUTPUT, "/Clumping/summary/",trait_name,".locus"), header=T)
for(j in 1:nrow(COJO_locus)){
  chr=as.numeric(gsub("chr","",COJO_locus$chr[j]))
  start=COJO_locus$start[j]
  end=COJO_locus$end[j]
  locus=COJO_locus$GWAS_LOCUS[j]
  lead_snp=COJO_locus$Lead_SNP[j]
  lead_snp_bp=COJO_locus$Lead_SNP_BP[j]

  index=which(result$chr==chr & result$start<=end & result$start>=start & result$end<=end & result$end>=start)
  result$GWAS_LOCUS[index]=locus
  result$Lead_SNP[index]=lead_snp
  result$Lead_SNP_BP[index]=lead_snp_bp
}


## ********************************************************************************************
# GAMMA L2G summary
## ********************************************************************************************
# L2G QTL name mapping list
L2G_QTL_name_mapping_list=fread(QTL_name_mapping_file,head=F,stringsAsFactors=F,data.table=F)

# *****************************
# 1. loading SMR results
# *****************************
# Remember to add SMR results data treatment script
# SMR_dir=paste0(OUTPUT,"/SMR/SMR_Portal/")
# SMR_summary=fread(paste0(SMR_dir, trait_name, "_SMR_plot.summary"))
# index=match(result$gene_name, SMR_summary$gene_name, nomatch=0)

SMR_DIRT=paste0(OUTPUT,"/SMR/summary/")
# file_list1=list.files(SMR_DIRT)
# file_list=file_list1[grep(paste0("^",trait_name,"_(.*)_chrALL\\.msmr$"),file_list1)]

qtl_name_list=SMR_QTL_name_list

result_b_SMR=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_se_SMR=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_z_SMR=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_p_SMR=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_p_HEIDI=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_probeID_SMR=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_SMR=data.frame()
# qtl_name_list=c()

for(i in 1:length(qtl_name_list)){  
  # print(i)
  qtl_name = qtl_name_list[i]
  infile = paste0(trait_name,"_",qtl_name,"_chrALL.msmr")
  INFILE=paste0(SMR_DIRT,infile)
  
  # qtl_name=gsub(paste0("^",trait_name,"_(.*)_chrALL\\.msmr$"), "\\1", infile)
  # qtl_name_list=c(qtl_name_list, qtl_name)

  ############################
  # check if file is empty
  if(!file.exists(INFILE)) {
    data <- data.frame()
  } else {
    file_content <- readLines(INFILE)
  
    if (length(file_content) == 0 || all(trimws(file_content) == "")) {
      print(paste0(INFILE, ": The file is empty or contains only whitespace."))
      data <- data.frame()
    } else {
      data <- fread(INFILE, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
    }
  }

  
  if(nrow(data)>0){
  data$QTL_name=qtl_name
  data=data %>% group_by(Gene) %>% filter(p_SMR==min(p_SMR)) %>% ungroup() %>% data.frame()

  index=match(result$gene_name,data$Gene,nomatch=0)
  result_b_SMR[which(index!=0),i]=data$b_SMR[index]
  result_se_SMR[which(index!=0),i]=data$se_SMR[index]
  result_z_SMR[which(index!=0),i]=data$b_SMR[index]/data$se_SMR[index]
  result_p_SMR[which(index!=0),i]=data$p_SMR[index]
  result_p_HEIDI[which(index!=0),i]=data$p_HEIDI[index]
  result_probeID_SMR[which(index!=0),i]=data$probeID[index]

  # For illustration
  data_tmp=data[index,c("QTL_name", "probeID","ProbeChr","Gene","b_SMR","se_SMR","p_SMR","p_SMR_multi","p_HEIDI")]
  index_sig=which(data_tmp$p_SMR <= 0.05/nrow(data_tmp))
  result_SMR=rbind(result_SMR, data_tmp[index_sig,])
  }
}

index=match(qtl_name_list, L2G_QTL_name_mapping_list$V1, nomatch=0)
qtl_name_list[which(index!=0)]=L2G_QTL_name_mapping_list$V2[index]

index=match(result_SMR$QTL_name, L2G_QTL_name_mapping_list$V1, nomatch=0)
result_SMR$QTL_name[which(index!=0)]=L2G_QTL_name_mapping_list$V2[index]

colnames(result_b_SMR)=colnames(result_se_SMR)=colnames(result_z_SMR)=colnames(result_p_SMR)=colnames(result_p_HEIDI)=colnames(result_probeID_SMR)=qtl_name_list
rownames(result_b_SMR)=rownames(result_se_SMR)=rownames(result_z_SMR)=rownames(result_p_SMR)=rownames(result_p_HEIDI)=rownames(result_probeID_SMR)=result$gene_name


# -------------------
# all 19054 genes
write.csv(result_b_SMR, paste0(OUTPUT,"/L2G/summary/",trait_name, "_b_SMR.csv"), row.names = TRUE)
write.csv(result_se_SMR, paste0(OUTPUT,"/L2G/summary/",trait_name, "_se_SMR.csv"), row.names = TRUE)
write.csv(result_z_SMR, paste0(OUTPUT,"/L2G/summary/",trait_name, "_z_SMR.csv"), row.names = TRUE)
write.csv(result_p_SMR, paste0(OUTPUT,"/L2G/summary/",trait_name, "_p_SMR.csv"), row.names = TRUE)
write.csv(result_p_HEIDI, paste0(OUTPUT,"/L2G/summary/",trait_name, "_p_HEIDI.csv"), row.names = TRUE)
write.csv(result_probeID_SMR, paste0(OUTPUT,"/L2G/summary/",trait_name, "_probeID_SMR.csv"), row.names = TRUE)

write.table(result_SMR, paste0(OUTPUT,"/L2G/summary/",trait_name, "_SMR.txt"), row=F,col=T,quo=F,sep="\t")

result$eSMR=result$sSMR=result$pSMR=NA;
result$eSMR_p_HEIDI=result$sSMR_p_HEIDI=result$pSMR_p_HEIDI=NA
result$eSMR_QTL_name=result$sSMR_QTL_name=result$pSMR_QTL_name=NA
result$eSMR_probeID=result$sSMR_probeID=result$pSMR_probeID=NA

# index1=grep("eQTL",qtl_name_list)
# index0=grep("DIRECT",qtl_name_list)
# eQTL_index=unique(c(index0,index1))
# eQTL_index=index1
eQTL_index=grep("eQTL",qtl_name_list)
sQTL_index=grep("sQTL",qtl_name_list)
pQTL_index=grep("pQTL",qtl_name_list)

for(i in 1:nrow(result)){
	# print(i)
  if(all(is.na(result_p_SMR[i,]))){
    next
  }else{
    if(!all(is.na(result_p_SMR[i,eQTL_index]))){
    result$eSMR[i]=min(result_p_SMR[i,eQTL_index, drop=FALSE],na.rm=T)
    result$eSMR_QTL_name[i]=colnames(result_p_SMR)[which.min(result_p_SMR[i,eQTL_index])]
    result$eSMR_p_HEIDI[i]=result_p_HEIDI[i, which.min(result_p_SMR[i,eQTL_index])]
    result$eSMR_probeID[i]=result_probeID_SMR[i, which.min(result_p_SMR[i,eQTL_index])]
    }
    if(!all(is.na(result_p_SMR[i,sQTL_index]))){
    result$sSMR[i]=min(result_p_SMR[i,sQTL_index, drop=FALSE],na.rm=T);
    result$sSMR_QTL_name[i]=colnames(result_p_SMR)[which.min(result_p_SMR[i,sQTL_index])]
    result$sSMR_p_HEIDI[i]=result_p_HEIDI[i, which.min(result_p_SMR[i,sQTL_index])]
    result$sSMR_probeID[i]=result_probeID_SMR[i, which.min(result_p_SMR[i,sQTL_index])]    
    }
    if(!all(is.na(result_p_SMR[i,pQTL_index]))){
    result$pSMR[i]=min(result_p_SMR[i,pQTL_index, drop=FALSE],na.rm=T)
    result$pSMR_QTL_name[i]=colnames(result_p_SMR)[which.min(result_p_SMR[i,pQTL_index])]
    result$pSMR_p_HEIDI[i]=result_p_HEIDI[i, which.min(result_p_SMR[i,pQTL_index])]
    result$pSMR_probeID[i]=result_probeID_SMR[i, which.min(result_p_SMR[i,pQTL_index])]    
    }
  }
}

SMR_feature=cbind(result_z_SMR,result_p_SMR,result_p_HEIDI)
colnames(SMR_feature)=c(paste0("z_SMR_", colnames(result_z_SMR)),paste0("p_SMR_", colnames(result_p_SMR)),paste0("p_HEIDI_SMR_", colnames(result_p_HEIDI)))

# *****************************
# 2. loading COLOC results
# *****************************

COLOC_DIRT=paste0(OUTPUT,"/COLOC/summary/")
# file_list1=list.files(COLOC_DIRT)
# file_list=file_list1[grep(paste0("^",trait_name,"_(.*)_chrALL\\.coloc$"),file_list1)]
# file_list=file_list[-grep("0.05", file_list)]

qtl_name_list=COLOC_QTL_name_list

result_pp4_COLOC=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_probeID_COLOC=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_COLOC=data.frame()
# qtl_name_list=c()

for(i in 1:length(qtl_name_list)){
  # print(i)
  qtl_name = qtl_name_list[i]
  infile = paste0(trait_name,"_",qtl_name,"_chrALL.coloc")
  INFILE=paste0(COLOC_DIRT,infile) 
  # print(INFILE)

  # qtl_name=gsub(paste0("^",trait_name,"_(.*)\\.coloc$"), "\\1", infile)
  # qtl_name_list=c(qtl_name_list, qtl_name)

  ############################
  # check if file is empty
  if(!file.exists(INFILE)) {
    data <- data.frame()
  } else {
    file_content <- readLines(INFILE)
  
    if (length(file_content) == 0 || all(trimws(file_content) == "")) {
      print(paste0(INFILE, ": The file is empty or contains only whitespace."))
      data <- data.frame()
    } else {
      data <- fread(INFILE, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
    }
  }

    
  if(nrow(data)>0){
  colnames(data)=c("probe", "gene_name", "pp0", "pp1", "pp2", "pp3", "pp4", "PP3_PP4")
  data$QTL_name=qtl_name
  data=data %>% group_by(gene_name) %>% filter(pp4==max(pp4)) %>% ungroup() %>% data.frame()
  index=match(result$gene_name,data$gene_name,nomatch=0)
  result_pp4_COLOC[which(index!=0),i]=data$pp4[index]
  result_probeID_COLOC[which(index!=0),i]=data$probe[index]
  
  data_tmp=data[index,]
  index_sig=which(data_tmp$pp4 >= 0.8)
  result_COLOC=rbind(result_COLOC, data_tmp[index_sig,])
  }
}

index=match(result_COLOC$gene_name, result$gene_name, nomatch=0)
result_COLOC$chr[which(index!=0)]=result$chr[index]
if(nrow(result_COLOC) > 0) {
    result_COLOC = result_COLOC[, c("QTL_name", "probe", "chr", "gene_name", "pp0", "pp1", "pp2", "pp3", "pp4")]
} else {
    print("The result_COLOC dataframe is empty.")
}
colnames(result_COLOC)[2]="probeID"

index=match(qtl_name_list, L2G_QTL_name_mapping_list$V1, nomatch=0)
qtl_name_list[which(index!=0)]=L2G_QTL_name_mapping_list$V2[index]

index=match(result_COLOC$QTL_name, L2G_QTL_name_mapping_list$V1, nomatch=0)
result_COLOC$QTL_name[which(index!=0)]=L2G_QTL_name_mapping_list$V2[index]

colnames(result_pp4_COLOC)=colnames(result_probeID_COLOC)=qtl_name_list
rownames(result_pp4_COLOC)=rownames(result_probeID_COLOC)=result$gene_name

write.csv(result_pp4_COLOC, paste0(OUTPUT,"/L2G/summary/",trait_name, "_pp4_COLOC.csv"), row.names = TRUE)
write.csv(result_probeID_COLOC, paste0(OUTPUT,"/L2G/summary/",trait_name, "_probeID_COLOC.csv"), row.names = TRUE)
write.table(result_COLOC, paste0(OUTPUT,"/L2G/summary/",trait_name, "_COLOC.txt"), row=F,col=T,quo=F,sep="\t")


# COLOC_index<-apply(result_COLOC,1,function(x) all(is.na(x)))
# GAMMA_xQTL result ------------------------------------------------------------
result$eCOLOC=result$sCOLOC=result$pCOLOC=NA
result$eCOLOC_QTL_name=result$sCOLOC_QTL_name=result$pCOLOC_QTL_name=NA
result$eCOLOC_probeID=result$sCOLOC_probeID=result$pCOLOC_probeID=NA

eQTL_index=grep("eQTL",qtl_name_list)
sQTL_index=grep("sQTL",qtl_name_list)
pQTL_index=grep("pQTL",qtl_name_list)

for(i in 1:nrow(result)){
	# print(i)
  if(all(is.na(result_pp4_COLOC[i,]))){
    next
  }else{ 
    if(!all(is.na(result_pp4_COLOC[i,eQTL_index]))){
    result$eCOLOC[i]=max(result_pp4_COLOC[i,eQTL_index, drop=FALSE],na.rm=T);
    result$eCOLOC_QTL_name[i]=colnames(result_pp4_COLOC)[which.max(result_pp4_COLOC[i,eQTL_index])]
    result$eCOLOC_probeID[i]=result_probeID_COLOC[i,which.max(result_pp4_COLOC[i,eQTL_index])]
    }
    if(!all(is.na(result_pp4_COLOC[i,sQTL_index]))){
    result$sCOLOC[i]=max(result_pp4_COLOC[i,sQTL_index, drop=FALSE],na.rm=T)
    result$sCOLOC_QTL_name[i]=colnames(result_pp4_COLOC)[which.max(result_pp4_COLOC[i,sQTL_index])]
    result$sCOLOC_probeID[i]=result_probeID_COLOC[i,which.max(result_pp4_COLOC[i,sQTL_index])]
    }
    if(!all(is.na(result_pp4_COLOC[i,pQTL_index]))){
    result$pCOLOC[i]=max(result_pp4_COLOC[i,pQTL_index, drop=FALSE],na.rm=T)
    result$pCOLOC_QTL_name[i]=colnames(result_pp4_COLOC)[which.max(result_pp4_COLOC[i,pQTL_index])]
    result$pCOLOC_probeID[i]=result_probeID_COLOC[i,which.max(result_pp4_COLOC[i,pQTL_index])]
    }
  }
}

COLOC_feature=result_pp4_COLOC
colnames(COLOC_feature)=paste0("pp4_COLOC_",colnames(COLOC_feature))

# *****************************
# 3. loading FUSION results
# *****************************

FUSION_DIRT=paste0(OUTPUT,"/FUSION/summary/")
# file_list1=list.files(FUSION_DIRT)
# file_list=file_list1[grep(paste0("^",trait_name,"\\."),file_list1)]

qtl_name_list=FUSION_QTL_name_list

result_z_FUSION=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_p_FUSION=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_probeID_FUSION=matrix(NA,nrow=nrow(result),ncol=length(qtl_name_list))
result_FUSION=data.frame()
# qtl_name_list=c()

for(i in 1:length(qtl_name_list)){
	# print(i)
  qtl_name = qtl_name_list[i]
  infile = paste0(trait_name,".GTExv8.ALL.",qtl_name,".chrALL.dat")
  INFILE=paste0(FUSION_DIRT,infile)
  
  # qtl_name=gsub(paste0("^",trait_name,"\\.GTExv8\\.ALL\\.","(.*)\\.chrALL\\.dat$"), "\\1", infile)
  # qtl_name_list=c(qtl_name_list,qtl_name)


  ############################
  # check if file is empty
  if(!file.exists(INFILE)) {
    data <- data.frame()
  } else {
    file_content <- readLines(INFILE)
  
    if (length(file_content) == 0 || all(trimws(file_content) == "")) {
      print(paste0(INFILE, ": The file is empty or contains only whitespace."))
      data <- data.frame()
    } else {
      data <- fread(INFILE, header = TRUE, stringsAsFactors = FALSE, data.table = FALSE)
    }
  }


  if(nrow(data)>0){
  data$QTL_name=qtl_name
  data=data %>% group_by(ID) %>% filter(TWAS.P==min(TWAS.P)) %>% ungroup() %>% data.frame()

  if (nrow(data) > 0){
  data$gene_ID=str_split_fixed(data$ID,"\\.",Inf)[,1]

  index=match(result$gene_id,data$gene_ID,nomatch=0)
  result_z_FUSION[which(index!=0),i]=data$TWAS.Z[index]
  result_p_FUSION[which(index!=0),i]=data$TWAS.P[index]
  result_probeID_FUSION[which(index!=0),i]=data$ID[index]

  # For illustration
  data_tmp=data[index,]
  index_sig=which(data_tmp$TWAS.P <= 0.05/nrow(data_tmp))
  result_FUSION=rbind(result_FUSION, data_tmp[index_sig,])
  }
  }
}

index=match(result_FUSION$gene_ID, result$gene_id, nomatch=0)
result_FUSION$gene_name[which(index!=0)]=result$gene_name[index]
if(nrow(result_FUSION) > 0) {
    result_FUSION=result_FUSION[,c("QTL_name", "ID","CHR","gene_name","TWAS.Z","TWAS.P")]
    colnames(result_FUSION)[2]="probeID"
} else {
    print("The result_FUSION dataframe is empty.")
}


index=match(qtl_name_list, L2G_QTL_name_mapping_list$V1, nomatch=0)
qtl_name_list[which(index!=0)]=L2G_QTL_name_mapping_list$V2[index]

index=match(result_FUSION$QTL_name, L2G_QTL_name_mapping_list$V1, nomatch=0)
result_FUSION$QTL_name[which(index!=0)]=L2G_QTL_name_mapping_list$V2[index]

colnames(result_z_FUSION)=colnames(result_p_FUSION)=colnames(result_probeID_FUSION)=qtl_name_list
rownames(result_z_FUSION)=rownames(result_p_FUSION)=rownames(result_probeID_FUSION)=result$gene_name

write.csv(result_z_FUSION, paste0(OUTPUT,"/L2G/summary/",trait_name, "_z_FUSION.csv"), row.names = TRUE)
write.csv(result_p_FUSION, paste0(OUTPUT,"/L2G/summary/",trait_name, "_p_FUSION.csv"), row.names = TRUE)
write.csv(result_probeID_FUSION, paste0(OUTPUT,"/L2G/summary/",trait_name, "_probeID_FUSION.csv"), row.names = TRUE)

write.table(result_FUSION, paste0(OUTPUT,"/L2G/summary/",trait_name, "_FUSION.txt"), row=F,col=T,quo=F,sep="\t")
# FUSION_index<-apply(result_FUSION,1,function(x) all(is.na(x)))

# GAMMA_xQTL result ------------------------------------------------------------
result$FUSION=NA;
result$FUSION_QTL_name=NA;
result$FUSION_probeID=NA;

for(i in 1:nrow(result)){
	# print(i)
  if(all(is.na(result_p_FUSION[i,]))){
    next
  }else{
    result$FUSION[i]=min(result_p_FUSION[i,], na.rm = T)
    result$FUSION_QTL_name[i]=colnames(result_p_FUSION)[which.min(result_p_FUSION[i,])]
    result$FUSION_probeID[i]=result_probeID_FUSION[i,which.min(result_p_FUSION[i,])]
  }
}

FUSION_feature=cbind(result_z_FUSION, result_p_FUSION)
colnames(FUSION_feature)=c(paste0("z_FUSION_",colnames(result_z_FUSION)),paste0("p_FUSION_",colnames(result_p_FUSION)))

# *****************************
# 4. loading GSMR results
# *****************************

gsmr_pQTL_file="/storage/yangjianLab/guoyazhou/GAMMA_git_data/GSMR/data/gsmr_pQTL.txt"
gsmr_pQTL=fread(gsmr_pQTL_file, head=F, stringsAsFactors = F, data.table = F)
gsmr_pQTL$gene_name=NA
pQTL_epi_file="/storage/yangjianLab/sharedata/molecular_QTL/pQTL/freeze1/merged_pQTL_plasma/merged_pQTL_plasma_5pQTL.epi"
epi=fread(pQTL_epi_file,head=F,stringsAsFactors=F,data.table=F)

index=match(gsmr_pQTL$V1, epi$V2, nomatch=0)
gsmr_pQTL$gene_name[which(index!=0)]=epi$V5[index]
 
# unique(str_split_fixed(gsmr_pQTL[which(is.na(gsmr_pQTL$gene_name)),"V2"],"/",Inf)[,7])
index=which(str_split_fixed(gsmr_pQTL$V2,"/",Inf)[,7]=="pQTL_DIRECT")
gsmr_pQTL$gene_name[index]=str_split_fixed(gsmr_pQTL$V1[index],"_",Inf)[,1]
index=which(str_split_fixed(gsmr_pQTL$V2,"/",Inf)[,7]=="pQTL_NIAGADS")
gsmr_pQTL$gene_name[index]=str_split_fixed(gsmr_pQTL$V1[index],"_",Inf)[,2]
index=which(str_split_fixed(gsmr_pQTL$V2,"/",Inf)[,7]=="UKB_PPP")
gsmr_pQTL$gene_name[index]=str_split_fixed(gsmr_pQTL$V1[index],"_",Inf)[,1]


# ------------------------------------------------
# GSMR results -----------------------------------
result_GSMR=data.frame()
result_GSMR_feature=matrix(NA,nrow=nrow(result),ncol=2)

GSMR_file=paste0(OUTPUT,"/GSMR/summary/",trait_name,".gsmr")

if (file.exists(GSMR_file)){

  data=fread(GSMR_file,head=T,stringsAsFactors=F,data.table=F)
  data=data[which(data$bxy!="NaN"),]
  # dim(data)
  # head(data)

  index=match(data$Exposure,gsmr_pQTL$V1,nomatch=0)
  data$gene_name=data$Exposure
  data$gene_name[which(index!=0)]=gsmr_pQTL$gene_name[index]

  data=data[which(!is.na(data$bxy)),]
  maxcpp=data.frame(data %>% group_by(gene_name) %>% filter(p==min(p)))

  index=match(result$gene_name,maxcpp$gene_name,nomatch=0)
  result_GSMR_feature[which(index!=0),1]=maxcpp$bxy[index]/maxcpp$se[index]
  result_GSMR_feature[which(index!=0),2]=maxcpp$p[index]

  result_GSMR=rbind(result_GSMR,maxcpp[index,])
  result_GSMR=result_GSMR[,c("Exposure","gene_name","Outcome","bxy","se","p","nsnp","multi_snp_based_heidi_outlier")]
}


colnames(result_GSMR_feature)=c("z_GSMR","p_GSMR")
rownames(result_GSMR_feature)=result$gene_name

write.csv(result_GSMR_feature, paste0(OUTPUT,"/L2G/summary/",trait_name, "_GSMR_feature.csv"), row.names = TRUE)
write.table(result_GSMR, paste0(OUTPUT,"/L2G/summary/",trait_name, "_GSMR.txt"), row=F,col=T,quo=F,sep="\t")


# GSMR_index<-apply(result_GSMR,1,function(x) all(is.na(x)))
# GAMMA_xQTL result ------------------------------------------------------------
result$GSMR=NA;

for(i in 1:nrow(result)){
  if(all(is.na(result_GSMR_feature[i,]))){
    next
  }else{
    result$GSMR[i]=result_GSMR_feature[i,2]
  }
}

# GAMMA_xQTL feature ------------------------------------------------------------
GSMR_feature=result_GSMR_feature


# # ------------------------------------------------
# # GSMR results -----------------------------------
# result_GSMR=data.frame()
# result_GSMR_feature=matrix(NA,nrow=nrow(result),ncol=2)

# GSMR_file=paste0(OUTPUT,"/MR/GSMR/results/",trait_name,".gsmr")

# if (file.exists(GSMR_file)){
#   data=fread(GSMR_file,head=T,stringsAsFactors=F,data.table=F)
#   data=data[which(data$bxy!="NaN"),]
#   data=data[which(!is.na(data$bxy)),]

#   maxcpp=data.frame(data %>% group_by(gene_name) %>% filter(pval==min(pval)))

#   index=match(result$gene_name,maxcpp$gene_name,nomatch=0)
#   result_GSMR_feature[which(index!=0),1]=maxcpp$bxy[index]/maxcpp$se[index]
#   result_GSMR_feature[which(index!=0),2]=maxcpp$pval[index]

#   result_GSMR=rbind(result_GSMR,maxcpp[index,])
#   result_GSMR=result_GSMR[,c("probeID","gene_name","bxy","se","pval","nsnps")]
# }


# colnames(result_GSMR_feature)=c("z_GSMR","p_GSMR")
# rownames(result_GSMR_feature)=result$gene_name

# write.csv(result_GSMR_feature, paste0(OUTPUT,"/L2G/summary/",trait_name, "_MR_GSMR_feature.csv"), row.names = TRUE)
# write.table(result_GSMR, paste0(OUTPUT,"/L2G/summary/",trait_name, "_MR_GSMR.txt"), row=F,col=T,quo=F,sep="\t")


# # GSMR_index<-apply(result_GSMR,1,function(x) all(is.na(x)))
# # GAMMA_xQTL result ------------------------------------------------------------
# result$GSMR=NA;

# for(i in 1:nrow(result)){
#   if(all(is.na(result_GSMR_feature[i,]))){
#     next
#   }else{
#     result$GSMR[i]=result_GSMR_feature[i,2]
#   }
# }

# # GAMMA_xQTL feature ------------------------------------------------------------
# GSMR_feature=result_GSMR_feature



# *****************************
# 5. loading MAGIC results
# *****************************
# -----------------------------------------------------
result$MAGIC=result$eMAGIC=result$sMAGIC=result$pMAGIC=result$mMAGIC=result$caMAGIC=result$hMAGIC=NA
result$eMAGIC_QTL_name=result$sMAGIC_QTL_name=result$pMAGIC_QTL_name=result$mMAGIC_QTL_name=result$caMAGIC_QTL_name=result$hMAGIC_QTL_name=NA
result$eMAGIC_probeID=result$sMAGIC_probeID=result$pMAGIC_probeID=result$mMAGIC_probeID=result$caMAGIC_probeID=result$hMAGIC_probeID=NA

MAGIC_file=paste0(OUTPUT,"/MAGIC/summary/",trait_name,"_MAGIC.txt")

if(file.exists(MAGIC_file)){
data=fread(MAGIC_file,head=T,stringsAsFactors=F,data.table=F)

index=match(result$gene_name, data$gene_name, nomatch=0)
result$MAGIC[which(index!=0)]=data$MAGIC[index]
result$eMAGIC[which(index!=0)]=data$eMAGIC[index]
result$sMAGIC[which(index!=0)]=data$sMAGIC[index]
result$pMAGIC[which(index!=0)]=data$pMAGIC[index]
result$mMAGIC[which(index!=0)]=data$mMAGIC[index]
result$caMAGIC[which(index!=0)]=data$caMAGIC[index]
result$hMAGIC[which(index!=0)]=data$hMAGIC[index]

result$eMAGIC_QTL_name[which(index!=0)]=data$eMAGIC_QTL_name[index]
result$sMAGIC_QTL_name[which(index!=0)]=data$sMAGIC_QTL_name[index]
result$pMAGIC_QTL_name[which(index!=0)]=data$pMAGIC_QTL_name[index]
result$mMAGIC_QTL_name[which(index!=0)]=data$mMAGIC_QTL_name[index]
result$caMAGIC_QTL_name[which(index!=0)]=data$caMAGIC_QTL_name[index]
result$hMAGIC_QTL_name[which(index!=0)]=data$hMAGIC_QTL_name[index]

result$eMAGIC_probeID[which(index!=0)]=data$eMAGIC_probeID[index]
result$sMAGIC_probeID[which(index!=0)]=data$sMAGIC_probeID[index]
result$pMAGIC_probeID[which(index!=0)]=data$pMAGIC_probeID[index]
result$mMAGIC_probeID[which(index!=0)]=data$mMAGIC_probeID[index]
result$caMAGIC_probeID[which(index!=0)]=data$caMAGIC_probeID[index]
result$hMAGIC_probeID[which(index!=0)]=data$hMAGIC_probeID[index]
}

MAGIC_feature=result[,c("MAGIC","eMAGIC", "sMAGIC", "pMAGIC", "mMAGIC", "caMAGIC", "hMAGIC")]

## ********************************************************************************************
# Construct GAMMA results xQTL score
## ********************************************************************************************

result1=data.frame(matrix(0,nrow=nrow(result),ncol=13))
result1[,1]=ifelse(result$eSMR < 0.05/length(which(!is.na(result$eSMR))), 1, 0)
result1[,2]=ifelse(result$sSMR < 0.05/length(which(!is.na(result$sSMR))), 1, 0)
result1[,3]=ifelse(result$pSMR < 0.05/length(which(!is.na(result$pSMR))), 1, 0)
result1[,4]=ifelse(result$eCOLOC >= 0.8, 1, 0)
result1[,5]=ifelse(result$sCOLOC >= 0.8, 1, 0)
result1[,6]=ifelse(result$pCOLOC >= 0.8, 1, 0)
result1[,7]=ifelse(result$FUSION < 0.05/length(which(!is.na(result$FUSION))), 1, 0)
result1[,8]=ifelse(result$GSMR < 0.05/length(which(!is.na(result$GSMR))), 1, 0)
result1[,9]=ifelse(result$MAGIC < 0.05/length(which(!is.na(result$MAGIC))), 1, 0)
result1[,10]=ifelse(result$eMAGIC < 0.05/length(which(!is.na(result$eMAGIC))), 1, 0)
result1[,11]=ifelse(result$sMAGIC < 0.05/length(which(!is.na(result$sMAGIC))), 1, 0)
result1[,12]=ifelse(result$pMAGIC < 0.05/length(which(!is.na(result$pMAGIC))), 1, 0)
result1[,13]=ifelse(result$mMAGIC < 0.05/length(which(!is.na(result$mMAGIC))), 1, 0)
result1[,14]=ifelse(result$caMAGIC < 0.05/length(which(!is.na(result$caMAGIC))), 1, 0)
result1[,15]=ifelse(result$hMAGIC < 0.05/length(which(!is.na(result$hMAGIC))), 1, 0)

result$GAMMA_xQTL=apply(result1,1,function(x) sum(x,na.rm=T))
feature=result[,c("Gene_ID","gene_id","gene_name","chr", "start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP",
                  "GAMMA_xQTL",
                  "eSMR","sSMR","pSMR","eCOLOC","sCOLOC","pCOLOC","FUSION","GSMR",
                  "MAGIC","eMAGIC", "sMAGIC", "pMAGIC", "mMAGIC", "caMAGIC", "hMAGIC")]

feature=cbind(feature, SMR_feature, COLOC_feature, FUSION_feature,GSMR_feature)
# feature=cbind(feature, SMR_feature, COLOC_feature, FUSION_feature,GSMR_feature, MAGIC_feature)
# feature$GAMMA_xQTL=apply(result1,1,function(x) sum(x,na.rm=T))
# SMR_featuer
# COLOC_feature
# FUSION_feature
# GSMR_feature
# MAGIC_feature

GAMMA_xQTL_name_list=c(## Gene information 
	"Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP",
	## GAMMA xQTL information
	"GAMMA_xQTL", "eSMR","sSMR","pSMR","eCOLOC","sCOLOC","pCOLOC","FUSION","GSMR", 
  "MAGIC","eMAGIC", "sMAGIC", "pMAGIC", "mMAGIC", "caMAGIC", "hMAGIC", 
	## GAMMA xQTL corresponding QTL name
	"eSMR_QTL_name","sSMR_QTL_name","pSMR_QTL_name", "eCOLOC_QTL_name", "sCOLOC_QTL_name", "pCOLOC_QTL_name","FUSION_QTL_name",
  "eMAGIC_QTL_name", "sMAGIC_QTL_name", "pMAGIC_QTL_name", "mMAGIC_QTL_name", "caMAGIC_QTL_name", "hMAGIC_QTL_name", 
  "eSMR_p_HEIDI","sSMR_p_HEIDI","pSMR_p_HEIDI","eSMR_probeID","sSMR_probeID","pSMR_probeID","eCOLOC_probeID","sCOLOC_probeID","pCOLOC_probeID","FUSION_probeID",
  "eMAGIC_probeID", "sMAGIC_probeID", "pMAGIC_probeID", "mMAGIC_probeID",  "caMAGIC_probeID", "hMAGIC_probeID")

write.table(result[,..GAMMA_xQTL_name_list],paste0(OUTPUT,"/L2G/score/",trait_name,"_GAMMA_xQTL.summary"),row=F,col=T,quo=F,sep="\t")
write.table(feature, paste0(OUTPUT,"/L2G/feature/",trait_name,"_GAMMA_xQTL.feature"),row=F,col=T,quo=F,sep="\t")



