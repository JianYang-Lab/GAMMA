args=commandArgs(TRUE)
trait_name=args[1]
OUTPUT=args[2]
magic_functions_file=args[3]
gencode_file=args[4]
CpG_link_file=args[5]
hQTL_link_file=args[6]
caQTL_link_file=args[7]
# ---- gwas locus data ----
GWAS_DATA=args[8]
reference_bim_file=args[9]
QTL_name_list_file=args[10]

suppressMessages({
library(data.table);
library(stringr);
library(dplyr);
library(qvalue);
library(tidyverse);
library(GenomicRanges);
library(regioneR);
library(Repitools);
library(mgcv)
})

source(magic_functions_file)

####################################################################################################
# GWAS and bim file
####################################################################################################
gwas=fread(GWAS_DATA,head=TRUE,stringsAsFactors=FALSE,data.table=FALSE)
colnames(gwas)=toupper(colnames(gwas))
gwas$P = as.numeric(gwas$P)

bim=fread(reference_bim_file,head=F,stringsAsFactors=F,data.table=F)

index=match(gwas$SNP, bim$V2, nomatch=0)
gwas$CHR=gwas$POS=NA
gwas$CHR[which(index!=0)]=bim$V1[index]
gwas$POS[which(index!=0)]=bim$V4[index]

####################################################################################################
# Gencode annotations
####################################################################################################
anot=fread(gencode_file, head=T,stringsAsFactors=F,data.table=F)
anot$gene_id=str_split_fixed(anot$gene_id,"\\.",Inf)[,1]
anot <- anot[anot$gene_type == "protein_coding", ]
anot <- anot[!anot[, 1] %in% c("chrX", "chrY", "chrM"), ]
anot <- anot[!duplicated(anot$gene_name), ]
anot$V1 <- as.numeric(sub("^chr", "", anot$V1))

result <- data.frame(chr=anot$V1,start=anot$V4,end=anot$V5,strand=anot$V7,gene_id=anot$gene_id,gene_name=anot$gene_name)
result$GWAS_LOCUS=NA;result$Lead_SNP=NA;result$Lead_SNP_BP=NA;
Locus_data=fread(paste0(OUTPUT, "/Clumping/summary/",trait_name,".locus"), header=T)
for(j in 1:nrow(Locus_data)){
  # print(j)
  chr=sub("^chr", "", Locus_data$chr[j])
  start=Locus_data$start[j]
  end=Locus_data$end[j]
  locus=Locus_data$GWAS_LOCUS[j]
  lead_snp=Locus_data$Lead_SNP[j]
  lead_snp_bp=Locus_data$Lead_SNP_BP[j]

  index=which(result$chr==chr & result$start<=end & result$start>=start & result$end<=end & result$end>=start)
  result$GWAS_LOCUS[index]=locus
  result$Lead_SNP[index]=lead_snp
  result$Lead_SNP_BP[index]=lead_snp_bp

  gwas_index=which(gwas$CHR==chr & gwas$POS<=end & gwas$POS>=start)
  gwas_locus_data=gwas[gwas_index,c("CHR","POS","P","SNP")]
  write.table(gwas_locus_data,paste0(OUTPUT,"/MAGIC/gwas/",trait_name,"_", locus,".txt"),row=F,col=T,quo=F,sep="\t")
}

####################################################################################################
# Link files
####################################################################################################

# mQTL/hQTL/caQTL-trait SMR results (probeID is not the same;)
# link mQTL/hQTL/caQTL probeIDs to genes
CpG_link=fread(CpG_link_file,head=F,stringsAsFactors=F,data.table=F)
hQTL_link=fread(hQTL_link_file,head=F,stringsAsFactors=F,data.table=F)
caQTL_link=fread(caQTL_link_file,head=F,stringsAsFactors=F,data.table=F)

####################################################################################################
# read SMR assocaition between molecular trait and complex trait
####################################################################################################
SMR_DIRT=paste0(OUTPUT,"/SMR/summary/")

# eQTL/sQTL/pQTL-trait SMR results
print("||===================================================================================")
eSMR <- read_smr_data1(SMR_DIRT,result,trait_name,qtl_type="eQTL")
sSMR <- read_smr_data1(SMR_DIRT,result,trait_name,qtl_type="sQTL")
pSMR <- read_smr_data1(SMR_DIRT,result,trait_name,qtl_type="pQTL")
edSMR <- read_smr_data1(SMR_DIRT,result,trait_name,qtl_type="edQTL") 

mSMR <- read_smr_data2(SMR_DIRT,result,trait_name,qtl_type="mQTL",QTL_link=CpG_link)
hSMR <- read_smr_data2(SMR_DIRT,result,trait_name,qtl_type="hQTL",QTL_link=hQTL_link)
caSMR <- read_smr_data2(SMR_DIRT,result,trait_name,qtl_type="caQTL",QTL_link=caQTL_link)

# select QTL list for MAGIC analysis
select_SMR_results_in_QTL_name_list <- function(xSMR, QTL_name_list) {
    SMR_p_ACAT <- xSMR$SMR_p_ACAT
    SMR_probeID <- xSMR$SMR_probeID
    SMR_results <- xSMR$SMR_results
    SMR_p_ACAT_tmp = SMR_p_ACAT[, which(colnames(SMR_p_ACAT) %in% QTL_name_list$V1), drop = FALSE]
    SMR_probeID_tmp = SMR_probeID[, which(colnames(SMR_probeID) %in% paste0("probeID_", QTL_name_list$V1)), drop = FALSE]
    SMR_results_tmp = SMR_results[which(SMR_results$QTL_name %in% QTL_name_list$V1), , drop = FALSE]
    xSMR$SMR_p_ACAT = SMR_p_ACAT_tmp
    xSMR$SMR_probeID = SMR_probeID_tmp
    xSMR$SMR_results = SMR_results_tmp
    return(xSMR)
}

QTL_name_list = fread(QTL_name_list_file, header = FALSE)

eSMR <- select_SMR_results_in_QTL_name_list(eSMR, QTL_name_list)
sSMR <- select_SMR_results_in_QTL_name_list(sSMR, QTL_name_list)
pSMR <- select_SMR_results_in_QTL_name_list(pSMR, QTL_name_list)
edSMR <- select_SMR_results_in_QTL_name_list(edSMR, QTL_name_list)
mSMR <- select_SMR_results_in_QTL_name_list(mSMR, QTL_name_list)
hSMR <- select_SMR_results_in_QTL_name_list(hSMR, QTL_name_list)
caSMR <- select_SMR_results_in_QTL_name_list(caSMR, QTL_name_list)

####################################################################################################
# Only keep within locus texts
####################################################################################################

SMR_results=rbind(eSMR$SMR_results, sSMR$SMR_results, pSMR$SMR_results, edSMR$SMR_results,
                  mSMR$SMR_results, hSMR$SMR_results, caSMR$SMR_results)
SMR_results$GWAS_LOCUS=NA
SMR_results$Lead_SNP=NA
SMR_results$Lead_SNP_BP=NA
for(j in 1:nrow(Locus_data)){
  # print(j)
  chr=sub("^chr", "", Locus_data$chr[j])
  start=Locus_data$start[j]
  end=Locus_data$end[j]
  locus=Locus_data$GWAS_LOCUS[j]
  lead_snp=Locus_data$Lead_SNP[j]
  lead_snp_bp=Locus_data$Lead_SNP_BP[j]

  index_tmp=which(SMR_results$ProbeChr==chr & SMR_results$Probe_bp<=end & SMR_results$Probe_bp>=start)
  SMR_results$GWAS_LOCUS[index_tmp]=locus
  SMR_results$Lead_SNP[index_tmp]=lead_snp
  SMR_results$Lead_SNP_BP[index_tmp]=lead_snp_bp
}

SMR_results_within_locus=SMR_results[which(!is.na(SMR_results$GWAS_LOCUS)),]
write.table(SMR_results_within_locus,paste0(OUTPUT,"/MAGIC/results/",trait_name,"_MAGIC_within_locus.txt"),row=F,col=T,quo=F,sep="\t")



####################################################################################################
# MAGIC calculation for each category
####################################################################################################
p_ACAT_eSMR=apply(eSMR$SMR_p_ACAT,1,ACAT)
p_ACAT_sSMR=apply(sSMR$SMR_p_ACAT,1,ACAT)
p_ACAT_pSMR=apply(pSMR$SMR_p_ACAT,1,ACAT)
p_ACAT_edSMR=apply(edSMR$SMR_p_ACAT,1,ACAT)
p_ACAT_mSMR=apply(mSMR$SMR_p_ACAT,1,ACAT)
p_ACAT_hSMR=apply(hSMR$SMR_p_ACAT,1,ACAT)
p_ACAT_caSMR=apply(caSMR$SMR_p_ACAT,1,ACAT)
p_ACAT=apply(cbind(eSMR$SMR_p_ACAT, sSMR$SMR_p_ACAT, pSMR$SMR_p_ACAT, edSMR$SMR_p_ACAT,mSMR$SMR_p_ACAT, hSMR$SMR_p_ACAT, caSMR$SMR_p_ACAT),1,ACAT)

ACAT_results=data.frame(p_ACAT,p_ACAT_eSMR,p_ACAT_sSMR,p_ACAT_pSMR,p_ACAT_edSMR,p_ACAT_mSMR,p_ACAT_hSMR,p_ACAT_caSMR)
# ACAT_results$gene_name=names(p_ACAT)
MAGIC_results=cbind(result, ACAT_results)
colnames(MAGIC_results) = c("chr", "start", "end", "strand", "gene_id", "gene_name", "GWAS_LOCUS", "Lead_SNP", "Lead_SNP_BP",
                            "MAGIC","eMAGIC", "sMAGIC", "pMAGIC", "edMAGIC", "mMAGIC", "hMAGIC", "caMAGIC")
# write.table(MAGIC_results,paste0(OUTPUT,"/MAGIC/summary/",trait_name,"_MAGIC.txt"),row=F,col=T,quo=F,sep="\t")


####################################################################################################
# MAGIC plots
####################################################################################################
# index <- which(!is.na(MAGIC_results$GWAS_LOCUS) & MAGIC_results$MAGIC < 0.05 / length(MAGIC_results$MAGIC))
index <- which(!is.na(MAGIC_results$GWAS_LOCUS))
MAGIC_plot <- MAGIC_results[index, ]

eSMR_plot <- cbind(result[index, , drop=FALSE], eSMR$SMR_p_ACAT[index, , drop=FALSE], eSMR$SMR_probeID[index, , drop=FALSE])
sSMR_plot <- cbind(result[index, , drop=FALSE], sSMR$SMR_p_ACAT[index, , drop=FALSE], sSMR$SMR_probeID[index, , drop=FALSE])
pSMR_plot <- cbind(result[index, , drop=FALSE], pSMR$SMR_p_ACAT[index, , drop=FALSE], pSMR$SMR_probeID[index, , drop=FALSE])
edSMR_plot <- cbind(result[index, , drop=FALSE], edSMR$SMR_p_ACAT[index, , drop=FALSE], edSMR$SMR_probeID[index, , drop=FALSE])
mSMR_plot <- cbind(result[index, , drop=FALSE], mSMR$SMR_p_ACAT[index, , drop=FALSE], mSMR$SMR_probeID[index, , drop=FALSE])
hSMR_plot <- cbind(result[index, , drop=FALSE], hSMR$SMR_p_ACAT[index, , drop=FALSE], hSMR$SMR_probeID[index, , drop=FALSE])
caSMR_plot <- cbind(result[index, , drop=FALSE], caSMR$SMR_p_ACAT[index, , drop=FALSE], caSMR$SMR_probeID[index, , drop=FALSE])

write.table(eSMR_plot,paste0(OUTPUT,"/MAGIC/plot/",trait_name,"_eSMR.summary"),row=F,col=T,quo=F,sep="\t")
write.table(sSMR_plot,paste0(OUTPUT,"/MAGIC/plot/",trait_name,"_sSMR.summary"),row=F,col=T,quo=F,sep="\t")
write.table(pSMR_plot,paste0(OUTPUT,"/MAGIC/plot/",trait_name,"_pSMR.summary"),row=F,col=T,quo=F,sep="\t")
write.table(edSMR_plot,paste0(OUTPUT,"/MAGIC/plot/",trait_name,"_edSMR.summary"),row=F,col=T,quo=F,sep="\t")
write.table(mSMR_plot,paste0(OUTPUT,"/MAGIC/plot/",trait_name,"_mSMR.summary"),row=F,col=T,quo=F,sep="\t")
write.table(hSMR_plot,paste0(OUTPUT,"/MAGIC/plot/",trait_name,"_hSMR.summary"),row=F,col=T,quo=F,sep="\t")
write.table(caSMR_plot,paste0(OUTPUT,"/MAGIC/plot/",trait_name,"_caSMR.summary"),row=F,col=T,quo=F,sep="\t")


# only genes in locus ----
QTL_types <- c("eSMR", "sSMR", "pSMR", "edSMR", "mSMR", "hSMR", "caSMR")
MAGIC_plot[paste0(QTL_types, "_name")] <- NA
MAGIC_plot[paste0(QTL_types, "_probeID")] <- NA


get_min_index_and_name <- function(gene_name, SMR_p_ACAT, SMR_probeID) {
    if (!gene_name %in% rownames(SMR_p_ACAT)) {
        return(list(name = NA, probeID = NA))
    }

    row_data <- SMR_p_ACAT[gene_name, , drop = FALSE]
    if (ncol(row_data) == 0) {
        return(list(name = NA, probeID = NA))
    }
    index <- which.min(row_data)
    name <- if(length(index) > 0 && !is.na(index)) colnames(SMR_p_ACAT)[index] else NA
    probeID <- if(length(index) > 0 && !is.na(index)) SMR_probeID[gene_name, index] else NA
    return(list(name = name, probeID = probeID))
}

for (qtl in QTL_types) {
    if (is.null(dimnames(get(qtl)$SMR_p_ACAT))) {
      MAGIC_plot[[paste0(qtl, "_name")]] <- NA
      MAGIC_plot[[paste0(qtl, "_probeID")]] <- NA
    }else{
      
    MAGIC_plot_list <- apply(MAGIC_plot[, "gene_name", drop = FALSE], 1, function(gene_name) {
        get_min_index_and_name(gene_name, get(qtl)$SMR_p_ACAT, get(qtl)$SMR_probeID)
    })
    MAGIC_plot[[paste0(qtl, "_name")]] <- sapply(MAGIC_plot_list, function(x) x$name)
    MAGIC_plot[[paste0(qtl, "_probeID")]] <- sapply(MAGIC_plot_list, function(x) x$probeID)
    }
}

colnames(MAGIC_plot) = c("chr", "start", "end", "strand", "gene_id", "gene_name", "GWAS_LOCUS", "Lead_SNP", "Lead_SNP_BP",
                            "MAGIC","eMAGIC", "sMAGIC", "pMAGIC", "edMAGIC","mMAGIC", "hMAGIC", "caMAGIC",
                            "eMAGIC_QTL_name", "sMAGIC_QTL_name", "pMAGIC_QTL_name", "edMAGIC_QTL_name", "mMAGIC_QTL_name", "hMAGIC_QTL_name", "caMAGIC_QTL_name",
                            "eMAGIC_probeID", "sMAGIC_probeID", "pMAGIC_probeID", "edMAGIC_probeID","mMAGIC_probeID", "hMAGIC_probeID", "caMAGIC_probeID")
write.table(MAGIC_plot,paste0(OUTPUT,"/MAGIC/plot/",trait_name,"_MAGIC_plot.summary"),row=F,col=T,quo=F,sep="\t")


# MAGIC results ----
# all genes 
QTL_types <- c("eSMR", "sSMR", "pSMR", "edSMR","mSMR", "hSMR", "caSMR")
MAGIC_results[paste0(QTL_types, "_name")] <- NA
MAGIC_results[paste0(QTL_types, "_probeID")] <- NA

get_min_index_and_name <- function(gene_name, SMR_p_ACAT, SMR_probeID) {
    if (!gene_name %in% rownames(SMR_p_ACAT)) {
        return(list(name = NA, probeID = NA))
    }

    row_data <- SMR_p_ACAT[gene_name, , drop = FALSE]
    if (ncol(row_data) == 0) {
        return(list(name = NA, probeID = NA))
    }
    index <- which.min(row_data)
    name <- if(length(index) > 0 && !is.na(index)) colnames(SMR_p_ACAT)[index] else NA
    probeID <- if(length(index) > 0 && !is.na(index)) SMR_probeID[gene_name, index] else NA
    return(list(name = name, probeID = probeID))
}

for (qtl in QTL_types) {
    if (is.null(dimnames(get(qtl)$SMR_p_ACAT))) {
      MAGIC_results[[paste0(qtl, "_name")]] <- NA
      MAGIC_results[[paste0(qtl, "_probeID")]] <- NA
    }else{
      
    MAGIC_results_list <- apply(MAGIC_results[, "gene_name", drop = FALSE], 1, function(gene_name) {
        get_min_index_and_name(gene_name, get(qtl)$SMR_p_ACAT, get(qtl)$SMR_probeID)
    })
    MAGIC_results[[paste0(qtl, "_name")]] <- sapply(MAGIC_results_list, function(x) x$name)
    MAGIC_results[[paste0(qtl, "_probeID")]] <- sapply(MAGIC_results_list, function(x) x$probeID)
    }
}



colnames(MAGIC_results) = c("chr", "start", "end", "strand", "gene_id", "gene_name", "GWAS_LOCUS", "Lead_SNP", "Lead_SNP_BP",
                            "MAGIC","eMAGIC", "sMAGIC", "pMAGIC", "edMAGIC","mMAGIC", "hMAGIC", "caMAGIC",
                            "eMAGIC_QTL_name", "sMAGIC_QTL_name", "pMAGIC_QTL_name", "edMAGIC_QTL_name","mMAGIC_QTL_name", "hMAGIC_QTL_name", "caMAGIC_QTL_name",
                            "eMAGIC_probeID", "sMAGIC_probeID", "pMAGIC_probeID", "edMAGIC_probeID","mMAGIC_probeID", "hMAGIC_probeID", "caMAGIC_probeID")

MAGIC_results=MAGIC_results[which(!is.na(MAGIC_results$MAGIC)),]
MAGIC_results=MAGIC_results[order(MAGIC_results$MAGIC), ]

write.table(MAGIC_results,paste0(OUTPUT,"/MAGIC/summary/",trait_name,"_MAGIC.txt"),row=F,col=T,quo=F,sep="\t")
