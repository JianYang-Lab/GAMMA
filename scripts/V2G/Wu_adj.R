args <- commandArgs(trailingOnly = TRUE)
GWAS_DATA=args[1]
trait_name=args[2]
OUTPUT=args[3]
locus_path=args[4]
PLINK=args[5]
reference_all_bim_file=args[6]
reference_freq_file=args[7]
reference_bfile=args[8]

# ------------------------------------------------------------
suppressMessages({
	library(data.table)
})

# Reference information --------------------------------------
reference_bim=fread(reference_all_bim_file, col.names=c("CHR","SNP","MOL","POS","A1","A2"))
reference_freq=fread(reference_freq_file,head=T,stringsAsFactors=F,data.table=F)
reference_bim$AF1=reference_bim$MAF=NA
index1=match(paste0(reference_bim$SNP,"_",reference_bim$A1), paste0(reference_freq$ID,"_",reference_freq$ALT), nomatch=0)
reference_bim$AF1[which(index1!=0)]=reference_freq$ALT_FREQS[index1]
reference_bim$MAF=ifelse(reference_bim$AF1>0.5,1-reference_bim$AF1,reference_bim$AF1)

# GWAS data ---------------------------------------------------
GWAS=fread(GWAS_DATA, header=T,col.names=c("SNP","A1","A2","AF1","BETA","SE","P","N"))
GWAS$P=as.numeric(GWAS$P)
index=which(GWAS$P == 0)
GWAS$P[index]=1e-300
GWAS$A1=toupper(GWAS$A1)
GWAS$A2=toupper(GWAS$A2)


GWAS$CHR=GWAS$POS=NA
index=match(GWAS$SNP, reference_bim$SNP, nomatch=0)
GWAS$CHR[which(index!=0)]=reference_bim$CHR[index]
GWAS$POS[which(index!=0)]=reference_bim$POS[index]

GWAS$Chisq <- (GWAS$BETA/GWAS$SE)^2
GWAS$MAF = ifelse(GWAS$AF1 > 0.5,1-GWAS$AF1,GWAS$AF1)


# COJO or Clumping results for Wu_adj method.
# locus_data=read.table(COJO_or_Clumping_file, head=T,stringsAsFactors=F)
# colnames(locus_data)=toupper(colnames(locus_data))

# locus_path="/storage/yangjianLab/guoyazhou/GAMMA_git_data/COJO/summary/T2D.locus"
locus_data=fread(locus_path, header=T)



CS_dataset=data.frame()
bed_dataset=data.frame()

for(j in 1:nrow(locus_data)){
	print(j)
	
	SNP=locus_data$Lead_SNP[j]
	CHR=locus_data$chr[j]
	CHR=gsub("chr","",CHR)
	BP=locus_data$Lead_SNP_BP[j]

	start=locus_data$start[j]
	end=locus_data$end[j]
	gwas_locus=locus_data$GWAS_LOCUS[j]

	index=which(GWAS$CHR==CHR & GWAS$POS<=end & GWAS$POS>=start)
	gwas_data=GWAS[index,]
	gwas_data$R2=NA
	gwas_data$CS=NA


	index=which(reference_bim$SNP == SNP)
	MAF=reference_bim$MAF[index]

	# 5kb, MAF 0.01, LD>=0.9, GWAS p<5e-8
	# 34kb, MAF 0.05, LD>=0.8
	# index=which(reference_bim$CHR==CHR & reference_bim$POS<=(BP+5000) & reference_bim$POS>=(BP-5000) & reference_bim$MAF<=(MAF+0.01) & reference_bim$MAF>=(MAF-0.01))
	index=which(reference_bim$CHR==CHR & reference_bim$POS<=(BP+34000) & reference_bim$POS>=(BP-34000) & reference_bim$MAF<=(MAF+0.05) & reference_bim$MAF>=(MAF-0.05))	

	if(length(index)>0){
		
		SNP_list=reference_bim[index,]
		SNP_file=paste0(OUTPUT,"/Wu_adj/LD/",trait_name,"_",SNP,"_",CHR,"_",BP,".snp.list")
		write.table(SNP_list$SNP,SNP_file,row=F,col=F,quo=F,sep="\t")

		command=paste0(PLINK," --bfile ",reference_bfile,"_chr",CHR," --extract ",SNP_file," --r square spaces --out ", OUTPUT,"/Wu_adj/LD/",trait_name,"_",SNP,"_",CHR,"_",BP)
		system(command)
		Ld_matrix <- fread(paste0(OUTPUT,"/Wu_adj/LD/",trait_name,"_",SNP,"_",CHR,"_",BP,".ld"), header=F)

		SNP_list$R2=as.numeric(Ld_matrix[which(SNP_list$SNP==SNP),])^2
		index=which(SNP_list$R2>=0.8)

		if(length(index)>0){
			CS=SNP_list[index,]
			bed=data.frame(chr=paste0("chr",CS$CHR),start=CS$POS,end=CS$POS+1,SNP=CS$SNP,R2=CS$R2)
			
			CS_dataset=rbind(CS_dataset,CS)
			bed_dataset=rbind(bed_dataset,bed)
			
			write.table(CS, paste0(OUTPUT,"/Wu_adj/Result/",trait_name,"_",SNP,"_",CHR,"_",BP,".CS"),row=F,col=T,quo=F,sep="\t")
			write.table(bed,paste0(OUTPUT,"/Wu_adj/Result/",trait_name,"_",SNP,"_",CHR,"_",BP,".bed"),row=F,col=F,quo=F,sep="\t")

			index=match(gwas_data$SNP, CS$SNP, nomatch=0)
			gwas_data$R2[which(index!=0)]=CS$R2[index]
			gwas_data$CS[which(index!=0)]=1
		}	

	
	file.remove(paste0(OUTPUT,"/Wu_adj/LD/",trait_name,"_",SNP,"_",CHR,"_",BP,".snp.list"))
	file.remove(paste0(OUTPUT,"/Wu_adj/LD/",trait_name,"_",SNP,"_",CHR,"_",BP,".ld"))
	file.remove(paste0(OUTPUT,"/Wu_adj/LD/",trait_name,"_",SNP,"_",CHR,"_",BP,".log"))
	}
	# print(paste0(OUTPUT,"/Wu_adj/Result/",trait_name,"_",gwas_locus,".txt.gz"))
	# fwrite(gwas_data, paste0(OUTPUT,"/Wu_adj/Result/",trait_name,"_",gwas_locus,".txt.gz"),sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")
	write.table(gwas_data, paste0(OUTPUT,"/Wu_adj/Result/",trait_name,"_",gwas_locus,".txt"),row=F,col=T,quo=F,sep="\t")
}


index=match(GWAS$SNP, CS_dataset$SNP, nomatch=0)
GWAS$R2=NA
GWAS$R2[which(index!=0)]=CS_dataset$R2[index]
GWAS_Wu_adj=GWAS[index,]

CS_dataset=CS_dataset[index,]

index=match(GWAS$SNP, bed_dataset$SNP, nomatch=0)
bed_dataset=bed_dataset[index,]

write.table(CS_dataset, paste0(OUTPUT,"/Wu_adj/summary/",trait_name,".CS"),row=F,col=T,quo=F,sep="\t")
write.table(bed_dataset,paste0(OUTPUT,"/Wu_adj/summary/",trait_name,".bed"),row=F,col=F,quo=F,sep="\t")
write.table(GWAS_Wu_adj,paste0(OUTPUT,"/Wu_adj/summary/",trait_name,".gwas"),row=F,col=T,quo=F,sep="\t")
