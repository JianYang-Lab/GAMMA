#### ------------------------
#### CARMA analysis ####
#### ------------------------
args <- commandArgs(trailingOnly = TRUE)
GWAS_DATA=args[1]
trait_name=args[2]
locus_num=as.numeric(args[3])
OUTPUT=args[4]
locus_path=args[5]
REFERENCE=args[6]
PLINK_1.9=args[7]
PLINK_2=args[8]

# PLINK_1.9="/storage/yangjianLab/guoyazhou/software/plink/plink1.90beta/plink_1.90_beta"
# PLINK_2="/storage/yangjianLab/guoyazhou/software/plink/plink2.00alpha/plink2"

suppressMessages({
	library(data.table)
	library(magrittr)
	library(dplyr)
	library(devtools)
	library(R.utils)
	library(CARMA)
	library(data.table)
})

GWAS=fread(GWAS_DATA, header=T,col.names=c("SNP","A1","A2","AF1","BETA","SE","P","N"))
GWAS$P=as.numeric(GWAS$P)
index=which(GWAS$P == 0)
GWAS$P[index]=1e-300

# to make the LD matrix smaller
# GWAS=GWAS[which(GWAS$P<0.0001),]
GWAS=GWAS[which(GWAS$P<0.01),]

GWAS$A1=toupper(GWAS$A1)
GWAS$A2=toupper(GWAS$A2)

GWAS$Z = GWAS$BETA / GWAS$SE

index=which(is.na(GWAS$Z))
if(length(index) > 0){
	GWAS$Z[index] = ifelse(GWAS$BETA[index] > 0,
						-qnorm(GWAS$P[index] / 2),
						qnorm(GWAS$P[index] / 2))									
	# z_abs = -qnorm(p / 2)
	# z = ifelse(beta > 0, z_abs, -z_abs)
}

COJO_locus=fread(locus_path, header=T)
# ------------------------
# Run CARMA analysis by COJO independent SNP (locus)
# ------------------------
chr = gsub("chr","",COJO_locus$chr[locus_num])
start = COJO_locus$start[locus_num]
end = COJO_locus$end[locus_num]
locus = COJO_locus$GWAS_LOCUS[locus_num]


# Reference ----------------
GRCh38_reference_bim=fread(paste0(REFERENCE,"_chr",chr,".bim"), col.names=c("CHR","SNP","MOL","POS","A1","A2"))
Reference=paste0(REFERENCE, "_chr",chr)

GWAS$CHR=NA
GWAS$POS=NA
index=match(GWAS$SNP, GRCh38_reference_bim$SNP, nomatch=0)
GWAS$CHR[which(index!=0)]=GRCh38_reference_bim$CHR[index]
GWAS$POS[which(index!=0)]=GRCh38_reference_bim$POS[index]

index=which(GRCh38_reference_bim$CHR == chr & GRCh38_reference_bim$POS >= start & GRCh38_reference_bim$POS <= end)
reference_locus=GRCh38_reference_bim[index,]
commSNP = intersect(GWAS$SNP, reference_locus$SNP)

if (length(commSNP) <= 5) {
	print(paste0("Not enough SNPs in GWAS locus ",locus,"!"))
	quit(save = "no", status = 0)
}

### Keep the SNP order as refFile
commSum <- GWAS[match(commSNP,GWAS$SNP),]
commRef <- reference_locus[match(commSNP,reference_locus$SNP),]
	

## Flip alleles to make consistent sumFile & LD
flipIndex <- which(commSum$A1 != commRef$A1)
if(length(flipIndex)>0){
A0 <- commSum$A1[flipIndex]
commSum$A1[flipIndex] <- commSum$A2[flipIndex]
commSum$A2[flipIndex] <- A0

commSum$Z[flipIndex] <- -commSum$Z[flipIndex]
commSum$BETA[flipIndex] <- -commSum$BETA[flipIndex]
commSum$AF1[flipIndex] <- 1-commSum$AF1[flipIndex]
}

keepIndex <- which(commSum$A1 == commRef$A1 & commSum$A2 == commRef$A2)
commSum <- commSum[keepIndex,]
commRef <- commRef[keepIndex,]

fwrite(as.data.frame(commRef$SNP),paste0(OUTPUT,"/CARMA/LD/",locus,".snplist"),quote=F,col.names=F,row.names=F,na="NA")


SNP_FILE=paste0(OUTPUT,"/CARMA/LD/",locus,".snplist")
LD_FILE=paste0(OUTPUT,"/CARMA/LD/",locus)

system(paste0(PLINK_2," --bfile ", Reference, " --extract ", SNP_FILE, " --make-bed --out ", LD_FILE))
system(paste0(PLINK_1.9," --bfile ", LD_FILE," --r square spaces --keep-allele-order --threads 12 --out ",LD_FILE))
# command=paste0(PLINK_1.9," --bfile ",Reference," --extract ",SNP_FILE," --r square spaces --out ",outdirt,"LD/chr",Chr,"_",bp_38)

## Get final inputs
bimFile <- fread(paste0(LD_FILE,".bim"),header=F,col.names=c("CHR","SNP","MOL","POS","A1","A2"))
ld <- as.matrix(fread(paste0(LD_FILE,".ld"),header=F))

# Note set NA or NaN to 0 ---------------------------------------------
ld[is.na(ld)] <- 0
diag(ld) <- 1

finalIndex <- which(complete.cases(ld))
if (length(finalIndex) <= 5) {
	print(paste0("Not enough SNP in GWAS locus: ",locus,"!"))
	quit(save = "no", status = 0)
}

ld <- ld[finalIndex,finalIndex]
refFinal <- bimFile[finalIndex,]
sumFinal <- commSum[match(refFinal$SNP,commSum$SNP),]



## -----------------------------------------------------------
## Run CARMA -------------------------------------------------
## -----------------------------------------------------------
setwd(paste0(OUTPUT,"/CARMA/tmp/"))

z.list<-list()
ld.list<-list()
lambda.list<-list()
annot.list<-list()

z.list[[1]]<-sumFinal$Z
ld.list[[1]]<-ld
# annot.list[[1]]<- as.matrix(annotFinal[,-c(1:4)])
lambda.list[[1]]<-1

CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,outlier.switch=T)
# CARMA.results<-CARMA(z.list,ld.list,lambda.list=lambda.list,w.list=annot.list,outlier.switch=T)

sumstat.result = sumFinal %>% mutate(PIP = CARMA.results[[1]]$PIPs, CS = 0)
if(length(CARMA.results[[1]]$`Credible set`[[2]])!=0){
	for(l in 1:length(CARMA.results[[1]]$`Credible set`[[2]])){
	sumstat.result$CS[CARMA.results[[1]]$`Credible set`[[2]][[l]]]=l
	}
}

warnings()
print(paste0(OUTPUT,"/CARMA/Result/",trait_name,"_",locus,".txt.gz"))
fwrite(sumstat.result,paste0(OUTPUT,"/CARMA/Result/",trait_name,"_",locus,".txt.gz"),sep = "\t", quote = F, na = "NA", row.names = F, col.names = T, compress = "gzip")

system(paste0("rm ", LD_FILE, ".*"))
