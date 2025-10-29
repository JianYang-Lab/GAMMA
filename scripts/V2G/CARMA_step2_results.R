#### ------------------------
#### CARMA analysis ####
#### ------------------------
args <- commandArgs(trailingOnly = TRUE)
trait_name=args[1]
OUTPUT=args[2]

suppressMessages({
	library(data.table)
})

# CARMA results ----------------------------------------
CARMA_result_dir=paste0(OUTPUT,"/CARMA/Result/")
CARMA_CS_dir=paste0(OUTPUT,"/CARMA/summary")

CARMA_file_list=list.files(CARMA_result_dir)
CARMA_files=CARMA_file_list[grep(paste0("^",trait_name,"_chr"),CARMA_file_list)]

CS_dataset=data.frame()
bed_dataset=data.frame()

if(length(CARMA_files)>0){

	for(file_index in CARMA_files){
	CARMA_result=fread(paste0(CARMA_result_dir,file_index))
	index=which(CARMA_result$CS>=1 & CARMA_result$PIP > 0.1 & CARMA_result$P < 5e-8)
	# index=which(CARMA_result$CS>=1)
	CS=CARMA_result[index,]

	if(nrow(CS)>0){
	CS_dataset=rbind(CS_dataset,CS)
	}
	}

if(nrow(CS_dataset) >0){
bed_dataset=data.frame(chr=paste0("chr",CS_dataset$CHR),start=CS_dataset$POS,end=CS_dataset$POS+1,SNP=CS_dataset$SNP,PIP=CS_dataset$PIP)
}
}

CS_dataset$POS <- as.integer(CS_dataset$POS)
bed_dataset$start <- as.integer(bed_dataset$start)
bed_dataset$end <- as.integer(bed_dataset$end)

write.table(CS_dataset,paste0(CARMA_CS_dir,"/",trait_name,"_CARMA.CS"),row=F,col=T,quo=F,sep="\t")
write.table(bed_dataset,paste0(CARMA_CS_dir,"/",trait_name,"_CARMA.bed"),row=F,col=F,quo=F,sep="\t")

