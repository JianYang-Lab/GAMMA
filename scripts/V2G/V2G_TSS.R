args <- commandArgs(trailingOnly = TRUE)

trait_name=args[1]
OUTPUT=args[2]
gencode_tss=args[3]
CARMA_bed_file=args[4]

# ------------------------------------------------
# CARMA fine-mapping cloest gene -----------------
# ------------------------------------------------
suppressMessages({
library(data.table)
})

gene_bed=fread(gencode_tss,head=F,stringsAsFactors=F,data.table=F)
if(is.null(CARMA_bed_file)){
  CARMA_bed_file=paste0(OUTPUT,"/CARMA/summary/",trait_name,"_CARMA.bed")
}

INFILE=CARMA_bed_file

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
    data <- fread(INFILE, header = F, stringsAsFactors = FALSE, data.table = FALSE)
  }
}



if(nrow(data)>0){
  data$TSS_Gene=NA;
for(j in 1:nrow(data)){
  # print(j)
  index=which(gene_bed$V1 == data$V1[j])
  bed_tmp=gene_bed[index,]
  # print(data$V2[j])
  # print(bed_tmp$V2)
  dis_TSS=abs(data$V2[j]-bed_tmp$V2)
  data$TSS_Gene[j]=bed_tmp$V4[which.min(dis_TSS)]
}
}

write.table(data,paste0(OUTPUT,"/V2G/summary/",trait_name,"_closestTSS.annot.summary"),row=F,col=F,quo=F,sep="\t")
