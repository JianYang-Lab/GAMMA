args <- commandArgs(trailingOnly = TRUE)
trait_name=args[1]
OUTPUT=args[2]
window_size=as.numeric(args[3])
gencode_file=args[4]

suppressMessages({
  library(GenomicRanges)
  library(regioneR)
  library(Repitools)
  library(data.table)
})

COJO=read.table(paste0(OUTPUT,"/COJO/summary/",trait_name,".jma.cojo"),head=T,stringsAsFactors=F)

COJO$Chr=paste0("chr", COJO$Chr)
COJO$LOCUS_START=COJO$bp-window_size
COJO$LOCUS_END=COJO$bp+window_size
region <- GRanges(COJO$Chr, IRanges(start = COJO$LOCUS_START, end = COJO$LOCUS_END))
merges=mergeRegions(region,region)
res=annoGR2DF(merges)

index=which(res$width >= 3*window_size)

if(length(index) > 0){
res_tmp=res[index,]
res_need=res[-index,]

res_tmp_split=data.frame()

for(region_index in 1:nrow(res_tmp)) {
  # print(paste0("region_index",region_index))
  chr_region = res_tmp$chr[region_index]
  start_region = res_tmp$start[region_index]
  end_region = res_tmp$end[region_index]

  # Subset COJO for SNPs within the current region
  snp_indices = which(COJO$Chr == chr_region & COJO$bp >= start_region & COJO$bp <= end_region)
  COJO_tmp = COJO[snp_indices, ]
  COJO_tmp = COJO_tmp[order(COJO_tmp$bp), ]  # Order by bp

# ---------------------------------------------------------------
# Make all the regions <= 3 Mb
# Rescale the LOCUS_START and LOCUS_END part
# The key here is that the distance between neighbor SNPs > 1Mb are separated
# < 1Mb are in the same locus.
#  ---------------------------------------------------------------

  # Initialize vectors for adjusted LOCUS_START and LOCUS_END
  adjusted_LOCUS_START = numeric(nrow(COJO_tmp))
  adjusted_LOCUS_END = numeric(nrow(COJO_tmp))
  adjusted_LOCUS_START[1] = COJO_tmp$LOCUS_START[1]
  
  # Adjust LOCUS_START and LOCUS_END based on SNP distances
  for(snp_index in 1:nrow(COJO_tmp)) {
    # print(paste0("snp_index",snp_index))
    if(snp_index < nrow(COJO_tmp)) {
      diff_bp = COJO_tmp$bp[snp_index + 1] - COJO_tmp$bp[snp_index]
      if(diff_bp > window_size) {
        mid_point = round((COJO_tmp$bp[snp_index] + COJO_tmp$bp[snp_index + 1]) / 2)
        adjusted_LOCUS_END[snp_index] = mid_point
        adjusted_LOCUS_START[snp_index + 1] = mid_point
      } else {
        adjusted_LOCUS_START[snp_index + 1] = adjusted_LOCUS_START[snp_index]
      }
    } else {
      adjusted_LOCUS_END[snp_index] = COJO_tmp$LOCUS_END[snp_index]
    }
  }
  # Ensure continuity for LOCUS_END in segments with adjacent SNPs closer than 1Mb
  for(snp_index in 1:(nrow(COJO_tmp) - 1)) {
    # print(paste0("snp_index_new",snp_index))
    if(adjusted_LOCUS_START[snp_index] == adjusted_LOCUS_START[snp_index + 1]) {
      last_index = max(which(adjusted_LOCUS_START == adjusted_LOCUS_START[snp_index]))
      adjusted_LOCUS_END[snp_index] = adjusted_LOCUS_END[last_index]
    }
  }

  # Update COJO_tmp with adjusted start/end
  COJO_tmp$LOCUS_START = adjusted_LOCUS_START
  COJO_tmp$LOCUS_END = adjusted_LOCUS_END

  # Aggregate and combine results
  res_tmp_split_tmp = unique(data.frame(chr = COJO_tmp$Chr, start = COJO_tmp$LOCUS_START, end = COJO_tmp$LOCUS_END, width = COJO_tmp$LOCUS_END - COJO_tmp$LOCUS_START))
  res_tmp_split = rbind(res_tmp_split, res_tmp_split_tmp)
}
res=rbind(res_need, res_tmp_split)
}

index=which(res$start < 0)
res$start[index]=0
res=res[order(res$chr, res$start), ]

res$GWAS_LOCUS=paste0(res$chr,":",res$start,":",res$end)
res$Lead_SNP=NA
res$Lead_SNP_BP=NA

# ----------------------------------------------------------------------------
# Once the locus are defined, we need to define the lead SNP within each locus
#  ---------------------------------------------------------------------------
for(j in 1:nrow(res)){
  chr=res$chr[j]
  start=res$start[j]
  end=res$end[j]  
  
  k=which(COJO$Chr==chr & COJO$bp <= end & COJO$bp >= start)
  if(length(k)>1){
    k_index=which.min(COJO$p[k]) 
    res$Lead_SNP[j]=COJO$SNP[k[k_index]]
    res$Lead_SNP_BP[j]=COJO$bp[k[k_index]]
  }else{
    res$Lead_SNP[j]=COJO$SNP[k]
    res$Lead_SNP_BP[j]=COJO$bp[k]
  }    
}

write.table(res, paste0(OUTPUT,"/COJO/summary/",trait_name,"_",format(2*window_size, scientific = FALSE),".locus"),row=F,col=T,quo=F,sep="\t")


### ******************************************************************************************** 
# genecode summary
## ********************************************************************************************
# 0. loading gene annotation
# gencode_file="/storage/yangjianLab/qiting/data/annotation/gencode/gencode.v40.GRCh38.gene.annotation.bed"
gencode=fread(gencode_file,head=T,stringsAsFactors=F,data.table=F)
gencode$V1 <- sub("^chr", "", gencode$V1)
gencode$gene_id <- sub("\\..*", "", gencode$gene_id)

index=which(gencode$gene_type=="protein_coding" & gencode$V1!="M" & gencode$V1!="X" & gencode$V1!="Y")
gencode1=gencode[index,]
gencode1$gene_tss = ifelse(gencode1$V7 == "+", gencode1$V4, gencode1$V5)

result=data.frame(gene_id=gencode1$gene_id,gene_name=gencode1$gene_name,chr=gencode1$V1,start=gencode1$V4,end=gencode1$V5,strand=gencode1$V7)
result$GWAS_LOCUS=NA;result$Lead_SNP=NA;result$Lead_SNP_BP=NA;

# 1. COJO locus 
COJO_locus=res

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

res$P=NA
index=match(res$Lead_SNP, COJO$SNP, nomatch=0)
res$P[which(index!=0)]=COJO$p[index]

write.table(result, paste0(OUTPUT,"/COJO/summary/",trait_name,"_",format(2*window_size, scientific = FALSE),"_gene.locus"),row=F,col=T,quo=F,sep="\t")
