##Title:trans qtl discovery
################################

#!/usr/bin/env Rscript --vanilla
args = commandArgs(trailingOnly=TRUE)
gene = args[1]
## select working directory
setwd("/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/trans_qtl_discover")
inputdir<-"/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/trans_qtl_discover"

## create outdirectory:
outdir<-("/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/trans_qtl_discover/trans_output/")

## load the R matrix for expression data
load("/dcs04/lieber/statsgen/sjahansi/voomSVA.RData",verbose=T)
## make a variable for gene expression dataset
expr <- v$E

## prepare expression file
## match the sample_ids with Pheno_pc file
Pheno_PC = read.delim("/dcs04/lieber/statsgen/shizhong/database/libd/genotype/postmortem/phenotype/pheno_PC", header = TRUE)
targets <- v$targets
sample_id <- match(targets$BrNum,Pheno_PC$BrNum)
sample_id <- Pheno_PC[sample_id,]
## then replace the columns of the Expression_data file
colnames(expr) <- c(sample_id$ID)
sample_ids<-colnames(expr)
write.table(sample_ids,file=paste0(outdir,"sample_ids.txt"),col.names=F,row.names=F,quote=F,sep="\t")


# prepare covariates file
#############################
covariates <-t(v$target)
colnames(covariates) <- c(sample_id$ID)
covariates<- covariates[c('Age','snpPC1','snpPC2','snpPC3','snpPC4','snpPC5','snpPC6','snpPC7','snpPC8','snpPC9','snpPC10'),]
cov<-t(covariates)
write.table(cov,file=paste0(outdir,"covariates"),col.names=F,row.names=F,quote=F,sep="\t")

##get correlated Genes by performing cor function
##########################
cormat <- cor(t(expr), method = "pearson")
## Get only correlated genes >0.6
correlated_genes_name <- apply(cormat, 2, function(x) names(which(x>0.6)))

#Prepeare phenotype file
###########################

correlated_genes<-correlated_genes_name[[gene]]
idx <- match(correlated_genes,rownames(expr))
pheno <- expr[idx,]
pheno<-t(pheno)
write.table(pheno,paste0(outdir,"pheno","_",gene),col.names=T,row.names=T,quote=F,sep="\t")


## collect correlated genes snps
################################


## Combine gene annotation file and new expression file to get the chr,start,end position,gene name information for correlated genes
## Load gene annotation file
geneAnnotationFile <-"/dcs04/lieber/statsgen/shizhong/GRN/data/gene_annotation.tsv"
geneAnnotation <- read.table(geneAnnotationFile,header=T)
head(geneAnnotation)

## match IDs
idx <- match(correlated_genes,geneAnnotation$names)
RelatedGenes <- geneAnnotation[idx,c(1:5)]
## combine two files
exprRelatedGenes <- cbind(RelatedGenes,correlated_genes)

## create range file
library(dplyr)
exprRelatedGenes <-exprRelatedGenes  %>% select(seqnames, start,end,Symbol)
#remove chr character from the column
exprRelatedGenes$ seqnames<-gsub("chr","",as.character(exprRelatedGenes$seqnames))

write.table(exprRelatedGenes,paste0(outdir,"exprRelatedGenes","_",gene,".txt"),col.names=F,row.names=F,quote=F)

# collect correlated genes snps
#############################
plink_bfile <-"/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/plink_file/EA_all"

command<-paste("plink --bfile",plink_bfile," --write-snplist","--extract range",paste0(outdir,"exprRelatedGenes","_",gene,".txt"), "--out",paste0(outdir,"snps_correlated","_",gene))

system(command)


## Run .glm with these correlated genes snps to get cis snps
#############################################################


covar_file <- "/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/trans_qtl_discover/cov"

command<-paste("plink2 --bfile",plink_bfile," --write-snplist","--extract ",paste0(outdir,"snps_correlated","_",gene,".snplist"), "--out",paste0(outdir,"cis_snps","_",gene), "--covar",covar_file, "iid-only", "--pheno ",paste0(outdir,"pheno","_",gene),"iid-only", "--pfilter 0.05", "--glm hide-covar", "--covar-variance-standardize", "--mac 5")
system(command)



## make a new genotype file with those cis snp
########################################################

command<-paste("plink2 --bfile",plink_bfile," --make-bed","--extract ",paste0(outdir,"cis_snps","_",gene,".snplist"), "--out",paste0(outdir,"cis_genotype","_",gene))
system(command)


##run .glm again to get trans snps with new genotype
########################################################

command<-paste("plink2 --bfile",paste0(outdir,"cis_genotype","_",gene)," --write-snplist","--extract ",paste0(outdir,"cis_snps","_",gene,".snplist"), "--out",paste0(outdir,"trans_snps","_",gene), "--covar",covar_file, "iid-only", "--pheno ",paste0(outdir,"pheno","_",gene),"iid-only", "--pfilter 0.05", "--glm hide-covar", "--covar-variance-standardize", "--mac 5")
system(command)