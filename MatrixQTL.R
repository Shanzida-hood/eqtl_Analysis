#!/usr/bin/env Rscript --vanilla
args = commandArgs(trailingOnly=TRUE)
gene = args[1]

## create outdir
##########################
outdir<-"/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/matrixQTL/"

setwd("/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/matrixQTL/")

test<-"/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/matrixQTL/test/"


##############################################################################

## create expression file
###########################

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
#create genotype file
########################
plink_bfile <-"/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/plink_file/EA_all"
command <- paste("plink2 --bfile",plink_bfile, "--make-bed", "--keep",paste0(outdir,"sample_ids.txt"), "--out",paste0(outdir,"sample"))
system(command)
## create final sample ids
#############################
sample_ID_fam<-read.table("/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/matrixQTL/sample.fam", header = F)
colnames(sample_ID_fam)
sample_id_fam<-(sample_ID_fam$V2)
final_sample_id<-intersect(sample_ids,sample_id_fam)
final_sample_id
write.table(final_sample_id,paste0(outdir,"final_sample_id.txt"),col.names=T,row.names=T,quote=F,sep="\t")

## write expression file 
expression<-expr[,final_sample_id]
write.table(expression,paste0(outdir,"expression.txt"),col.names=T,row.names=T,quote=F,sep="\t")

## create gene_location file
geneAnnotationFile <-"/dcs04/lieber/statsgen/shizhong/GRN/data/gene_annotation.tsv"
geneAnnotation <- read.table(geneAnnotationFile,header=T)
## select columns gene name,chr,end,start
gene_location <- geneAnnotation[,c(1:4)]
library(dplyr)
## change columns name
gene_location <- gene_location %>% 
  rename("geneid" = "names",
         "chr" = "seqnames",
         "s1" = "start","s2"="end")
write.table(gene_location,"gene_location.txt", sep = "\t",
            row.names=F,col.names=T,quote=F)


## create snp_location file
# read bim file
SNP_location<-read.table("/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/matrixQTL/sample.bim", header = F)
#select columns
SNP_location <- SNP_location[,c(1:4)]
#reorder columns position
SNP_location <- SNP_location[, c(2, 1, 4)]
#add chr string
SNP_location$V1<-sub("^","chr",SNP_location$V1)
SNP_location <- SNP_location %>% 
  rename("snp" = "V2",
         "chr" = "V1",
         "pos" = "V4")
rownames(SNP_location) <-SNP_location$snp
write.table(SNP_location,paste0(outdir,"SNP_location.txt"), sep = "\t",
            row.names=F,col.names=T,quote=F)


## create covariates file

covariates <-t(v$target)
colnames(covariates) <- c(sample_id$ID)
covariates<- covariates[c('Age','snpPC1','snpPC2','snpPC3','snpPC4','snpPC5','snpPC6','snpPC7','snpPC8','snpPC9','snpPC10'),]
covariates<-covariates[,final_sample_id]
write.table(covariates,paste0(outdir,"covariates.txt"),col.names=T,row.names=T,quote=F,sep="\t")

##Create genotype file
##########################
library('BEDMatrix',lib='/users/sjahansi/R/4.2.x')
m <- BEDMatrix("/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/matrixQTL/sample")
dim(m)
m[1:3, 1:5]
library('genio',lib='/users/sjahansi/R/4.2.x')
time_read_bedmatrix_2 <- system.time(
  X_BEDMatrix_Rmat <- as.matrix(m)
)
genotype<-t(X_BEDMatrix_Rmat)


## Remove special string from row names
rownames(genotype)<-gsub("_.*","",rownames(genotype))

colnames(genotype)<-c(final_sample_id)
write.table(genotype,paste0(outdir,"genotype.txt"),col.names=T,row.names=T,quote=F,sep="\t")


##get correlated Genes by performing cor function
##########################
cormat <- cor(t(expression), method = "pearson")
## Get only correlated genes >0.6
correlated_genes_name <- apply(cormat, 2, function(x) names(which(x>0.6)))

#Prepeare expression file
###########################

correlated_genes<-correlated_genes_name[[gene]]
idx <- match(correlated_genes,rownames(expression))
gene_expression <- expression[idx,]
write.table(gene_expression,paste0(outdir,"correlated_gene_expression","_",gene,".txt"),col.names=T,row.names=T,quote=F,sep = "\t")
## create gene location file 
################################
## Combine gene annotation file and new expression file to get the chr,start,end position,gene name information for correlated genes
## Load gene annotation file
geneAnnotationFile <-"/dcs04/lieber/statsgen/shizhong/GRN/data/gene_annotation.tsv"
geneAnnotation <- read.table(geneAnnotationFile,header=T)
head(geneAnnotation)

## match IDs
idx <- match(correlated_genes,geneAnnotation$names)
RelatedGenes <- geneAnnotation[idx,c(1:4)]
library(dplyr)
## change columns names of gene location file 
gene_location <-RelatedGenes  %>% 
  rename("geneid" = "names",
         "chr" = "seqnames",
         "s1" = "start","s2"="end")
rownames(gene_location) <- gene_location$geneid
write.table(gene_location,paste0(outdir,"gene_location","_",gene,".txt"),col.names=T,row.names=F,quote=F,sep = "\t")
gene_expression_2<-read.table(paste0(outdir,"correlated_gene_expression","_",gene,".txt"))
colnames(gene_expression_2)<-gsub("X","",colnames(gene_expression_2))
# collect correlated genes snps
#############################
for(i in 1:nrow(gene_location)){
  id<-  gene_location$ geneid[i]
  
  ## Run matrixqtl with these correlated genes snps to get cis snps
  #############################################################
  ## create gene location file   
  gene_location_2<-as.data.frame(gene_location[id,], drop=false)  
  write.table(gene_location_2,paste0(test,"gene_location_2","_",id),row.names = F,col.names = T,quote = F,sep="\t")
  
  ## create gene expression file
  gene_expression_3<-gene_expression_2[id,]
  
  write.table(gene_expression_2,paste0(test,"gene_expression_2","_",id),row.names = T,col.names = T,quote = F,sep="\t")
  
  library('MatrixEQTL',lib='/users/sjahansi/R/4.2.x')
  
  ## Settings
  
  # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  
  
  # Genotype file name
  SNP_file_name = paste0(outdir,"genotype.txt");
  snps_location_file_name = paste0(outdir,"SNP_location.txt");
  
  # Gene expression file name
  expression_file_name = paste0(test,"gene_expression_2","_",id);
  gene_location_file_name = paste0(test,"gene_location_2","_",id);
  
  # Covariates file name
  # Set to character() for no covariates
  covariates_file_name = paste0(outdir,"covariates.txt");
  
  # Output file name
  output_file_name_cis <-paste0(test,"cis","_",id);
  output_file_name_tra <-paste0(test,"trans","_",id);
  
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis = 1e-2;
  pvOutputThreshold_tra = 1e-2;
  
  # Error covariance matrix
  # Set to numeric() for identity.
  errorCovariance = numeric();
  # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
  
  # Distance for local gene-SNP pairs
  cisDist = 1e6;
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = "\t";      # the TAB character
  cvrt$fileOmitCharacters = "NA"; # denote missing values;
  cvrt$fileSkipRows = 1;          # one row of column labels
  cvrt$fileSkipColumns = 1;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  
  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  
  
  ###################################################  
  ##run matrixqtl again to get trans snps with new genotype
  ########################################################  
  ## step_1:make a new genotype file with those cis snp
  ########################################################
  cis_df<- read.table(output_file_name_cis, header = TRUE, stringsAsFactors = FALSE); 
  new_snp_id<-cis_df$SNP
  new_genotype<-genotype[new_snp_id,]
  write.table(new_genotype,paste0(test,"_","new_genotype","_",id,".txt"),col.names=T,row.names=T,quote=F,sep="\t")
  ## step2:create new snps location file with cis snps
  ################################################
  new_snp_location<-SNP_location[new_snp_id,]
  write.table(new_snp_location,paste0(test,"_","new_snp_location","_",id,".txt"),sep ="\t",
              row.names = F, col.names = T,quote=F)
  # Genotype file name
  SNP_file_name = paste0(test,"_","new_genotype","_",id,".txt");
  snps_location_file_name = paste0(test,"_","new_snp_location","_",id,".txt");
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  # Output file name
  output_file_name_cis <-paste0(outdir,"cis_new_genotype",id);
  output_file_name_tra <-paste0(outdir,"trans_new_genotype",id);
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis = 4e-2;
  pvOutputThreshold_tra = 1e-2;
  
  # Error covariance matrix
  # Set to numeric() for identity.
  errorCovariance = numeric();
  # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
  
  # Distance for local gene-SNP pairs
  cisDist = 1e6;
  
  ## Run the analysis
  snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
}

#remove temporary directory
######################
unlink(output_file_name_tra)
unlink("/dcs04/lieber/statsgen/sjahansi/Project_Drug_Discovery/Data/matrixQTL/test/"
       ,recursive=TRUE)











[compute-124 /dcs04/lieber/statsgen/sjahansi]$ 
  