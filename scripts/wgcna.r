#!/usr/bin/env Rscript
#Check if the optparse library is installed, required to parse arguments
`%notin%` <- function(x,y) !(x %in% y) 

if(suppressPackageStartupMessages(require("optparse"))){
    print("optparse is loaded correctly");
} else {
    print ("trying to install optparse")
    install.packages("optparse",repos="http://cran.rstudio.com/");
    if(suppressPackageStartupMessages(require("optparse"))){
	print ("optparse installed and loaded");
    } else{
	stop("could not install optparse");
    }
}



options(stringsAsFactors = FALSE);

#Perform the argument checking before spending time in loading DESeq2
option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false", dest="verbose", help="Print little output"),
    make_option(c("-f", "--file"), type="character",  help="The input count table", dest="inputFile"),
    make_option(c("-t", "--thread"), type="integer", default=10, help="Number of thread to use [default \"%default\"]", metavar="number"),
    make_option(c("-s", "--softpower"), type="integer", default=0, help="Selected softpower. if 0, will automatically calculate the softpower and select one at the first peak whereas an input of -1 will indicates the use of estimated power from WGCNA [default \"%default\"]", metavar="number"),
    make_option(c("-o", "--out"), type="character", help="Output prefix [default \"count table name\"]", metavar="output")
)

 
opt <- parse_args(OptionParser(option_list=option_list))

if(is.null(opt$inputFile)){
    stop("Require count table input");
}
if(is.null(opt$out)){
    opt$out = opt$inputFile;
}
if(file.access(opt$inputFile) ==-1){
    stop(sprintf("Count table (%s) does not exist", opt$inputFile));
} else{
    sprintf("reading from %s", opt$inputFile);
    inputTable = read.csv(opt$inputFile, row.names=1, header=T);
}
print("Count table read");

if(suppressPackageStartupMessages(require("WGCNA"))){
    print("WGCNA is loaded correctly");
} else {
    print ("trying to install WGCNA")
    source("http://bioconductor.org/biocLite.R")
    biocLite();
    biocLite("WGCNA");
    if(suppressPackageStartupMessages(require("WGCNA"))){
	print ("WGCNA installed and loaded");
    } else{
	stop("could not install WGCNA");
    }
}
allowWGCNAThreads(opt$thread);

datExpr0 = as.data.frame(t(inputTable));
names(datExpr0) = rownames(inputTable);
rownames(datExpr0) = names(inputTable);

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK){
# Optionally, print the gene and sample names that were removed:

if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


sampleTrees = hclust(dist(datExpr), method = "average")

#Plot the cluster
pdf(file = paste(opt$out,"sampleClustering.pdf", sep="_"), width = 12, height = 12);
    par(mfrow=c(2,1))
    par(mar = c(0, 4, 2, 0))
    plot(sampleTrees, main = "Sample clustering on all genes",xlab="", sub="", cex = 0.7);
dev.off();

powers = c(c(1:10), seq(from = 12, to=30, by=1))
#If it crash, it is most likely in this part. From my experience, it is usually because it has used too much cpu
sft=NULL;
if(opt$softpower <=0){
    sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
}
prev = NULL;
prevPower=NULL
#should allow users to select their power. three options, user number, estimated power, and our slope calculation model
selectedPower = NULL;
if(opt$softpower == -1){
    selectedPower=sft$powerEstimate
} else if(opt$softpower == 0){
    for(i in 2:length(sft$fitIndices$SFT.R.sq)){
	if(sft$fitIndices$slope[i] <0){
	    if(is.null(prev)){
		prev = (sft$fitIndices$SFT.R.sq[i]-sft$fitIndices$SFT.R.sq[i-1])/(sft$fitIndices$Power[i]-sft$fitIndices$Power[i-1])
		prevPower =sft$fitIndices$Power[i-1];
	    } else{
		newSlop = (sft$fitIndices$SFT.R.sq[i]-sft$fitIndices$SFT.R.sq[i-1])/(sft$fitIndices$Power[i]-sft$fitIndices$Power[i-1]);
		if(newSlop <=0){
		    selectedPower=prevPower
		    break;
		}
		prev = newSlop;
		prevPower = sft$fitIndices$Power[i]
	    }
	}
    }
} else{
    selectedPower = opt$softpower;
}

sprintf("Selected soft power = %i", selectedPower);
mod= blockwiseModules(datExpr, corType="bicor", power=selectedPower, nThreads=opt$thread, maxBlockSize=nGenes, quickCor=0, impute=TRUE, minModuleSize=30,saveTOMs =TRUE)
ADJ1 = abs(bicor(datExpr))^selectedPower
TOM = TOMsimilarity(ADJ1);
Alldegrees1=intramodularConnectivity(ADJ1, mod$colors)
MEs=mod$MEs
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));

geneInfo = data.frame(GeneID=names(datExpr), moduleColor=mod$colors, geneModuleMembership ,Alldegrees1)
write.csv(geneInfo, paste(opt$out,".csv", sep=""));

