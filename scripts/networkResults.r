
library(WGCNA)
library(ggplot2)
source("~/workspace/inheritance_estimate/scripts/samFunction.r")
#Data input
#Now it will be much more complicated for we have many regions. Most importantly, I want to separate the samples into different time periods
# Myelination 29weeks to end
# Neural migration, neurugenesis from week 4 to birth
# Synaptogenesis 20 week to 4 years of age
# Setting the reference files
# Because of sample size, let's do it prenatal and postnatal

brainSpan = "~/workspace/inheritance_estimate/network/expression_matrix.csv"
brainSpanRow="~/workspace/inheritance_estimate/network/rows_metadata.csv"
brainSpanCol="~/workspace/inheritance_estimate/network/columns_metadata.csv"
humanMetaFile = "~/workspace/inheritance_estimate/network/Human_Brain_Seq_Stages.csv"

# Read required files
colmeta = read.csv(brainSpanCol,row.names=1)
rowmeta=read.csv(brainSpanRow, row.names=1)
count = read.csv(brainSpan, row.names=1, header=F)
row.names(count) = rowmeta$ensembl_gene_id
humanMeta = read.csv(humanMetaFile)

networkFiles = list.files(pattern="network.csv")
significantModules = data.frame(region=character(), period=character(), module=character(), cor=numeric(), pvalue=numeric());

printNetworkGraph <- function(region, period, module, sampleInfo, selected, moduleGenes){
	# Get the count data
	library(ggplot2)
	countData = read.csv(paste(region, ".log.count", sep=""),row.names=1)
	maxX = max(plotData$Age);
	minX = min(plotData$Age);
	xCor = (maxX-minX)/10+minX
	maxY = max(plotData$expression);
	minY = min(plotData$expression);
	yCor = maxY-(maxY-minY)/10
	plotData=data.frame(Age = sampleInfo[selected,]$month, expression=colMeans((countData[row.names(countData)%in%moduleGenes,selected]+1)), age=sampleInfo[selected,]$age)
	#r2 = (cor(plotData$expression, predict(loess(plotData$expression~plotData$Age))))
	thesisPng(paste(region, "_",period,"_",module,".png",sep=""))
	gg=ggplot(plotData, aes(x=Age, y=expression))+geom_point()+labs(x = "Age (in months)", y = "Mean log2 RPKM")+theme_classic(base_size = 12)+ theme(axis.text.x = element_text(angle = 50, hjust = 1))+ theme(plot.title = element_text(lineheight=.8, face="bold",size=24),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+ stat_smooth(method = "loess", formula = y ~ x, size = 0.5,colour="black", alpha = 0.2)+scale_x_discrete(labels = plotData$Age, breaks=plotData$age)
	#+geom_text(x =xCor, y = yCor, label = paste("R^2: ", r2,sep=""),parse=TRUE)
	#+ ggtitle("Mean Expression of Genes in \nthe Co-expression network")
	print(gg)
	dev.off();
}

for(i in networkFiles){
	#For each network
	region = strsplit(i, split="\\.")[[1]][1]
	regionResult = i;
	regionCount = paste(region,".count",sep="")
	regionSample = paste(region, ".sampleInfo", sep="")
	modules = read.csv(regionResult, row.names=1)
	expr = read.csv(regionCount, row.names=1)
	datExpr = t(expr)
	names(datExpr) = row.names(expr)
	row.names(datExpr) = names(expr)
	sampleInfo = read.csv(regionSample, row.names=1)
	prenatalSamples =  row.names(datExpr)%in%as.character(sampleInfo$donor_name[sampleInfo$month < 12])
	postnatalSamples =  !(row.names(datExpr)%in%as.character(sampleInfo$donor_name[sampleInfo$month < 12]))
	PreMEs0 = moduleEigengenes(datExpr[prenatalSamples,], modules$moduleColor, excludeGrey=T)$eigengenes
	PreMEs = orderMEs(PreMEs0)
	PostMEs0 = moduleEigengenes(datExpr[postnatalSamples,], modules$moduleColor, excludeGrey=T)$eigengenes
	PostMEs = orderMEs(PostMEs0)
	AllMEs0 = moduleEigengenes(datExpr, modules$moduleColor, excludeGrey=T)$eigengenes
	AllMEs = orderMEs(AllMEs0)
	
	PremoduleTraitCor = cor(PreMEs, sampleInfo[sampleInfo$month < 12,]$month, use = "p");
	nSamples = sum(sampleInfo$month < 12)
	PremoduleTraitPvalue = corPvalueStudent(PremoduleTraitCor, nSamples);
	if(sum(PremoduleTraitPvalue[,1] < 0.5/nrow(PremoduleTraitPvalue))>0){
		sigModule = PremoduleTraitPvalue[,1] < 0.5/nrow(PremoduleTraitPvalue);
		for(sigNum in 1:sum(sigModule)){
			modName =gsub('ME','',row.names(AllmoduleTraitPvalue)[sigModule][sigNum]);
			significantModules = rbind(significantModules, data.frame(region=region, period = "Prenatal", module=modName, cor=(PremoduleTraitCor[sigModule,])[sigNum], pvalue=(PremoduleTraitPvalue[sigModule,])[sigNum]))
			#printNetworkGraph(region, "Prenatal", modName, sampleInfo, sampleInfo$month<12, as.character(subset(modules, modules$moduleColor==modName)$GeneID))
		}
	}
	
	PostmoduleTraitCor = cor(PostMEs, sampleInfo[!sampleInfo$month < 12,]$month, use = "p");
	nSamples = sum(!sampleInfo$month < 12)
	PostmoduleTraitPvalue = corPvalueStudent(PostmoduleTraitCor, nSamples);
	if(sum(PostmoduleTraitPvalue[,1] < 0.5/nrow(PostmoduleTraitPvalue))>0){
		sigModule = PostmoduleTraitPvalue[,1] < 0.5/nrow(PostmoduleTraitPvalue);
		for(sigNum in 1:sum(sigModule)){
			modName =gsub('ME','',row.names(AllmoduleTraitPvalue)[sigModule][sigNum]);
			significantModules = rbind(significantModules, data.frame(region=region, period = "Postnatal", module=modName, cor=(PostmoduleTraitCor[sigModule,])[sigNum], pvalue=(PostmoduleTraitPvalue[sigModule,])[sigNum]))
			#printNetworkGraph(region, "Postnatal", modName, sampleInfo, !sampleInfo$month<12, as.character(subset(modules, modules$moduleColor==modName)$GeneID))
		}
	}
	AllmoduleTraitCor = cor(AllMEs, sampleInfo$month, use = "p");
	nSamples = length(sampleInfo$month)
	AllmoduleTraitPvalue = corPvalueStudent(AllmoduleTraitCor, nSamples);
	if(sum(AllmoduleTraitPvalue[,1] < 0.5/nrow(AllmoduleTraitPvalue))>0){
		sigModule = AllmoduleTraitPvalue[,1] < 0.5/nrow(AllmoduleTraitPvalue);
		for(sigNum in 1:sum(sigModule)){
			modName =gsub('ME','',row.names(AllmoduleTraitPvalue)[sigModule][sigNum]);
			significantModules = rbind(significantModules, data.frame(region=region, period = "All", module=modName, cor=(AllmoduleTraitCor[sigModule,])[sigNum], pvalue=(AllmoduleTraitPvalue[sigModule,])[sigNum]))
			#printNetworkGraph(region, "All", modName, sampleInfo, rep(TRUE, nrow(sampleInfo)), as.character(subset(modules, modules$moduleColor==modName)$GeneID))
		}
	}
}