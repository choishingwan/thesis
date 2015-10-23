
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
significantModules = data.frame(region=character(), period=character(), module=character(), cor=numeric(), pvalue=numeric(),lmR2=numeric(), lmP=numeric(), modSize=numeric());
TotalModules=0;

#Get the month for each age
colmeta$month = 0;
for(i in 1:nrow(colmeta)){
	temp = as.character(colmeta$age)[i]
	temp = strsplit(temp, split=" ")[[1]]
	if(temp[2] == "pcw"){
		colmeta$month[i] = as.numeric(as.character(temp[1]))/4
	}else if(temp[2]=="mos"){
		colmeta$month[i] = as.numeric(as.character(temp[1]))+10
	}else if (temp[2]=="yrs"){
		colmeta$month[i] = as.numeric(as.character(temp[1]))*12+10
	}
}

printNetworkGraph <- function(region, period, module, colmeta, moduleGenes){
	print(paste(region, period, module))
	
	# Get the count data
	library(ggplot2)
	countData = read.csv(paste(region, ".log.count", sep=""),row.names=1)
	samples = subset(colmeta, colmeta$structure_acronym==region)
	colnames(countData) = samples$donor_name
	if(period =="Prenatal"){
		selected = samples$month <= 10;
	}else if(period=="Postnatal"){
		selected = samples$month > 10;
	}else{
		selected = rep(TRUE, nrow(samples))
	}
	
	plotData=data.frame(Age = samples[selected,]$month, expression=colMeans((countData[row.names(countData)%in%moduleGenes,selected]+1)), age=samples[selected,]$age)
	maxX = max(plotData$Age);
	minX = min(plotData$Age);
	xCor = (maxX-minX)/10+minX
	maxY = max(plotData$expression);
	minY = min(plotData$expression);
	yCor = maxY-(maxY-minY)/10
	r2 = (cor(plotData$expression, predict(loess(plotData$expression~plotData$Age))))
	thesisPng(paste(region, "_",period,"_",module,".png",sep=""))
	gg=ggplot(plotData, aes(x=Age, y=expression))+geom_point()+labs(x = "Age (in months)", y = "Mean log2 RPKM")+theme_classic(base_size = 12)+ theme(axis.text.x = element_text(angle = 50, hjust = 1))+ theme(plot.title = element_text(lineheight=.8, face="bold",size=24),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+ stat_smooth(method = "loess", formula = y ~ x, size = 0.5,colour="black", alpha = 0.2)+geom_text(x =xCor, y = yCor, label = paste("R^2: ", r2,sep=""),parse=TRUE)
	#+scale_x_discrete(labels = plotData$Age, breaks=plotData$age)
	#+ ggtitle("Mean Expression of Genes in \nthe Co-expression network")
	print(gg)
	dev.off();
}

for(i in networkFiles){
	#For each network
	region = strsplit(i, split="\\.")[[1]][1]
	regionResult = i;
	regionCount = paste(region,".log.count",sep="")
	regionSample = paste(region, ".sampleInfo", sep="")
	modules = read.csv(regionResult, row.names=1)
	TotalModules = TotalModules+length(unique(modules$moduleColor))-1
	expr = read.csv(regionCount, row.names=1)
	sampleInfo = read.csv(regionSample, row.names=1)
	
	datExpr = t(expr)
	names(datExpr) = row.names(expr)
	row.names(datExpr) = names(expr)
	datExpr = datExpr[as.character(row.names(datExpr))%in%as.character(sampleInfo$donor_name),]
	prenatalSamples =  row.names(datExpr)%in%as.character(sampleInfo$donor_name[sampleInfo$month <= 10])
	postnatalSamples =  !(row.names(datExpr)%in%as.character(sampleInfo$donor_name[sampleInfo$month <= 10]))
	PreMEs0 = moduleEigengenes(datExpr[prenatalSamples,], modules$moduleColor, excludeGrey=T)$eigengenes
	PreMEs = orderMEs(PreMEs0)
	PostMEs0 = moduleEigengenes(datExpr[postnatalSamples,], modules$moduleColor, excludeGrey=T)$eigengenes
	PostMEs = orderMEs(PostMEs0)
	AllMEs0 = moduleEigengenes(datExpr, modules$moduleColor, excludeGrey=T)$eigengenes
	AllMEs = orderMEs(AllMEs0)
	
	PremoduleTraitCor = cor(PreMEs, sampleInfo[sampleInfo$month <= 10,]$month, use = "p");
	nSamples = sum(sampleInfo$month <= 10)
	PremoduleTraitPvalue = corPvalueStudent(PremoduleTraitCor, nSamples);
	if(sum(PremoduleTraitPvalue[,1] < 0.5/nrow(PremoduleTraitPvalue))>0){
		sigModule = PremoduleTraitPvalue[,1] < 0.5/nrow(PremoduleTraitPvalue);
		for(sigNum in 1:sum(sigModule)){
			modName =gsub('ME','',row.names(PremoduleTraitPvalue)[sigModule][sigNum]);
			#significantModules = rbind(significantModules, data.frame(region=region, period = "Prenatal", module=modName, cor=(PremoduleTraitCor[sigModule,])[sigNum], pvalue=(PremoduleTraitPvalue[sigModule,])[sigNum],modSize=length(as.character(subset(modules, modules$moduleColor==modName)$GeneID))))
			#printNetworkGraph(region, "Prenatal", modName, colmeta, as.character(subset(modules, modules$moduleColor==modName)$GeneID))
		}
	}
	
	PostmoduleTraitCor = cor(PostMEs, sampleInfo[!sampleInfo$month <= 10,]$month, use = "p");
	nSamples = sum(!sampleInfo$month <= 10)
	PostmoduleTraitPvalue = corPvalueStudent(PostmoduleTraitCor, nSamples);
	if(sum(PostmoduleTraitPvalue[,1] < 0.5/nrow(PostmoduleTraitPvalue))>0){
		sigModule = PostmoduleTraitPvalue[,1] < 0.5/nrow(PostmoduleTraitPvalue);
		for(sigNum in 1:sum(sigModule)){
			modName =gsub('ME','',row.names(PostmoduleTraitPvalue)[sigModule][sigNum]);
			#significantModules = rbind(significantModules, data.frame(region=region, period = "Postnatal", module=modName, cor=(PostmoduleTraitCor[sigModule,])[sigNum], pvalue=(PostmoduleTraitPvalue[sigModule,])[sigNum],modSize=length(as.character(subset(modules, modules$moduleColor==modName)$GeneID))))
			#printNetworkGraph(region, "Postnatal", modName, colmeta, as.character(subset(modules, modules$moduleColor==modName)$GeneID))
		}
	}
	AllmoduleTraitCor = cor(AllMEs, sampleInfo$month, use = "p");
	nSamples = length(sampleInfo$month)
	AllmoduleTraitPvalue = corPvalueStudent(AllmoduleTraitCor, nSamples);
	if(sum(AllmoduleTraitPvalue[,1] < 0.5/nrow(AllmoduleTraitPvalue))>0){
		sigModule = AllmoduleTraitPvalue[,1] < 0.5/nrow(AllmoduleTraitPvalue);
		for(sigNum in 1:sum(sigModule)){
			modName =gsub('ME','',row.names(AllmoduleTraitPvalue)[sigModule][sigNum]);
			significantModules = rbind(significantModules, data.frame(region=region, period = "All", module=modName, cor=(AllmoduleTraitCor[sigModule,])[sigNum], pvalue=(AllmoduleTraitPvalue[sigModule,])[sigNum],modSize=length(as.character(subset(modules, modules$moduleColor==modName)$GeneID))))
			#printNetworkGraph(region, "All", modName, colmeta, as.character(subset(modules, modules$moduleColor==modName)$GeneID))
		}
	}
}

significantModules$fdrS = p.adjust(significantModules$pvalue, method="BH", n=TotalModules*3)#Because prenatal + postnatal + all
trueSig = subset(significantModules, significantModules$fdrS<0.05)
for(i in 1:nrow(trueSig)){
	mod = read.csv(paste(trueSig$region[i], ".network.csv", sep=""), row.names=1)
	modName = as.character(trueSig$module[i])
	region = as.character(trueSig$region[i])
	printNetworkGraph(region, trueSig$period[i], modName, colmeta, as.character(subset(mod, mod$moduleColor==modName)$GeneID))
	geneset = rowmeta[rowmeta$ensembl_gene_id%in%as.character(subset(mod, mod$moduleColor==modName)$GeneID),]$entrez_id	
	geneset = unique(geneset)
	geneset = geneset[!is.na(geneset)]
	write.table(geneset, paste(trueSig$region[i], "_",modName, "_", trueSig$period[i], ".temp",sep=""),row.names=F, col.names=F, quote=F);	
}


ls *.temp | awk -F "." '{print "awk \x27{printf $0\"\\t\"}END{print \"\"}\x27 "$0"> "$1".tempSet"}' | bash
ls *.temp | awk -F "." '{print "echo -n \""$1" \"| cat - "$1".tempSet > "$1".set"}'   | bash

rm finalSet
for i in `ls *.set`; do cat $i >> finalSet; done
rm *.tempSet

#Old Scripts

#moduleResult = "Amy.network.csv"
moduleResult = "hippocampus/Hippo.longshot.network.csv"
#networkinput = "Amy.normalize.csv"
networkinput = "hippocampus/Hippo.normalized.count.csv"
#hipMetaFile = "amyMeta.csv"
hipMetaFile = "hippocampus/hipMeta.csv"
#acronym = "AMY"
acronym = "HIP"
hipAge = c(1,2,3,3,3,4,4,4,5,5,6,7,8,9,10,11,11,12,13,14,15,15,16,17,18,19,20,21,22,23,24,25)
amyAge = c(1,1,2,2,2,3,3,3,4,4,5,6,7,8,9,9,9,10,11,11,12,13,13,14,15,16,17,18,19,20,21,22,23)
modules = read.csv(moduleResult, row.names=1)
expr = read.csv(networkinput, row.names=1)
#Calculate the eigen-genes
datExpr = t(expr)
names(datExpr) = row.names(expr)
row.names(datExpr) = names(expr)
#The softpower won't be used as we don't need to use the hub gene. But we still enter the softpower we used for building the network
MEs0 = moduleEigengenes(datExpr, modules$moduleColor, softPower=28,excludeGrey=T)$eigengenes
MEs = orderMEs(MEs0)
#Perform correlation between the age and the module eigen-genes
traits = read.csv(hipMetaFile, row.names=1)
#The age is coded as month starts from conception
modTraits = data.frame(row.names=traits$donor_name, gd=traits$Month, gender = traits$gender)
moduleTraitCor = cor(MEs, modTraits$gd, use = "p");
nSamples = nrow(modTraits)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sigModuleColors = gsub('ME','',row.names(subset(moduleTraitPvalue, moduleTraitPvalue[,1] < 0.05/nrow(moduleTraitCor))))
count.hippo=count[,colmeta$structure_acronym==acronym]
hipSample = subset(colmeta, colmeta$structure_acronym==acronym)
colnames(count.hippo) = hipSample$donor_name
#Now for each significant modules, we draw the trend
for(i in sigModuleColors){
count.hippo.select= count.hippo[row.names(count.hippo)%in%row.names(modules)[modules$moduleColor%in%i],]
plotData=data.frame(age = factor(hipSample$age,levels=c(unique(as.character(hipSample$age)))), expression=colMeans(log2(count.hippo.select+1)))
plotData$Age = 0
plotData$Age=hipAge
#plotData$Age=amyAge
png(paste("Network_Expression_",i,".png",sep=""), width = 8, height = 8, units = 'in', res = 300)
print(cor(plotData$expression, predict(loess(plotData$expression~plotData$Age))))
gg=ggplot(plotData, aes(x=Age, y=expression))+geom_point()+labs(x = "Age", y = "Mean log2 PRKM")+theme_classic(base_size = 12)+ theme(axis.text.x = element_text(angle = 50, hjust = 1))+ theme(plot.title = element_text(lineheight=.8, face="bold",size=24),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+ stat_smooth(method = "loess", formula = y ~ x, size = 0.5,colour="black", alpha = 0.2)+scale_x_discrete(labels = plotData$age, breaks=plotData$Age)#+ ggtitle("Mean Expression of Genes in \nthe Co-expression network")
print(gg)
dev.off()
}

for(i in 1:nrow(hipMeta)){
	temp = as.character(hipMeta$age)[i]
	temp = strsplit(temp, split=" ")[[1]]
	if(temp[2] == "pcw"){
		hipMeta$Month[i] = as.numeric(as.character(temp[1]))/4
	}else if(temp[2]=="mos"){
		hipMeta$Month[i] = as.numeric(as.character(temp[1]))+10
	}else if (temp[2]=="yrs"){
		hipMeta$Month[i] = as.numeric(as.character(temp[1]))*12+10
	}
}
