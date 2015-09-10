
library(WGCNA)
library(ggplot2)
#Data input

moduleResult = "Amy.network.csv"
#moduleResult = "hippocampus/Hippo.longshot.network.csv"
networkinput = "Amy.normalize.csv"
#networkinput = "hippocampus/Hippo.normalized.count.csv"
brainSpan = "expression_matrix.csv"
brainSpanRow="rows_metadata.csv"
brainSpanCol="columns_metadata.csv"
hipMetaFile = "amyMeta.csv"
#hipMetaFile = "hippocampus/hipMeta.csv"
humanMetaFile = "Human_Brain_Seq_Stages.csv"
acronym = "AMY"
#acronym = "HIP"
hipAge = c(1,2,3,3,3,4,4,4,5,5,6,7,8,9,10,11,11,12,13,14,15,15,16,17,18,19,20,21,22,23,24,25)
amyAge = c(1,1,2,2,2,3,3,3,4,4,5,6,7,8,9,9,9,10,11,11,12,13,13,14,15,16,17,18,19,20,21,22,23)
modules = read.csv(moduleResult, row.names=1)
expr = read.csv(networkinput, row.names=1)
colmeta = read.csv(brainSpanCol,row.names=1)
rowmeta=read.csv(brainSpanRow, row.names=1)
count = read.csv(brainSpan, row.names=1, header=F)
row.names(count) = rowmeta$ensembl_gene_id
humanMeta = read.csv(humanMetaFile)
#Calculate the eigen-genes
datExpr = t(expr)
names(datExpr) = row.names(expr)
row.names(datExpr) = names(expr)
#The softpower won't be used as we don't need to use the hub gene. But we still enter the softpower we used for building the network
MEs0 = moduleEigengenes(datExpr, modules$moduleColor, softPower=15,excludeGrey=T)$eigengenes
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
#plotData$Age=hipAge
plotData$Age=amyAge
png(paste("Network_Expression_",i,".png",sep=""), width = 8, height = 8, units = 'in', res = 300)
gg=ggplot(plotData, aes(x=Age, y=expression))+geom_point()+labs(x = "Age", y = "Mean log2 PRKM")+theme_classic(base_size = 12)+ theme(axis.text.x = element_text(angle = 50, hjust = 1))+ theme(plot.title = element_text(lineheight=.8, face="bold",size=24),axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))+ stat_smooth(method = "loess", formula = y ~ x, size = 0.5,colour="black", alpha = 0.2)+scale_x_discrete(labels = plotData$age, breaks=plotData$Age)#+ ggtitle("Mean Expression of Genes in \nthe Co-expression network")
print(gg)
dev.off()
}

