sampleSelect <-function(rnaMeta, sampleMeta, region){
	#Change the age label such that they match eachother
	temp = merge(rnaMeta, sampleMeta, by.x=c("donor_name","structure_acronym"),by.y=c("Allen.Institute.ID","Region.Area"))
	temp = subset(temp, temp$structure_acronym==region & temp$Dissection.score>3& temp$RIN >=7.0)
	if(nrow(temp)==0){
		return(NULL)
	}else{
		temp = temp[order(temp$Age, temp$Dissection.score, temp$RIN),]
		select = rep(FALSE, nrow(temp));
		prev = NULL;
		temp$month = 0;
		for(i in 1:nrow(temp)){
			if(is.null(prev) || as.character(prev)!=as.character(temp[i,]$age)){
				select[i] = TRUE;
				prev = as.character(temp[i,]$age)
				Age = strsplit(prev, split=" ")[[1]]
				if(Age[2] == "pcw"){
					temp$month[i] = as.numeric(as.character(Age[1]))/4
				}else if(Age[2]=="mos"){
					temp$month[i] = as.numeric(as.character(Age[1]))+10
				}else if (Age[2]=="yrs"){
					temp$month[i] = as.numeric(as.character(Age[1]))*12+10
				}
			}
		}
		#We only want sample of age less than or equal to 30
		result = temp[select,]
		result = subset(result, result$month <= 370)
		return(result[order(result$month),])
	}
}
#Data info
brainSpan = "~/workspace/inheritance_estimate/scripts/network/expression_matrix.csv"
brainSpanRow="~/workspace/inheritance_estimate/scripts/network/rows_metadata.csv"
brainSpanCol="~/workspace/inheritance_estimate/scripts/network/columns_metadata.csv"
humanMetaFile = "~/workspace/inheritance_estimate/scripts/network/Human_Brain_Seq_Stages.csv"

# Reading files
colmeta = read.csv(brainSpanCol,row.names=1)
rowmeta=read.csv(brainSpanRow, row.names=1)
count = read.csv(brainSpan, row.names=1, header=F)
row.names(count) = rowmeta$ensembl_gene_id
humanMeta = read.csv(humanMetaFile)

# Start process
# As long as the brain region have more than 10 samples each with different age, we will try to construct the network
# Because Schizophrenia seems to be related to almost any brain regions...
# Most important regions are STR, HIP, AMY and DFC

okRegions= NULL
for(i in unique(colmeta$structure_acronym)){
	temp = sampleSelect(colmeta, humanMeta, i)
	if(!is.null(temp)){
		if(dim(temp)[1] > 10){
			okRegions=rbind(okRegions, i)
			select = colmeta$donor_name %in% temp$donor_name & colmeta$structure_acronym %in% temp$structure_acronym
			select.count = count[,select]
			colnames(select.count) = colmeta$donor_name[select]
			select.count.filter = select.count[rowMeans(select.count)>1,]
			select.count.log2 = log2(select.count.filter+1)
			select.count.sd = apply(select.count.log2, 1, sd)
			select.count.norm = (select.count.log2-rowMeans(select.count.log2))/select.count.sd
			write.csv(select.count.norm, paste(i, ".count", sep=""))
			write.csv(select.count.log2, paste(i, ".log.count", sep=""))
			write.csv(temp,paste(i,".sampleInfo", sep=""))
			#system(paste("asub -c \"Rscript ~/workspace/inheritance_estimate/scripts/network/wgcna.r -f ",i,".count -o ",i,".network -s 0 -t 10 \" -e choishingwan@gmail.com -m 40gb -n nodes=1:ppn=12 -j ",i,".network",sep=""));
			# Due to server restrictions and some random reasons, the above command will only generate all the required cmd files but not submitting the job to server.
		}
	}
}

