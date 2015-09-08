args = commandArgs(trailingOnly = TRUE)
input=args[1];
error = FALSE;
if(is.na(input)){
    error = TRUE;
    print("Please provide the input file");
}
output=args[2];
if(is.na(output)){
    error = TRUE;
    print("Please provide the output name");
}
if(error) stop();

data =read.csv(input, row.names=1);
#Remove low expression
data.select = data[rowMeans(data)>1,]
#Log transform the data
data.log =log2(data.select+1)
#Normalized the data
data.sd = apply(data.log, 1, sd)
data.norm = (data.log-rowMeans(data.log))/data.sd
write.csv(data.norm, output);

