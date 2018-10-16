# Author : DE OLIVEIRA, GONCALVES CLARO 
# Date : 10/10/2018
# Miaou Project 2018

# Read table
data <- read.table("../donnees/P_EXPR.txt", header = T)

# Get GOI (gene of interest) data
GOI <- read.table("../donnees/P_GOI.txt", header=T)
list_GOI <- as.vector(GOI[,"id"])

# Filter the data
filter <- data[,"id"] %in% list_GOI
data_filtered <- data[filter,]

# Get id as rownames and del first col
rownames(data_filtered) <- data_filtered[,"id"]
data_filtered <- data_filtered[,-1]

# Transpose rownames and colnames (get gene as attributes and measurement as element)
data_filtered <- t(data_filtered)

# Process correlation between GOI 
result_correlation <- cor(data_filtered)

#Use library reshape to convert
library("reshape2")
result_correlation=melt(result_correlation)
result_correlation

# Get only values > 0.8
result_correlation <- result_correlation[abs(result_correlation["value"])>0.8,]


# Sort duplicate, redundant data
result_correlation["alphabetic"]<-as.character(result_correlation[,"Var1"])<as.character(result_correlation[,"Var2"])
result_correlation=result_correlation[result_correlation[,4]==TRUE,]

# Keep sens of correlation
result_dir_correlation <- result_correlation[["value"]] > 0
result_dir_correlation

# Delete columns Value and alphabetic
result_correlation <- result_correlation[,c(-4,-3)]
result_correlation["type"] <- "correlation_pair"

# Download PPI 
data_PPI <- read.table("../donnees/P_PPI.txt",sep = "\t", header = T)

# Change column name to get identical colnames PPI and correlation table
colnames(result_correlation)[1]="id1"
colnames(result_correlation)[2]="id2"
colnames(data_PPI)[3]="type"

# Filter gene name PPI with gene name correlation
id1<-as.vector(result_correlation[,"id1"])
id2<-as.vector(result_correlation[,"id2"])
result_net_genes <- c(id1,id2)
result_net_genes <- unique(result_net_genes)
for (i in c("id1","id2")) {
  data_PPI_filter<- data_PPI[,i] %in% result_net_genes
  data_PPI <- data_PPI[data_PPI_filter,]
}

# Bind the two table with extra column (direction of correlation)
LISTE <- rbind(result_correlation, data_PPI)
LISTE[["cor_dir"]] <- "not_relevant"
for(i in 1:length(result_dir_correlation)){
  if(result_dir_correlation[i]){
    LISTE[i,"cor_dir"] <- "correlated"
  }else if(!result_dir_correlation[i]){
    LISTE[i,"cor_dir"] <- "anti_correlated"
  }
}

# Create file to save the result
write.table (LISTE,"../resultats/LISTE.txt", sep="\t",row.names=FALSE,quote=F)
