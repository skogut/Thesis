# some light exploartion of other directions for thesis research
# objectives: does the pred mod work on non caucasain females?
# do we have transcript details for TNFa and INFy? = yes
genedetails= read.csv("/Users/skogut/Desktop/Thesis/annotations/gene_details.csv")
genedetails$gene= as.character(genedetails$gene)
View(genedetails)
which(genedetails$gene == " IFNG")
genedetails[4396,2]
which(genedetails$gene == " HERVK")
genedetails[69554,]

deltarsq= read.csv("/Users/skogut/Desktop/Thesis/deltaRsq.csv")
View(deltarsq)
deltarsq$Transcript
for ( i in deltarsq$Transcript[1:210]){
  print(i)
}


# -------------------------------------------------------------------------
# full genome LM vs weighted ranks (iteration 2 will involve normalized unranked ERVs)
# files 
count= read.csv("/Users/skogut/Desktop/Thesis/All Cornell Counts and Results/status_sex_cornell_normalized_count_caucasian_genes.csv")
View(count)
ERVs= read.csv("/Users/skogut/Desktop/Thesis/weightedrank2.csv")
View(ERVs)
pheno= read.csv("/Users/skogut/Desktop/Thesis/cornell rna-seq data/cornell_pheno.csv")
pheno$sample_id= gsub("-", "_", pheno$sample_id)
caucfemale= which(pheno$sex== "F" & pheno$race== "w")
caucfemalepheno= pheno[caucfemale,]
IDcaucfemale=caucfemalepheno$s_name

# -------------------------------------------------------------------------
# creating the dataset
rm(count2)
names(count)= gsub("X", "", names(count))
count2= cbind(count[,1], count[,names(count)  %in% IDcaucfemale])
row.names(count2)= count2[,1]
count2= count2[,-1]
count2= t(count2)
View(count2)
write.csv(count2, "/Users/skogut/Desktop/Thesis/transposedcount2.csv" )
# now have normalized counts= combine
ERVrankedALL= cbind(ERVs,count2)
View(ERVrankedALL)
ERVrankedALL <- ERVrankedALL[,-2]
write.csv(ERVrankedALL, "/Users/skogut/Desktop/Thesis/ERVvgenome.csv")
x= read.csv("/Users/skogut/Desktop/Thesis/ERVvgenome.csv")

Results= as.data.frame(matrix(nrow=141*ncol(ERVrankedALL), ncol=5))
LMresults= as.data.frame(matrix(nrow=141, ncol= 5))
colnames(LMresults)= c("ERV", "Transcript", "Multiple R-squared", "Adjusted R-squared", "Slope Estimate")
colnames(Results)= c("ERV", "Transcript", "Multiple R-squared", "Adjusted R-squared", "Slope Estimate")
View(LMresults)


##### HAVE YOU SAVED RECENTLY???????? 
for (j in 143:ncol(ERVrankedALL)){
  for (i in 2:142) {
    print(colnames(ERVrankedALL)[j])
    print(colnames(ERVrankedALL)[i])
    EC= lm(ERVrankedALL[,j]~ ERVrankedALL[,i])
    
    LMresults$ERV[i-1]= colnames(ERVrankedALL)[i]
    LMresults$Transcript[i-1]= colnames(ERVrankedALL)[j]
    LMresults$`Multiple R-squared`[i-1]=summary(EC)$r.squared
    LMresults$`Adjusted R-squared`[i-1]= summary(EC)$adj.r.squared
    LMresults$`Slope Estimate`[i-1]= EC$coefficients[2]
  }
  row1= (i-1)*141+1
  Results[row1:(row1+140),]= LMresults
}
View(Results)
na.omit(Results)




