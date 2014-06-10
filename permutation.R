### permutation for drug sensitivity data at correaltion cutoff = 0.79 ###

i=79

for (perm in 2:1000000) {
  ### permute drug labels ###
  drug.data.sample <- drug.data[sample(nrow(drug.data)),]
  dst.rank <- cor(t(drug.data.sample),method='spearman',use='pairwise.complete.obs')
  	
  # calculate within groups cor coef that pass the threshold tt (within groups being the cor coef along the diagonal line)
  withinGroup.cpt <- wig.cpt(tt.vec[i],dst.rank,group) 
  
  # calculate all cor coef that pass the thrshold tt
  acpt <- all.cpt(tt.vec[i],dst.rank)
  
  # calculate drug number in each drug group
  drug.num <- drugs(group)
  
  # calculate number of cor coef for each between group (all) and put in data matrix 
  all.btwg.cor.matrix <- all.btwg(drug.num)
  
  # calculate the number of cor coef that pass the threshold and put in data matrix 
  btwg.cpt.matrix <- btwg.cpt(tt.vec[i],dst.rank)
  
  # calculate the number of btwg.cpt without the within groups cpt 
  all.btwg.cpt <- acpt - withinGroup.cpt
  
  # calculate cor coef for each between group
  all.btwg.cor <- ((sum(all.btwg.cor.matrix) - 116)/2) - sum((drug.num^2 -drug.num)/2)
  
  # calculate hypergeometric p-value and put in data matrix
  # substitute the phyper() with the values from the matrices
  # phyper(q,m,n,k) 
  # q = # corr coef between groups that pass threshold
  # m = # corr coef between groups 
  # n = # corr coef between groups of the REST 
  # k = # ALL other corr coef betwen groups that pass the threshold tt
  phyper.matrix <- matrix(nrow=length(group),ncol=length(group))
  phyper.matrix <- phyper(btwg.cpt.matrix,all.btwg.cor.matrix,all.btwg.cor-all.btwg.cor.matrix,all.btwg.cpt,F)
  
  ### calculate hypergeometric p-value of the within-groups ###
  for (j in 1:length(group)) {
  	q <- (btwg.cpt.matrix[j,j] - drug.num[j])/2
  	m <- (all.btwg.cor.matrix[j,j] - drug.num[j])/2
  	n <- all.btwg.cor-m
  	k <- acpt
  	phyper.matrix[j,j] <- phyper(q,m,n,k,F)
  	}
   
  phyper.matrix.list[[perm]] <- phyper.matrix
}
