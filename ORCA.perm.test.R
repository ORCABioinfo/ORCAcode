## drug sensitivity data ## 
load('drug.sensitivity')
load('no.member')
source('orca.functions.R')

dst.rank <- cor(t(drug.data),method='spearman',use='pairwise.complete.obs')

### calculate hypergeometric p-values : within-groups correlations  ###

group <- create.group(no.member)

tt.tick <- seq(0.1,by=0.1,0.9)
tt.vec <- seq(0.01,by=0.01,0.99)

shannon.entropy.vec <- vector()
phyper.list <- list()
phyper.matrix.list <- list()

source('phyperg.actual.R')
source('permutation.R')

p.permutation <- phyper.matrix.list
save(p.permutation,file='p.permutation')

########################################################################
## extract p-values from matrix ## 

require(gtools)

groups.comb <- combinations(length(group),2,seq(1,length(group),1))

pvec <- vector(mode='numeric',length=length(phyper.matrix.list))
p.list <- list()
for (i in 1:length(groups.comb[,1])){
	for (j in 1:length(phyper.matrix.list)){
		pvec[j] <- p.permutation[[j]][groups.comb[i,][1],groups.comb[i,][2]]
		}
	p.list[[i]] <- pvec
	}

save(p.list,file='p.list')
