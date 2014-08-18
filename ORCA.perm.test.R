############################################################################################################
## OVER-REPRESENTATION OF CORRELATION (ORCA): A METHOD FOR IDENTIFYING ASSOCIATIONS BETWEEN VARIABLE SETS ## 
############################################################################################################

##############################################################################################
## 																							##
## The input data for ORCA can be any numerical data that have enough data points for  		##
## correlation coefficient calculation. The actual analysis uses two input data types: 		##
## 		1. Correlation matrix of all the variables calculated from the raw numerical data   ##
## 		2. Groupings or categories of the variables in the correlation matrix				##
## 																							##
## The input data for the examples we provided in the manuscripts are as follow:			##
## 		1. Drug sensitivity dataset															##
##			1.1 Correlation matrix between growth inhibition at 50% concentrations 			##
## 				(GI50 values) of 116 chemotherapeutic drugs on 58 cell lines from NCI60 		##
##				cell panle. 																	##
##			1.2 Groupings are the mechanisms of action of the drugs.							##
##		2. microRNA datasets																	##
## 			2.1 Correlation matrix between baseline expression of microRNAs on 58 cell  		## 
## 				lines from NCI60 cell panel. 												##
## 			2.2 Groupings are the microRNA clusters based on miRConnect study 				##
## 				(Hua et al. 2011).															##
## 																							##
## 			The output from ORCA is a hypergeometric p-value matrix calculated by 			##
## 			hypergeometric test between every group-pair. 									##
## 																							##
##############################################################################################

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
