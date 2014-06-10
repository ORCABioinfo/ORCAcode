### functions for ORCA ##

## creat.group is a function that produce a list containing the two index numbers for each group. take a vector of numbers of the member in each group
create.group <- function(vec){
	temp <- list()
	x <- 0
	y <- 0
	for (i in 1:length(vec)){
		temp[[i]] <- c(x+1,y+vec[i])
		x <- x+vec[i]
		y <- y+vec[i]
		}
	temp
	}

create.xvec <- function(vec){
	temp <- vector()
	x <- 0
	for (i in 1:length(vec)){
		temp[i] <- x+1
		x <- x+vec[i]
		}
	temp
	}

create.yvec <- function(vec){
	temp <- vector()
	y <- 0
	for (i in 1:length(vec)){
		temp[i] <- y+vec[i]
		y <- y+vec[i]
		}
	temp
	}

## drugs is a function that calculates the number of drugs in each group according to the number of drugs in each group

drugs <- function(group){
	drug.num <- vector()
	for (i in 1:length(group)){
		drug.num[i] <- (group[[i]][2] - group[[i]][1]) + 1
		}
	drug.num
}
## all.btwg is a function that counts corr coef between two drug groups
# make 9 x 9 matrix 

all.btwg <- function(drug.num){
	all.btwg.cor.matrix <- matrix(nrow=length(group),ncol=length(group))
	for (i in 1:length(group)){
		all.btwg.cor.matrix[i,] <- drug.num * t(drug.num[i])
		}
	all.btwg.cor.matrix
}

# all.cpt is a function that calculates the number of cor coef which pass the threshold, theta (tt)

all.cpt <- function(tt,cormatrix) (length(which(abs(cormatrix) > tt)) - length(cormatrix[1,]))/2

# wig.cpt is a function that calculates the number of cor coef which pass the threshold, theta (tt), 
# WITHIN a given between group of correlation 

wig.cpt <- function(tt,cormatrix,group) {
  temp <- vector()
  for (i in 1:length(group)){
    temp[i] <- length(which(abs(cormatrix[group[[i]][1]:group[[i]][2],group[[i]][1]:group[[i]][2]]) > tt))
    temp[i] <- (temp[i] - ((group[[i]][2] - group[[i]][1])+1))/2
  }
  sum(temp)
}

## btwg.cpt is a fucntion that counts corr coef that PASS THRESHOLD between two drug groups
# BeTWeenGroup.cpt
# make 9 x 9 matrix

btwg.cpt <- function(tt,cormatrix){
	btwg.cpt.matrix <- matrix(nrow=length(group),ncol=length(group))

	for (i in 1:length(group)) {
		for (j in 1:length(group)){
			btwg.cpt.matrix[i,j] <- length(which(abs(cormatrix[group[[i]][1]:group[[i]][2],group[[j]][1]:group[[j]][2]]) > tt))
			}
		}
	btwg.cpt.matrix
}

## shannon's entropy calculation function 
## from https://stat.ethz.ch/pipermail/r-help/2008-July/167112.html 

shannon.entropy <- function(p){
	if (min(na.omit(p)) < 0 || sum(na.omit(p)) <= 0)
		return(NA)
	p.norm <- p[p>0]/sum(p)
	-sum(log2(p.norm)*p.norm)
	}


