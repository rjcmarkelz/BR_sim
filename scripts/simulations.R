#simulations.R
# quick code to play with some QTL assumptions for a simulated brassica rapa population.
library(qtl)

#some of this involves random draws
set.seed(123567)
# simulate brassica cross
# this is for version 1.2 NOT 1.5 of the map
n.ind <- 124
genlen <- c(87, 110, 137, 61, 85, 91, 102, 75, 144, 102)
brass_markers <- c(136, 100, 212, 89, 72, 134, 139, 71, 177, 143)

mymap <- sim.map(len = genlen, n.mar = brass_markers,
	eq.spacing = FALSE, include.x = FALSE, anchor.tel=TRUE)
plot(mymap)
summary(mymap)

?sim.cross
br_test_cross <- sim.cross(map = mymap, n.ind = 124, type = "riself")
summary(br_test_cross)

genotypes <- pull.geno(br_test_cross)
head(genotypes)
str(genotypes)
genotypes[genotypes == 2] <- 0
head(genotypes)
geno.names <- dimnames(genotypes)[[2]]

# Demo values using brassica marker simulated data
# sample significant markers for our 4 simulated phenotypes
# or define them as you wish
m1 <- sample(geno.names, 3, replace = FALSE)
m2 <- sample(geno.names, 2, replace = FALSE)
m3 <- sample(geno.names, 2, replace = FALSE)
m4 <- sample(geno.names, 1, replace = FALSE)

## get marker genotypes
g11 <- (genotypes[,m1[1]])
g11
g12 <- (genotypes[,m1[2]]) 
g13 <- (genotypes[,m1[3]])
g21 <- (genotypes[,m2[1]])
g22 <- (genotypes[,m2[2]])
g31 <- (genotypes[,m3[1]])
g32 <- (genotypes[,m3[2]])
g41 <- (genotypes[,m4[1]]) 

# function (orginal by Micheal Kuhn) to convert 0, 1 genotype assignments at markers into QTL with 
# tunable effect size and variance
br_traits <- function(geno, effect = 1, variance = 0.1){
 variance <- variance*abs(effect)
 geno <- as.logical(geno)
 pheno <- numeric(length(geno))
 pheno[!geno] = rnorm(sum(!geno), 0, variance)
 pheno[geno] = rnorm(sum(geno), effect, variance)
 return(pheno)
}

y1 <- br_traits(g11, effect = 2, variance = 0.1) + br_traits(g12, effect = 2, variance = 0.1)
y1

# replace trait columns what trait names
br_test_cross$pheno$trait1 <- br_traits(g11, effect = 2, variance = 0.1) + br_traits(g12, effect = 2, variance = 0.1)
br_test_cross$pheno


scanout <- scanone(br_test_cross, pheno.col = 2)
plot(scanout)
# this works to build whatever QTL model you want!


# another way to make correlated phenotypes
# generate correlated phenotypes, this is where you can generate any number of phenotypes that are or 
# are not correlated with one another
# 
# y1 <- runif(3,0.5,1)[g11] + runif(3,0.5,1)[g12] + rnorm(n.ind)
# y2 <- runif(3,0.5,1)[g21] + runif(3,0.5,1)[g22] + rnorm(n.ind)
# y3 <- runif(1,0.5,1) * y1 + runif(1,0.5,1) * y2 + runif(3,0.5,1)[g31] + runif(3,0.5,1)[g32] + rnorm(n.ind)
# y4 <- runif(1,0.5,1) * y3 + runif(3,0.5,1)[g41] + rnorm(n.ind)

# g11traits
# plot(g11traits)
# # take a quick look
# cor(y1,y2)
# cor(y2,y3)
# cor(y3,y4)
# cor(y2,y4)
# cor(y1,y3)

# # start with empty list
# allqtls <- list()
# m1.pos <- find.markerpos(br_test_cross, m1)
# m2.pos <- find.markerpos(br_test_cross, m2)
# m3.pos <- find.markerpos(br_test_cross, m3)
# m4.pos <- find.markerpos(br_test_cross, m4)
# ?makeqtl

# # you can choose to do more draws etc.
# br_test_cross <- sim.geno(br_test_cross, n.draws = 8, step = 2, err = 0.001)

# # you can put all of these into a list, or do one at a time
# # if you do one at a time, then the regular RQTL functions will work
# # leave it up to you to decide
# allqtls[[1]] <- makeqtl(br_test_cross, chr = m1.pos[,"chr"], pos = m1.pos[,"pos"])
# allqtls[[2]] <- makeqtl(br_test_cross, chr = m2.pos[,"chr"], pos = m2.pos[,"pos"])
# allqtls[[3]] <- makeqtl(br_test_cross, chr = m3.pos[,"chr"], pos = m3.pos[,"pos"])
# allqtls[[4]] <- makeqtl(br_test_cross, chr = m4.pos[,"chr"], pos = m4.pos[,"pos"])
# names(allqtls) <- c("trait1","trait2","trait3","trait4")

# # or if you want to fit multi qtl models directly
# multiqtl <- makeqtl(br_test_cross, c(chr = m1.pos[,"chr"], m2.pos[,"chr"]),
#                     pos = c(m1.pos[,"pos"], m2.pos[,"pos"]))
# summary(multiqtl)
# plot(multiqtl)