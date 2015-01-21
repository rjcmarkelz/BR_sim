#simulations.R
# quick code to play with some QTL assumptions for a simulated brassica rapa population.
library(qtl)

# simulate brassica cross
# this is for version 1.2 NOT 1.5 of the map
n.ind <- 124
genlen <- c(87, 110, 137, 61, 85, 91, 102, 75, 144, 102)
brass_markers <- c(136, 100, 212, 89, 72, 134, 139, 71, 177, 143)

mymap <- sim.map(len = genlen, n.mar = brass_markers,
	eq.spacing = FALSE, include.x = FALSE, anchor.tel=TRUE)
plot(mymap)
summary(mymap)

br_test_cross <- sim.cross(map = mymap, n.ind = 124, type = "riself")
summary(br_test_cross)

genotypes <- pull.geno(br_test_cross)
head(genotypes)
geno.names <- dimnames(genotypes)[[2]]

# Demo values using brassica marker simulated data
# sample significant markers for our 5 simulated phenotypes
m1 <- sample(geno.names, 3, replace = FALSE)
m2 <- sample(geno.names, 2, replace = FALSE)
m3 <- sample(geno.names, 2, replace = FALSE)
m4 <- sample(geno.names, 1, replace = FALSE)

## get marker genotypes
g11 <- genotypes[,m1[1]]; g12 <- genotypes[,m1[2]]; g13 <- genotypes[,m1[3]]
g21 <- genotypes[,m2[1]]; g22 <- genotypes[,m2[2]]
g31 <- genotypes[,m3[1]]; g32 <- genotypes[,m3[2]]
g41 <- genotypes[,m4[1]]; g42 <- genotypes[,m4[2]]

# generate correlated phenotypes, this is where you can generate any number of phenotypes that are or 
# are not correlated with one another
y1 <- runif(3,0.5,1)[g11] + runif(3,0.5,1)[g12] + rnorm(n.ind)
y2 <- runif(3,0.5,1)[g21] + runif(3,0.5,1)[g22] + rnorm(n.ind)
y3 <- runif(1,0.5,1) * y1 + runif(1,0.5,1) * y2 + runif(3,0.5,1)[g31] + runif(3,0.5,1)[g32] + rnorm(n.ind)
y4 <- runif(1,0.5,1) * y3 + runif(3,0.5,1)[g41] + runif(3,0.5,1)[g42] + rnorm(n.ind)

# take a quick look
cor(y1,y2)
cor(y2,y3)
cor(y3,y4)
cor(y2,y4)
cor(y1,y3)