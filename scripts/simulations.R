#simulations.R
# quick code to play with some QTL assumptions for a simulated brassica rapa population.
library(qtl)

# simulate brassica cross
n.ind <- 124
?sim.map
genlen <- c(87, 110, 137, 61, 85, 91, 102, 75, 144, 102)
brass_markers <- c(136, 100, 212, 89, 72, 134, 139, 71, 177, 143)

mymap <- sim.map(len = genlen, n.mar = brass_markers,
	eq.spacing = FALSE, include.x = FALSE, anchor.tel=TRUE)
plot(mymap)
summary(mymap)

mycross <- sim.cross(map = mymap, n.ind = 124, type = "riself")
summary(mycross)


genotypes <- pull.geno(mycross)
head(genotypes)
geno.names <- dimnames(genotypes)[[2]]