breed = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_breedresident_comdat.csv", row.names = "X")
nonbreed = read.csv("data/derivative_data/gridcell_data/birds_comdat/birds_nonbreedresident_comdat.csv", row.names = "X")

breed.r = rast(breed, type = "xyz", crs = "+proj=cea +datum=WGS84")
nonbreed.r = rast(nonbreed, type = "xyz", crs = "+proj=cea +datum=WGS84")

plot(breed.r$ses.vert - nonbreed.r$ses.vert, main = "SES verticality:\nbreeding + resident - non-breeding + resident")
plot(breed.r$vert - nonbreed.r$vert, main = "Mean verticality:\nbreeding + resident - non-breeding + resident")
plot(breed.r$p.arb - nonbreed.r$p.arb, main = "% arboreal:\nbreeding + resident - non-breeding + resident")
plot(breed.r$ses.body.size - nonbreed.r$ses.body.size, main = "SES body size:\nbreeding + resident - non-breeding + resident")
