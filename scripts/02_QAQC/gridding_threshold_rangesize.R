library(sf)
amph = st_read("./../data/ranges/AMPHIBIANS/AMPHIBIANS.shp")
colnames(amph)

# What number of ranges cover 30% of a 111km x 111km gridcell
sum(amph$SHAPE_Area >= (111/3*111/3))

