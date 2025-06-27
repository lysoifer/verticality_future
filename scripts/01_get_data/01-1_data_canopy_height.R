# Aggregate and average Lang et al. 2023 canopy height raster
# excluding non-forested areas/human-modified (e.g., oil palm)
# Non-forested areas are masked by Hansen et al. 2022 Global land use and land cover
# Canopy height data available from https://langnico.github.io/globalcanopyheight/
# Landcover data available from https://glad.umd.edu/dataset/global-land-cover-land-use-v1

# This script must be run on hipergator as canopy height and land cover data is stored there

library(terra)
library(numform)
library(doParallel)
library(foreach)

# canopy height indexed by lower left corner (ymin, xmin) (NXXEXXX_Map.tif)
ch = list.files("/orange/scheffers/dklinges/soiltemp/SoilTemp-big-data/vegetation/canopy_height/lang_et_al_2023_10m/", full.names = T)
ch = grep("Map.tif", ch, value = T)

# lc indexed by upper left corner (ymax, xmin) (XXNXXXE.tif)
lc = list.files("/orange/scheffers/lydia.soifer/big_data/land_use/hansen2022/tiles_strata/", full.names = T)
lctest = rast(lc[1])

# path to landcover tiles
lcpath = "/orange/scheffers/lydia.soifer/big_data/land_use/hansen2022/tiles_strata/GLCLU_2019_strata_"

# split list of canopy height tiles into 10 to parallelize
breaksNE = paste0("N", f_pad_zero(seq(0,87,3), 2), "E")
breaksSE = paste0("S", f_pad_zero(seq(0,87,3), 2), "E")
breaksNW = paste0("N", f_pad_zero(seq(0,87,3), 2), "W")
breaksSW = paste0("S", f_pad_zero(seq(0,87,3), 2), "W")

tilesubsNE = lapply(breaksNE, grep, ch, value = T)
tilesubsSE = lapply(breaksSE, grep, ch, value = T)
tilesubsNW = lapply(breaksNW, grep, ch, value = T)
tilesubsSW = lapply(breaksSW, grep, ch, value = T)

tilesubs = c(tilesubsNE, tilesubsNW, tilesubsSE, tilesubsSW)
lengthsubs = sapply(tilesubs, length)
tilesubs = tilesubs[lengthsubs>0] # remove groups with no tiles

# weird tile - N00W075 - there was an issue with the download of this tile in Dave's folder, so I have redownloaded it in the project data
t = 29
i = 9
t29i9 = rast("data/original/canopy_height_Lang2023/ETH_GlobalCanopyHeight_10m_2020_N00W075_Map.tif")

# try running 34-95 to see if I can catch the error in these
# tilesubs = tilesubs[34-95]
# 
# # subsetting tiles to the ones that errored
# tilesubs = tilesubs[c(2,1,3,5,6,28,29,30,41,42,48,49,
#                       53,54,55,56,50,58,59,57,71,72,
#                       73,74,77,75,80,76,78,79,87,88,89,83,
#                       84,92,93,94)]

cl = makeCluster(7)
registerDoParallel(cl)
foreach(t = 1:length(tilesubs), .packages = c("terra", "numform"), .errorhandling = "pass", .verbose = T) %dopar% {
  fname = paste0("/blue/scheffers/lydia.soifer/verticality_future/tempdir_chagg_nobareground/", f_pad_zero(t, 3), ".tif")
  
  if(!file.exists(fname)) {
    ttiles = tilesubs[[t]]
    tagg = list()
    for(i in 1:length(ttiles)) {
      # load canopy height raster and aggregate
      cat(sprintf("Running tile %d, index %d\n", t, i), 
          file = "/blue/scheffers/lydia.soifer/verticality_future/error-log2.txt", append = TRUE)
      c = rast(ttiles[i])
      if(t == 29 & i == 9) {
        c = rast("/blue/scheffers/lydia.soifer/verticality_future/data/original/canopy_height_Lang2023/ETH_GlobalCanopyHeight_10m_2020_N00W075_Map.tif")
      }
      c = aggregate(c, 3, fun = "mean")
      
      # test if raster is all NA values and skip if it is
      if(global(c, fun = "notNA")[1,1] != 0) {
        # identify bounds for reading in relevant landcover tiles
        xmin = ext(c)[1]
        xmax = ext(c)[2]
        ymin = ext(c)[3]
        ymax = ext(c)[4]
        
        xmin2 = floor(xmin/10)*10
        xmax2 = floor(xmax/10)*10
        ymin2 = ceiling(ymin/10)*10
        ymax2 = ceiling(ymax/10)*10
        
        xmin2.ew = ifelse(xmin2<0, "W", "E")
        xmin2 = ifelse(xmin2<0, xmin2*-1, xmin2)
        xmax2.ew = ifelse(xmax2<0, "W", "E")
        xmax2 = ifelse(xmax2<0, xmax2*-1, xmax2)
        ymin2.ns = ifelse(ymin2<0, "S", "N")
        ymin2 = ifelse(ymin2<0, ymin2*-1, ymin2)
        ymax2.ns = ifelse(ymax2<0, "S", "N")
        ymax2 = ifelse(ymax2<0, ymax2*-1, ymax2)
        
        # get tile names
        tiles = paste0(lcpath, f_pad_zero(ymax2,2), ymax2.ns, "_", f_pad_zero(xmin2,3), xmin2.ew, ".tif")
        if((xmin2 != xmax2) & (ymin2 == ymax2)) {tiles = c(tiles, paste0(lcpath, f_pad_zero(ymax2,2), ymax2.ns, "_", f_pad_zero(xmax2,3), xmin2.ew, ".tif"))}
        if((xmin2 == xmax2) & (ymin2 != ymax2)) {tiles = c(tiles, paste0(lcpath, f_pad_zero(ymin2,2), ymin2.ns, "_", f_pad_zero(xmin2,3), xmin2.ew, ".tif"))}
        if((xmin2 != xmax2) & (ymin2 != ymax2)) {tiles = c(tiles, paste0(lcpath, f_pad_zero(ymax2,2), ymax2.ns, "_", f_pad_zero(xmax2,3), xmax2.ew, ".tif"),
                                                           paste0(lcpath, f_pad_zero(ymin2,2), ymin2.ns, "_", f_pad_zero(xmin2,3), xmin2.ew, ".tif"),
                                                           paste0(lcpath, f_pad_zero(ymin2,2), ymin2.ns, "_", f_pad_zero(xmax2,3), xmax2.ew, ".tif"))}
        
        # check to make sure tiles exist (boundary case near lat lon extremes)
        tiles = tiles[file.exists(tiles)]
        # read in landcover tiles that cover canopy height tile
        # then crop to canopy height tile and merge if using multiple tiles
        
        if(length(tiles != 0)) {
          if(length(tiles)>1){tiles = sprc(tiles)} else {tiles = rast(tiles)}
          tiles = tryCatch({crop(tiles, c)}, 
                           error = function(e){
                             message("crop failed")
                             "failed"})
          #tiles = crop(tiles, c)
          if(is(tiles, "character")) { # deals with lat/lon edge cases
            cmask = c
          } else {
            if(is(tiles, "SpatRasterCollection")) {tiles = merge(tiles)}
            
            # mask values that are not natural land cover - see legend values in hansen2022 folder
            # redo excluding bare land as well
            #mv = c(6,7,13,14,15,16,17,18,19,20) # masking out human land uses
            mv = c(1,2,3,6,7,8,9,10,13,14,15,16,17,18,19,20) # mask out everything except in tact woody veg (> 3m tall)- this is better at capturing all deforested land
            if(ext(tiles) == ext(c)) {cmask = mask(c, tiles, maskvalue = mv, updatevalue = NA)} else {cmask = c} # prevent errors in north where lc is maps weren't produced
          }
        } else {
          cmask = c
        }
        tagg[[i]] = aggregate(cmask, 0.45/0.00025, fun = "median", na.rm = T)
        
        # clean diskspace
        rm(cmask, tiles)
        gc()
        
      }
      
      rm(c)
      gc()
      
    }
    
    # mosaic tiles and save
    tagg = tagg[!sapply(tagg, is.null)]
    mos = mosaic(sprc(tagg))
    x = round(xmin(mos))
    xew = ifelse(x<0, "W", "E")
    x = paste0(xew, f_pad_zero(abs(x),3))
    y = round(ymin(mos))
    yns = ifelse(y<0, "S", "N")
    y = paste0(yns, f_pad_zero(abs(y),2))
    
    
    writeRaster(mos, fname, overwrite = T)
  }
}

stopCluster(cl)

# test = list.files("/blue/scheffers/lydia.soifer/verticality_future/tempdir2/", full.names = T)
# chagg = sprc(list.files("/blue/scheffers/lydia.soifer/verticality_future/tempdir2/", full.names = T))
# chagg = mosaic(chagg, filename = "/blue/scheffers/lydia.soifer/verticality_future/data/derivative_data/canopy_height/canopy_height_agg.tif")
# 
# ch2 = rast("/orange/scheffers/lydia.soifer/big_data/vegetation/canopy_height_lang2023/canopy_height_reduceRes_mosaic.tif")
# 
test = list.files("/blue/scheffers/lydia.soifer/verticality_future/tempdir_chagg_nobareground/", full.names = T)
chagg = sprc(list.files("/blue/scheffers/lydia.soifer/verticality_future/tempdir_chagg_nobareground/", full.names = T))
# this removes all landcovers that don't have trees
chagg = mosaic(chagg, filename = "/blue/scheffers/lydia.soifer/verticality_future/data/derivative_data/canopy_height/canopy_height_agg_ForestOnly.tif", overwrite = T)
# use the layer that masks out bare ground and only model in forests because there is a lot of bare ground that looks like crops yet is not classified as so
# this should be a better representation of historical biogeographic patterns

# compare to canopy height computed with bare ground
chagg1 = rast("/blue/scheffers/lydia.soifer/verticality_future/data/derivative_data/canopy_height/canopy_height_agg.tif")

dif = chagg - chagg1
plot(dif)


