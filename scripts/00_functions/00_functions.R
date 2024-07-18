library(BAMMtools)
library(fmsb)
library(moments)
library(mFD)

# Verticality -------------------------------------------------------------

# Add verticality data to presence absence data in site x spp matrix

# get_vert_dat = function(verticality, occr) {
#   # verticality: dataframe with verticality data, must have columns named sciname and Verticality
#   # occr: site x species matrix - first two colums should have x and y coordinates
#   # trait = name of trait in datafrom whose values are used to populate the prex/abs raster
#   
#   # make dataframe to store cell x species verticality matrix
#   occr_vert = occr %>% as.data.frame()
#   colnames(occr_vert) = gsub(" ", "_", colnames(occr_vert))
#   
#   for(i in 1:nrow(spp)) {
#     # get species
#     sp = spp[i,"sciname"]
#     # get verticality score for the species
#     v = verticality[which(verticality$sciname == sp), "Verticality"]
#     # replace occurrence 1 with verticality score
#     sp_occr = occr_vert[,sp]
#     sp_occr[sp_occr==1] = v
#     occr_vert[,sp] = sp_occr
#   }
#   
#   colnames(occr_vert)[1:2] = c("x", "y")
#   
#   occr_vert = as.data.table(occr_vert)
#   
#   #richness per cell
#   occr_vert[,richness:= apply(.SD, 1, function(x)sum(x>0)), .SDcols = 3:nrow(spp)]
#   
#   # mean verticality per cell
#   occr_vert[,vert_mean:= rowSums(.SD), .SDcols = 3:nrow(spp)]
#   occr_vert[,vert_mean:= vert_mean/richness]
#   
#   # number of semi-arboreal and arboreal species in each cell
#   occr_vert[,richness_arb:= apply(.SD, 1, function(x)sum(x>2)), .SDcols = 3:nrow(spp)]
#   
#   # proportion of semi-arboreal and arboreal species per cell
#   occr_vert[,prop_arb:= richness_arb/richness]
#   
#   # number of semi-fossorial and fossorial species per cell
#   occr_vert[,richness_fos:= apply(.SD, 1, function(x)sum(x<2 & x>0)), .SDcols = 3:nrow(spp)]
#   
#   # proportion of semi-fossorial and fossorial species per cell
#   occr_vert[,prop_fos := richness_fos/richness]
#   
#   # number of terrestrial species per cell
#   occr_vert[,richness_ter:= apply(.SD, 1, function(x)sum(x==2)), .SDcols = 3:nrow(spp)]
#   
#   # proportion of terrestrial species per cell
#   occr_vert[,prop_ter := richness_ter/richness]
#   
#   vert_dat = occr_vert[,.SD, .SDcols = c(1:2,(ncol(occr)+1):(ncol(occr)+8))]
#   
#   return(vert_dat)
#   
# }


# subset to species in occurrence data


# Get species level data --------------------------------------------------

get_sp_dat = function(trait_dat, filter_col=NA, filter_val=NA, occ, env, traitbase = NA){
  # trait_dat: dataframe; species trait information. Columns related to diet must include the word "Diet"
  # occ: dataframe with presence absence data and xy values for cells
  # env: env data with values for xy in occ dataframe. xy values in occ and env MUST match
  
  if(!is.na(filter_col)){
    sub = which(trait_dat[,filter_col] == filter_val)
    trait_dat = trait_dat[sub,]
  }
  
  
  # range size (n grid cells occupied)
  trait_dat$range_size = colSums(subset(occ, select = trait_dat$sciname))
  
  # food plasticity: number of different food categories eaten by species
  food.cols = grep("Diet", colnames(trait_dat))
  if(length(food.cols != 0)) {
    trait_dat$food.plast = apply(trait_dat[,food.cols], 1, function(x) length(which(x>0)))
  }
  
  # vertical plasticity: number of different vertical strata in which a species occurs
  if(traitbase == "elton") {
    trait_dat$vert.plast = apply(subset(trait_dat, select = ForStrat.watbelowsurf:ForStrat.aerial), 
                                 1, function(x) length(which(x>0)))
  } else {
    trait_dat$vert.plast = apply(subset(trait_dat, select = c("Fos", "Ter", "Aqu", "Arb", "Aer")), 1, function(x) length(which(x>0)))
  }
  
  
  # mean temperature within species range
  trait_dat$mat = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "mat"], na.rm = T))
  
  # mean min temperature of coldest month within species range
  trait_dat$tmin_cold = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "tmin_cold"], na.rm = T))
 
  # max temperature within species range
  trait_dat$tmax_warm = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "tmax_warm"], na.rm = T))
  
  # species` temperature niche width (max of max temp warmest month - min of min temp coldest month)
  trait_dat$temp_width = apply(subset(occ, select = trait_dat$sciname), 2, 
                               function(x) max(env[which(x==1), "tmax_warm"], na.rm = T) - min(env[which(x==1), "tmin_cold"], na.rm = T))
  
  # mean seasonality in temperature within species range
  trait_dat$temp_sea = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "temp_sea"], na.rm = T))
  
  # mean annual precip within species range
  trait_dat$precip_ann = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "precip_ann"], na.rm = T))
  
  # mean precip dryest month within species range
  trait_dat$precip_dry = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "precip_dry"], na.rm = T))
  
  # mean precip wettest month within species range
  trait_dat$precip_wet = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "precip_wet"], na.rm = T))
  
  # mean seasonality in precip within species range
  trait_dat$precip_sea = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "precip_sea"], na.rm = T))

  # mean humidity within species range
  trait_dat$hurs_mean = apply(subset(occ, select = trait_dat$sciname), 2, function(x) min(env[which(x==1), "hurs_mean"], na.rm = T))
  
  # mean annual range in humidity within species range
  trait_dat$hurs_range = apply(subset(occ, select = trait_dat$sciname), 2, function(x) min(env[which(x==1), "hurs_range"], na.rm = T))
  
  # mean vpd within species range
  trait_dat$vpd_mean = apply(subset(occ, select = trait_dat$sciname), 2, function(x) min(env[which(x==1), "vpd_mean"], na.rm = T))
  
  # mean annual range in vpd within species range
  trait_dat$vpd_range = apply(subset(occ, select = trait_dat$sciname), 2, function(x) min(env[which(x==1), "vpd_range"], na.rm = T))
  
  # max elevation within species range
  # warning occurs when all values are NA - silence warning
  trait_dat$elev_max = suppressWarnings(apply(subset(occ, select = trait_dat$sciname), 2, 
                                              function(x) max(env[which(x==1), "elev"])))
  
  # min elevation within species range
  trait_dat$elev_min = suppressWarnings(apply(subset(occ, select = trait_dat$sciname), 2, 
                                              function(x) min(env[which(x==1), "elev"], na.rm = T)))
  
  # mean elevation within species range
  trait_dat$elev_mean = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "elev"], na.rm = T))
  
  # range in elevation within species range
  trait_dat$elev_range = trait_dat$elev_max - trait_dat$elev_min
  
  # mean climate change velocity within species range
  trait_dat$clim_velocity = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "clim_velocity"], na.rm = T))
  
  # mean tree height within species` range
  trait_dat$canopy_height = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "canopy_height"], na.rm = T))
  
  # mean tree height within species` range
  trait_dat$veg_den = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "veg_den"], na.rm = T))
  
  # mean veg complexity within species` range
  trait_dat$veg_complexity = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "veg_complexity"]))
  
  # mean latitude within species range
  trait_dat$lat.mean = apply(subset(occ, select = trait_dat$sciname), 2, function(x) mean(env[which(x==1), "y"]))
  
  # when all values for a variable are NA, returns Inf, -Inf, or NaN
  # convert Inf, -Inf, or NaN to NA
  trait_dat = trait_dat %>% 
    mutate_all(~ifelse(is.nan(.), NA, .))
  trait_dat[trait_dat == -Inf] = NA
  trait_dat[trait_dat == Inf] = NA
  
  return(trait_dat)
}

# Calculate Functional Regularity
fro_calc = function(tr, occ) {
  # calculate functional regularity from presence absence matrix
  
  # check that species names are same and in same order in tr and occ
  check = sum(rownames(tr) != colnames(occ))
  if (check > 0) {
    stop("Species names in occurrence dataframe and trait dataframe do not match")
  }
  
  # which species are present in each grid cell
  sp = apply(occ, 1, function(x) which(x > 0))
  
  # get traits of species that are present
  t = lapply(sp, function(x) tr[x,])
  
  # order traits from smallest to largest
  t = lapply(t, sort)
  
  # calculate EW vector (assume abundances are all 1)
  dif = lapply(t, diff)
  ew = lapply(dif, "/", 2)
  
  # calculate PEW
  pew = lapply(ew, function(x) x/sum(x))
  
  # calculate fro
  fro = sapply(pew, function(x)sum(pmin(x,1/(length(x)))))
  
  return(fro)
}


# Get community level data for each grid cell -----------------------------

## Testing data
# rich = rowSums(occr[,3:ncol(occr)])
# head(which(rich > 5))
# which(rich>5)[1000:1005]
# 
# which(occr[836,]>0)
# which(occr[837,]>0)
# which(occr[838,]>0)
# which(occr[839,]>0)
# which(occr[840,]>0)
# which(occr[841,]>0)
# which(occr[2945,3:ncol(occr)]>0)
# # 131,426,1014,1391,1570,1573,1704,1705,2435,2982,2983
# 
# # testing data
# trait_dat = trait[c(131,426,1014,1391,1570,1573,1704,1705,2435,2982,2983),]
# occ = occr[c(836:841, 2945:2947), c(1:2, c(131,426,1014,1391,1570,1573,1704,1705,2435,2982,2983)+2)]
# env = read.csv("data/derivative_data/env_data.csv")
# occtest = occr[c(836:841, 2945:2947),]
# occ = occtest
# 
# trait_sub = trait_dat
# trait_dat = trait
# # for testing purposes, subset occ data AFTER having obtained trait data
# occ = occ[c(836:841, 2945:2947),]

get_gridcell_dat = function(trait_dat, filter_col=NA, filter_val=NA, occ, env, rich_min = 5) {
  # trait_dat: datatable or dataframe with taxonomy, vertical position, and traits
  ##     species names should be in column named "sciname"  
  ##     This should be the whole taxonomic group so that we can get taxonomic subset with respect to whole group
  ##     Must have columns diu and noc with 1 indicating yes for each
  # filter_col = character; column to filter data taxonomically
  # filter_val = character; value to filter to
  # occ: gridcell x species dataframe (presence = 1, absence = 0)
  ##       first two columns are x, y; subsequent columns are species names
  ##       species names must be in same order as scinames in trait_dat
  # env: dataframe with x,y coordinates that match occ and environmental data
  # rich_min: numeric; remove gridcells with fewer than this number of species - will not be included in analysis
  
  # test if species names in trait_dat are in same order as occ
  if(sum(colnames(occ)[3:ncol(occ)] != trait_dat$sciname) > 0) {
    stop("Species names in occurrence dataframe and trait dataframe do not match")
  }
  
  # if filter arguments have values, then filter the trait data to the subset
  if(!is.na(filter_col)) {
    sub = which(trait_dat[,filter_col] == filter_val)
    trait_sub = trait_dat[sub,]
  } else{
    trait_sub = trait_dat
  }
  
   # test if occ is a datatable and convert if not
  if (class(occ)[1] != "data.table"){
    occ = as.data.table(occ)
  }
  
  # add env data to occurrence data
  occ = left_join(occ, env, by = c("x", "y"))
  
  # Richness
  occ[,richness := apply(.SD, 1, sum), .SDcols = which(colnames(occ) %in% trait_dat$sciname)]
  # remove rows with richness less than rich min
  occ = occ[richness >= 5,]

  # subset occurrence data columns to only include rows with species of interest
  # community data will be calculated only including species of interest
  
  # which columns in occ contain species of interest
  # if filter_val is NA, then trait_sub will equal trait_dat
  sp_cols = which(colnames(occ) %in% trait_sub$sciname)
  
  # calculate richness of species of interest
  occ[,richness_sub := apply(.SD, 1, sum), .SDcols = sp_cols]
  
  # which columns in occ contain env data and richness data
  env_cols = which(colnames(occ) %in% c(colnames(env), "richness", "richness_sub"))
  
  # subset rows with richness of group of interest > 0 and columns with species in group of interest
  occ_sub = occ[richness_sub>0, .SD, .SDcols = c(env_cols, sp_cols)]

  # subset occurrence data to species of interest
  #occ = subset(occ, select = c("x","y", trait_dat$sciname))
  
  #xy = occ[,1:2]
  #occ = occ[,3:ncol(occ)]
  
  
  # Calculate community data for the subset of species requested
  # get column indices for species in occ_sub
  sp_cols = which(colnames(occ_sub) %in% trait_sub$sciname)
  
  # Calculate richness of arboreal, terrestrial, and fossorial species per cell
  # If species have more than one classification, count them in both
  # Number arboreal species
  # .SDcols selects column index of species in trait_dat
  occ_sub[,arb := apply(.SD, 1, function(x) length(which(trait_sub$Arb[which(x==1)] == 1))), .SDcols = sp_cols]
  
  # Number terrestrial species
  occ_sub[,ter := apply(.SD, 1, function(x) length(which(trait_sub$Ter[which(x==1)] == 1))), .SDcols = sp_cols]
  
  # Number fossorial species
  occ_sub[,fos := apply(.SD, 1, function(x) length(which(trait_sub$Fos[which(x==1)] == 1))), .SDcols = sp_cols]
  
  # Number diurnal species
  occ_sub[,diu := apply(.SD, 1, function(x) length(which(trait_sub$Diu[which(x==1)] == 1))), .SDcols = sp_cols]
  
  # Number nocturnal species
  occ_sub[,noc := apply(.SD, 1, function(x) length(which(trait_sub$Noc[which(x==1)] == 1))), .SDcols = sp_cols]
  
  # Calculate proportion of arboreal, terrestrial, and fossorial species
  # sum of the three proportions can be greater than 1 because species can occupy more than one vertical stratum
  # Proportion arboreal species
  occ_sub[,p.arb := arb/richness_sub]
  
  # Proportion terrestrial species
  occ_sub[,p.ter := ter/richness_sub]
  
  # Proportion fossorial species
  occ_sub[,p.fos := fos/richness_sub]
  
  # Proportion diurnal species
  occ_sub[,p.diu := diu/richness_sub]
  
  # Proportion nocturnal species
  occ_sub[,p.noc := noc/richness_sub]
  
  # Proportion of species in group of interest out of species in whole group
  occ_sub[,p.richsub := richness_sub/richness]
  
  
  
  # Mean community verticality
  occ_sub[,vert := apply(.SD, 1, function(x) mean(trait_sub$Verticality[which(x==1)])), .SDcols = sp_cols]
  
  # SD community verticality
  occ_sub[,vert.sd := apply(.SD, 1, function(x) sd(trait_sub$Verticality[which(x==1)])), .SDcols = sp_cols]
  
  # skewness community verticality
  occ_sub[,vert.skew := apply(.SD, 1, function(x) skewness(trait_sub$Verticality[which(x==1)])), .SDcols = sp_cols]
  
  # kurtosis community verticality
  occ_sub[,vert.kurtosis := apply(.SD, 1, function(x) kurtosis(trait_sub$Verticality[which(x==1)])), .SDcols = sp_cols]
  
  
  # Mean community body size
  occ_sub[,body.size := apply(.SD, 1, function(x) mean(trait_sub$body_size[which(x==1)])), .SDcols = sp_cols]
  occ_sub[,log.body.size := apply(.SD, 1, function(x) mean(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_cols]
  
  # SD body size
  occ_sub[,body.size.sd := apply(.SD, 1, function(x) sd(trait_sub$body_size[which(x==1)])), .SDcols = sp_cols]
  occ_sub[,log.body.size.sd := apply(.SD, 1, function(x) sd(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_cols]
  
  # kurtosis body size
  occ_sub[,log.body.size.kurtosis := apply(.SD, 1, function(x) kurtosis(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_cols]
  
  # kurtosis body size
  occ_sub[,log.body.size.skew := apply(.SD, 1, function(x) skewness(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_cols]
  
  # Mean community range size
  occ_sub[,range.size := apply(.SD, 1, function(x) mean(trait_sub$range_size[which(x==1)])), .SDcols = sp_cols]
  
  # Mean community food plasticity
  if(sum(colnames(trait_sub) %in% "food.plast")==1) {
    occ_sub[,food.plast := apply(.SD, 1, function(x) mean(trait_sub$food.plast[which(x==1)])), .SDcols = sp_cols]
  }
  
  # Mean vertical plasticity
  occ_sub[,vert.plast := apply(.SD, 1, function(x) mean(trait_sub$vert.plast[which(x==1)])), .SDcols = sp_cols]
  
  # Mean community thermal niche width
  occ_sub[,temp.width := apply(.SD, 1, function(x) mean(trait_sub$temp_width[which(x==1)])), .SDcols = sp_cols]
  
  # Mean range size
  occ_sub[,range.size := apply(.SD, 1, function(x) mean(trait_sub$range_size[which(x==1)])), .SDcols = sp_cols]
  
  # Mean elevation range
  occ_sub[,elev.range := apply(.SD, 1, function(x) mean(trait_sub$elev_range[which(x==1)])), .SDcols = sp_cols]
  
  # Univariate functional diversity
  # verticality
  verttraits = trait_sub[,c("Arb", "Fos", "Ter", "Aqu")]
  rownames(verttraits) = trait_sub$sciname
  
  # get functional entities for verticality
  fe.vert.fuzzy = sp.to.fe(sp_tr = verttraits,
                           tr_cat = data.frame(trait_name = c("Arb", "Fos", "Ter", "Aqu"), 
                                               trait_type = c("F", "F", "F", "F"), fuzzy_name = c("vert", "vert", "vert", "vert")),
                           fe_nm_type = "fe_rank", check_input = T)
  
  # FUNCTIONAL ENTITIES FOR DIEL ACTIVITY
  dieltraits = trait_sub[, c("Noc", "Diu")] 
  rownames(dieltraits) = trait_sub$sciname
  
  fe.diel.fuzzy = sp.to.fe(sp_tr = dieltraits %>% mutate_all(.funs = function(x){as.numeric(x)}),
                           tr_cat = data.frame(trait_name = c("Noc", "Diu"), 
                                               trait_type = c("F", "F"), fuzzy_name = c("diel", "diel")),
                           fe_nm_type = "fe_rank", check_input = T)
  
  
  
  occ_sub.df = occ_sub[,.SD, .SDcols = sp_cols]
  occ_sub.df = as.data.frame(occ_sub.df)
  
  alpha.fd.vert.fuzzy = alpha.fd.fe(asb_sp_occ = occ_sub.df, sp_to_fe = fe.vert.fuzzy,
                                    ind_nm = c("fored", "fred"), check_input = T, details_returned = T)
  
  alpha.fd.diel.fuzzy = alpha.fd.fe(asb_sp_occ = occ_sub.df, sp_to_fe = fe.diel.fuzzy,
                                    ind_nm = c("fored", "fred"), check_input = T, details_returned = T)
  
  
  # BODY SIZE DATA FRAME
  # calculate log body size
  traits_size = trait_sub %>% dplyr::select(body_size) %>% 
    mutate(log_size = log(body_size)) %>% 
    dplyr::select(!body_size)
  rownames(trait_size) = trait_sub$sciname
  
  # functional regularity of log(body size)
  fro_size = fro_calc(tr = traits_size, occ = occ_sub.df)
  
  # SES calculations
  # get species pools for different realms
  realms = unique(comdat$realm)
  realms = realms[which(!is.na(realms))]
  
  
  
  
  
  
  
  # Standardized effect size (SES) of metrics with respect to all species in full trait_dat dataframe
  # 1. randomly pick n species from the same REALM and calculate the null metric
  # 2. repeat 100x
  # 3. calculate the mean and sd null verticality
  
  # occ_sub only contains species of interest (filter_col == filter_val)
  ncell = nrow(occ_sub)
  
  if(sum(colnames(trait_sub) %in% "food.plast")==1){
    ses.vert = ses.vert.sd = ses.vert.plast = ses.food.plast = ses.body.size = ses.log.body.size = ses.log.body.size.sd = NA
    if(!is.na(filter_val)) {
      ses.vert.sub = ses.vert.sd.sub = ses.vert.plast.sub = ses.food.plast.sub = ses.body.size.sub = ses.log.body.size.sub = ses.log.body.size.sd.sub = NA
    }
  } else {
    ses.vert = ses.vert.sd = ses.vert.plast = ses.body.size = ses.log.body.size = ses.log.body.size.sd = NA
    if(!is.na(filter_val)) {
      ses.vert.sub = ses.vert.sd.sub = ses.vert.plast.sub = ses.body.size.sub = ses.log.body.size.sub = ses.log.body.size.sd.sub = NA
  }}
  

  # takes approx. 92 min
  for (i in 1:ncell){
    cat("\r", i, "of", ncell)
    
    #cells.sp <- occ[i,which(colnames(occ) %in% trait_sub$sciname)] # index of species in trait_sub that occur in the cell i
    realm_i = occ_sub[i, realm] # identify realm of the given cell
    
    if (is.na(realm_i)) {
      ifelse((sum(colnames(trait_sub) %in% "food.plast")==1),
             (ses.vert[i] = ses.vert.sd[i] = ses.food.plast[i] = ses.vert.plast[i] = ses.body.size[i] = ses.log.body.size[i] = ses.log.body.size.sd[i] = NA), 
             (ses.vert[i] = ses.vert.sd[i] = ses.vert.plast[i] = ses.body.size[i] = ses.log.body.size[i] = ses.log.body.size.sd[i] = NA)) 
      
      if(!is.na(filter_val)){
        ifelse((sum(colnames(trait_sub) %in% "food.plast")==1), 
               (ses.vert.sub[i] = ses.vert.sd.sub[i] = ses.food.plast.sub[i] = ses.vert.plast.sub[i] = ses.body.size.sub[i] = ses.log.body.size.sub[i] = ses.log.body.size.sd.sub[i] = NA),
               (ses.vert.sub[i] = ses.vert.sd.sub[i] = ses.vert.plast.sub[i] = ses.body.size.sub[i] = ses.log.body.size.sub[i] = ses.log.body.size.sd.sub[i] = NA))
      }
    } else{
      # subset occurrence dataframe to realm_i
      occ.realm = occ[realm == realm_i]
      
      # get index of all species in full trait dataset that occur in realm_i
      ## sum columns tells which species occur in occ.realm
      ## then columns with colsums > 0
      ## colnames in occ.realm of species == trait_dat$sciname, so we can use the pool to index species in trait_dat
      pool = which(colSums(occ.realm[,.SD, .SDcols = which(colnames(occ.realm) %in% trait_dat$sciname)]) > 0) 
      
      # get verticality scores of species in realm_i: trait_dat$vert[pool]
      # get random sample of rich[i] verticality scores in realm i (rich[i] = cell richness)
      ## pull sample from full pool of all species, not just species in trait subsest (i.e., order or family subset)
      # take mean of the random sample of verticality scores
      # repeat this 100 times to create null distribution
      # then calculate SES
      
      null.vert <- replicate(100, mean(sample(trait_dat$Verticality[pool], occ_sub[i, richness_sub], replace = F)))
      ses.vert[i] = (occ_sub[i, vert] - mean(null.vert)) / sd(null.vert)
      
      null.vert.sd <- replicate(100, sd(sample(trait_dat$Verticality[pool], occ_sub[i, richness_sub], replace = F)))
      ses.vert.sd[i] = (occ_sub[i, vert.sd] - mean(null.vert.sd)) / sd(null.vert.sd)
      
      if(sum(colnames(trait_sub) %in% "food.plast")==1) {  
          null.food.plast <- replicate(100, mean(sample(trait_dat$food.plast[pool], occ_sub[i, richness_sub], replace = F)))
          ses.food.plast[i] = (occ_sub[i, food.plast] - mean(null.food.plast)) / sd(null.food.plast)
      }
      
      null.vert.plast <- replicate(100, mean(sample(trait_dat$vert.plast[pool], occ_sub[i, richness_sub], replace = F)))
      ses.vert.plast[i] = (occ_sub[i, vert.plast] - mean(null.vert.plast)) / sd(null.vert.plast)
        
      null.body.size <- replicate(100, mean(sample(trait_dat$body_size[pool], occ_sub[i, richness_sub], replace = F)))
      ses.body.size[i] = (occ_sub[i, body.size] - mean(null.body.size)) / sd(null.body.size)
      
      null.log.body.size <- replicate(100, mean(log(sample(trait_dat$body_size[pool], occ_sub[i, richness_sub], replace = F))))
      ses.log.body.size[i] = (occ_sub[i, log.body.size] - mean(null.log.body.size)) / sd(null.log.body.size)
        
      null.log.body.size.sd <- replicate(100, sd(log(sample(trait_dat$body_size[pool], occ_sub[i, richness_sub], replace = F))))
      ses.log.body.size.sd[i] = (occ_sub[i, log.body.size.sd] - mean(null.log.body.size.sd)) / sd(null.log.body.size.sd)
      
      if(!is.na(filter_val)) {
        # repeat using pool of species from the group of interest (i.e. SES with respect to all species in group of interest)
       
        # subset occurrence dataframe to realm_i
        # occ_sub contains only gridcells with species belonging to group of interest
        occ.realm.sub = occ_sub[realm == realm_i]
        
        # get index of all species in subset data (filter_col==filter_val) that occur in realm_i
        ## sum columns tells which species occur in occ.realm
        ## then columns with colsums > 0
        ## colnames in occ.realm of species == trait_dat$sciname, so we can use the pool to index species in trait_dat
        pool.sub = which(colSums(occ.realm.sub[,.SD, .SDcols = which(colnames(occ.realm.sub) %in% trait_sub$sciname)]) > 0) 
        
        # get verticality scores of species in realm_i: trait_dat$vert[pool]
        # get random sample of rich[i] verticality scores in realm i (rich[i] = cell richness)
        ## pull sample from full pool of all species, not just species in trait subsest (i.e., order or family subset)
        # take mean of the random sample of verticality scores
        # repeat this 100 times to create null distribution
        # then calculate SES
        
        null.vert.sub <- replicate(100, mean(sample(trait_sub$Verticality[pool.sub], occ_sub[i, richness_sub], replace = F)))
        ses.vert.sub[i] = (occ_sub[i, vert] - mean(null.vert.sub)) / sd(null.vert.sub)
        
        null.vert.sd.sub <- replicate(100, sd(sample(trait_sub$Verticality[pool.sub], occ_sub[i, richness_sub], replace = F)))
        ses.vert.sd.sub[i] = (occ_sub[i, vert.sd] - mean(null.vert.sd.sub)) / sd(null.vert.sd.sub)
        
        if(sum(colnames(trait_sub) %in% "food.plast")==1) {
            null.food.plast.sub <- replicate(100, mean(sample(trait_sub$food.plast[pool.sub], occ_sub[i, richness_sub], replace = F)))
            ses.food.plast.sub[i] = (occ_sub[i, food.plast] - mean(null.food.plast.sub)) / sd(null.food.plast.sub)
        }
        
        null.vert.plast.sub <- replicate(100, mean(sample(trait_sub$vert.plast[pool.sub], occ_sub[i, richness_sub], replace = F)))
        ses.vert.plast.sub[i] =  (occ_sub[i, vert.plast] - mean(null.vert.plast.sub)) / sd(null.vert.plast.sub)
        
        null.body.size.sub <- replicate(100, mean(sample(trait_dat$body_size[pool.sub], occ_sub[i, richness_sub], replace = F)))
        ses.body.size.sub[i] =  (occ_sub[i, body.size] - mean(null.body.size.sub)) / sd(null.body.size.sub)
        
        null.log.body.size.sub <- replicate(100, mean(log(sample(trait_dat$body_size[pool.sub], occ_sub[i, richness_sub], replace = F))))
        ses.log.body.size.sub[i] = (occ_sub[i, log.body.size]) - mean(null.log.body.size.sub) / sd(null.log.body.size.sub)
        
        null.log.body.size.sd.sub <- replicate(100, sd(log(sample(trait_dat$body_size[pool.sub], occ_sub[i, richness_sub], replace = F))))
        ses.log.body.size.sd.sub[i] = (occ_sub[i, log.body.size.sd]) - mean(null.log.body.size.sd.sub) / sd(null.log.body.size.sd.sub)
        
      }
      
    }
  }
  
  occ_sub[,ses.vert := ses.vert]
  occ_sub[,ses.vert.sd := ses.vert.sd]
  if(sum(colnames(trait_sub) %in% "food.plast")==1){
    occ_sub[,ses.food.plast := ses.food.plast]
  }
  occ_sub[,ses.vert.plast := ses.vert.plast]
  occ_sub[,ses.body.size := ses.body.size]
  occ_sub[,ses.log.body.size := ses.log.body.size]
  occ_sub[,ses.log.body.size.sd := ses.log.body.size.sd]
  
  if(!is.na(filter_val)) {
    occ_sub[,ses.vert.sub := ses.vert.sub]
    occ_sub[,ses.vert.sd.sub := ses.vert.sd.sub]
    if(sum(colnames(trait_sub) %in% "food.plast")==1){
      occ_sub[,ses.food.plast.sub := ses.food.plast.sub]
    }
    occ_sub[,ses.vert.plast.sub := ses.vert.plast.sub]
    occ_sub[,ses.body.size.sub := ses.body.size.sub]
    occ_sub[,ses.log.body.size.sub := ses.log.body.size.sub]
    occ_sub[,ses.log.body.size.sd.sub := ses.log.body.size.sd.sub]
    
  }
  occ_sub = select(occ_sub, !trait_sub$sciname)
  return(occ_sub)

}
  


# get equal breaks from vector --------------------------------------------

### function to get equal binds from a vector
func_splint <- function(x,interval=4) {
  require(ggplot2)
  is.odd <- function(x) x %% 2 != 0
  
  a <- levels(cut_interval(x,interval-1))
  b <- unlist(strsplit(a,','))
  c <- gsub('[','',b,fixed="TRUE")
  d <- gsub(']','',c,fixed="TRUE")
  e <- gsub('(','',d,fixed="TRUE")
  f <- gsub(')','',e,fixed="TRUE")
  return(as.numeric(c(unique(f))))
}

# get_plot_breaks ------------------------------------------------------

# Get breaks for plotting legend
get_plot_breaks = function(com_vars, n) {
  # com_var: dataframe with gridcell metrics
  # Get class intervals RICH
  richness <- BAMMtools::getJenksBreaks(com_vars$richness, n)
  # avoid identical breaks
  richness <- unique(richness)
  if(length(richness)<n) {richness[(length(richness)+1):n] = NA}
  
  # Get class intervals VERT
  # For vert, I need 20 classes, but 0 must be in the middle
  # As a solution created Jenks for positive and Jenks for netative vert values
  vert1 <- BAMMtools::getJenksBreaks(com_vars$ses.vert[which(com_vars$ses.vert>0)], n/2)
  vert2 <- BAMMtools::getJenksBreaks(com_vars$ses.vert[which(com_vars$ses.vert<0)], n/2)
  ses.vert <- sort(c(vert1,vert2))
  # avoid identical breaks
  ses.vert <- unique(ses.vert)
  if(length(ses.vert)<n) {ses.vert[(length(ses.vert)+1):n] = NA}
  
  # Get class intervals ARB
  p.arb <- BAMMtools::getJenksBreaks(com_vars$p.arb, n)
  # avoid identical breaks
  p.arb <- unique(p.arb)
  if(length(p.arb)<n) {p.arb[(length(p.arb)+1):n] = NA}
  
  # Get class intervals TER
  p.ter <- BAMMtools::getJenksBreaks(com_vars$p.ter, n)
  # avoid identical breaks
  p.ter <- unique(p.ter)
  if(length(p.ter)<n) {p.ter[(length(p.ter)+1):n] = NA}
  
  # Get class intervals FOS
  p.fos <- BAMMtools::getJenksBreaks(com_vars$p.fos, n)
  # avoid identical breaks
  p.fos <- unique(p.fos) 
  p.fos[length(p.fos)] <- 1
  if(length(p.fos)<n) {p.fos[(length(p.fos)+1):n] = NA}
  
  # foscale <- round(func_splint(com_vars$p.fos,5),1)
  # foscale[5] <- 1
  
  # Get class intervals for body size and log body size
  body.size = BAMMtools::getJenksBreaks(com_vars$body.size, n)
  body.size = unique(body.size)
  if(length(body.size)<n) {body.size[(length(body.size)+1):n] = NA}
  
  log.body.size = BAMMtools::getJenksBreaks(com_vars$log.body.size, n)
  log.body.size = unique(log.body.size)
  if(length(log.body.size)<n) {log.body.size[(length(log.body.size)+1):n] = NA}
  
  # Get class intervals for SES body size
  ses.body.size = BAMMtools::getJenksBreaks(com_vars$ses.body.size, n)
  ses.body.size = unique(ses.body.size)
  if(length(ses.body.size)<n) {ses.body.size[(length(ses.body.size)+1):n] = NA}
  
  # return dataframe with breaks
  return(data.frame(richness = richness, ses.vert = ses.vert, p.arb = p.arb,
                    p.ter = p.ter, p.fos = p.fos, body.size = body.size, 
                    log.body.size = log.body.size, ses.body.size = ses.body.size))
}


# RGB plot ----------------------------------------------------------------

plot_vert_RGB = function(com_vars) {
  # com_vars: dataframe with cell metrics
  
  # subset to cells with richness at least 5
  com_vars = com_vars %>% filter(richness >= 5)
  
  # subset layers to x,y,p.arb, p.ter, and p.fos
  rstack = com_vars %>% dplyr::select(x, y, p.fos, p.arb, p.ter)
  # convert to raster
  rstack = rast(com_vars, type = "xyz")
  
  # Coordinates of the triangle
  tri <- rbind(sin(0:2*2/3*pi), cos(0:2*2/3*pi))

  # Function for calculating the color of a set of points `pt`
  # in relation to the triangle
  tricol <- function(pt, sharpness=2){
    require(splancs)
    RGB <- sapply(1:3, function(i){
      a <- sweep(pt, 2, tri[,i])
      b <- apply(tri[,-i], 1, mean) - tri[,i]
      sharpness*((a %*% b) / sum(b^2))-sharpness+1
    })
    RGB[-inpip(pt,t(tri)),] <- 1    # Color points outside the triangle white
    do.call(rgb, unname(as.data.frame(pmin(pmax(RGB, 0), 1))))
  }

  # Add triangle
  res <- 1000                         # Resolution
  xi <- seq(-1, 1, length=res)        # Axis points
  yi <- seq(-.8, 1.2, length=res)
  x <- xi[1] + cumsum(diff(xi))       # Midpoints between axis points
  y <- yi[1] + cumsum(diff(yi))
  xy <- matrix(1:(length(x)*length(y)), length(x))
  
  
  plot(s.mammals$Rich,main='Mammals strategies',cex.main=1.5,axes=F,box=F, col='white',legend=F)
  terra::plotRGB(rstack, r = "p.fos", g = "p.arb", b = "p.ter", stretch='hist', scale = 1, zlim = c(-1,1))
  map(mundi,add=T,cex=.5)
  boxed.labels(-13000000, 9000000, "B)", cex=1.5,bg="white",border=NA)
  
  # THIS IS WORKING!
  test = cell.df %>% filter(richness >= 5)
  test$vert.rgb = test$p.fos*256^2 + test$p.arb*256 + test$p.ter
  test$vert.rgbhex = rgb(red = test$p.fos, green = test$p.arb, blue = test$p.ter, maxColorValue = 1, alpha = 1)
  test.r = rast(test[,c("x", "y", "vert.rgb")], type = "xyz", crs = "+proj=cea +datum=WGS84")
  coltab(test.r, layer = 1) = test[,c("vert.rgb", "vert.rgbhex")]
  plot(test.r)
  # but still not showing stretched values
  # also why are there no fossorial mammals
  
    

}





# VIF variable selection --------------------------------------------------

#stepwise VIF function used below (Author: Brunno Oliveira)
vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  for(val in names(in_frame)){
    form_in<-formula(paste(val,' ~ .'))
    vif_init<-rbind(vif_init,c(val,VIF(lm(form_in,data=in_frame,...))))
  }
  
  vif_max = max(as.numeric(vif_init[,2]))
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(names(in_frame))
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      
      # calculate VIF for all variables
      for(val in names(in_dat)){
        form_in<-formula(paste(val,' ~ .'))
        vif_add<-VIF(lm(form_in,data=in_dat))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      
      # get row with max vif
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))
      
      # get max VIF
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      # remove variable with max vif from dataframe
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}




# GET GRIDCELL DATA (PARALLEL) --------------------------------------------



get_gridcell_dat_parallel = function(trait_dat, filter_col=NA, filter_val=NA, occ, env, rich_min = 5, ncore, nsim, eltonbirds = FALSE) {
  # trait_dat: datatable or dataframe with taxonomy, vertical position, and traits
  ##     species names should be in column named "sciname"  
  ##     This should be the whole taxonomic group so that we can get taxonomic subset with respect to whole group
  ##     Must have columns diu and noc with 1 indicating yes for each
  # filter_col = character; column to filter data taxonomically
  # filter_val = character; value to filter to
  # occ: gridcell x species dataframe (presence = 1, absence = 0)
  ##       first two columns are x, y; subsequent columns are species names
  ##       species names must be in same order as scinames in trait_dat
  # env: dataframe with x,y coordinates that match occ and environmental data
  # rich_min: numeric; remove gridcells with fewer than this number of species - will not be included in analysis
  
  # test if species names in trait_dat are in same order as occ
  if(sum(colnames(occ)[3:ncol(occ)] != trait_dat$sciname) > 0) {
    stop("Species names in occurrence dataframe and trait dataframe do not match")
  }
  
  # if filter arguments have values, then filter the trait data to the subset
  if(!is.na(filter_col)) {
    sub = which(trait_dat[,filter_col] == filter_val)
    trait_sub = trait_dat[sub,]
  } else{
    trait_sub = trait_dat
  }
  
  # test if occ is a datatable and convert if not
  if (class(occ)[1] != "data.table"){
    occ = as.data.table(occ)
  }
  
  
  # add env data to occurrence data
  occ = left_join(occ, env, by = c("x", "y"))
  
  # filter out realm == NA
  occ = occ %>% filter(!is.na(realm))
  
  # Richness
  occ[,richness := apply(.SD, 1, sum), .SDcols = which(colnames(occ) %in% trait_dat$sciname)]
  # remove rows with richness less than rich min
  occ = occ[richness >= 5,]
  
  
  # subset occurrence data columns to only include rows with species of interest
  # community data will be calculated only including species of interest
  
  # which columns in occ contain species of interest
  # if filter_val is NA, then trait_sub will equal trait_dat
  sp_cols = which(colnames(occ) %in% trait_sub$sciname)
  
  # calculate richness of species of interest
  occ[,richness_sub := apply(.SD, 1, sum), .SDcols = sp_cols]
  
  # which columns in occ contain env data and richness data
  env_cols = which(colnames(occ) %in% c(colnames(env), "richness", "richness_sub"))
  
  # subset rows with richness of group of interest > 0 and columns with species in group of interest
  occ_sub = occ[richness_sub>0, .SD, .SDcols = c(env_cols, sp_cols)]
  
  # subset occurrence data to species of interest
  #occ = subset(occ, select = c("x","y", trait_dat$sciname))
  
  #xy = occ[,1:2]
  #occ = occ[,3:ncol(occ)]
  
  
  # Calculate community data for the subset of species requested
  # get column indices for species in occ_sub
  #sp_cols = which(colnames(occ_sub) %in% trait_sub$sciname)
  
  # get number of cells each species occurs in
  sp_occ = apply(occ[,..sp_cols], 2, sum)
  sp_abs = which(sp_occ == 0)
  sp_abs = names(sp_abs)
  
  # remove species that do not occur in any cell
  # this could happen for rare species that only occurred in cells with fewer than 5 species
  occ_sub = occ_sub[,(sp_abs):=NULL]
  
  # which species occur in at least one cell
  # this list of species will be used to calculate community data on
  sp_pres = which(sp_occ >0)
  sp_pres = names(sp_pres)
  
  # remove species from trait dataframe that are no longer in occ dataframe
  trait_sub = trait_sub %>% filter(sciname %in% sp_pres)
  
  # testing to make sure order of species is still the same
  # sum(trait_sub$sciname != sp_pres)
  
  
  # Calculate richness of arboreal, terrestrial, and fossorial species per cell
  # If species have more than one classification, count them in both
  # Number arboreal species
  # .SDcols selects column index of species in trait_dat
  
  occ_sub[,arb := apply(.SD, 1, function(x) length(which(trait_sub$Arb[which(x==1)] == 1))), .SDcols = sp_pres]
  
  # Number terrestrial species
  occ_sub[,ter := apply(.SD, 1, function(x) length(which(trait_sub$Ter[which(x==1)] == 1))), .SDcols = sp_pres]
  
  # Number fossorial species
  occ_sub[,fos := apply(.SD, 1, function(x) length(which(trait_sub$Fos[which(x==1)] == 1))), .SDcols = sp_pres]
  
  # Number diurnal species
  occ_sub[,diu := apply(.SD, 1, function(x) length(which(trait_sub$Diu[which(x==1)] == 1))), .SDcols = sp_pres]
  
  # Number nocturnal species
  occ_sub[,noc := apply(.SD, 1, function(x) length(which(trait_sub$Noc[which(x==1)] == 1))), .SDcols = sp_pres]
  
  # Calculate proportion of arboreal, terrestrial, and fossorial species
  # sum of the three proportions can be greater than 1 because species can occupy more than one vertical stratum
  # Proportion arboreal species
  occ_sub[,p.arb := arb/richness_sub]
  
  # Proportion terrestrial species
  occ_sub[,p.ter := ter/richness_sub]
  
  # Proportion fossorial species
  occ_sub[,p.fos := fos/richness_sub]
  
  # Proportion diurnal species
  occ_sub[,p.diu := diu/richness_sub]
  
  # Proportion nocturnal species
  occ_sub[,p.noc := noc/richness_sub]
  
  # Proportion of species in group of interest out of species in whole group
  occ_sub[,p.richsub := richness_sub/richness]
  
  
  
  # Mean community verticality
  occ_sub[,vert.mean := apply(.SD, 1, function(x) mean(trait_sub$Verticality[which(x==1)])), .SDcols = sp_pres]
  
  # SD community verticality
  occ_sub[,vert.sd := apply(.SD, 1, function(x) sd(trait_sub$Verticality[which(x==1)])), .SDcols = sp_pres]
  
  # skewness community verticality
  occ_sub[,vert.skew := apply(.SD, 1, function(x) skewness(trait_sub$Verticality[which(x==1)])), .SDcols = sp_pres]
  
  # kurtosis community verticality
  occ_sub[,vert.kurtosis := apply(.SD, 1, function(x) kurtosis(trait_sub$Verticality[which(x==1)])), .SDcols = sp_pres]
  
  # Mean community nocturnality
  occ_sub[,nocturnality.mean := apply(.SD, 1, function(x) mean(trait_sub$Nocturnality[which(x==1)])), .SDcols = sp_pres]
  
  # Mean community body size
  occ_sub[,size.mean := apply(.SD, 1, function(x) mean(trait_sub$body_size[which(x==1)])), .SDcols = sp_pres]
  occ_sub[,log.size.mean := apply(.SD, 1, function(x) mean(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_pres]
  
  # SD body size
  occ_sub[,size.sd := apply(.SD, 1, function(x) sd(trait_sub$body_size[which(x==1)])), .SDcols = sp_pres]
  occ_sub[,log.size.sd := apply(.SD, 1, function(x) sd(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_pres]
  
  # kurtosis body size
  occ_sub[,log.size.kurtosis := apply(.SD, 1, function(x) kurtosis(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_pres]
  
  # kurtosis body size
  occ_sub[,log.size.skew := apply(.SD, 1, function(x) skewness(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_pres]
  
  # Mean community range size
  occ_sub[,range.size := apply(.SD, 1, function(x) mean(trait_sub$range_size[which(x==1)])), .SDcols = sp_pres]
  
  # Mean community food plasticity
  if(sum(colnames(trait_sub) %in% "food.plast")==1) {
    occ_sub[,food.plast := apply(.SD, 1, function(x) mean(trait_sub$food.plast[which(x==1)])), .SDcols = sp_pres]
  }
  
  # Mean vertical plasticity
  occ_sub[,vert.plast := apply(.SD, 1, function(x) mean(trait_sub$vert.plast[which(x==1)])), .SDcols = sp_pres]
  
  # Mean community thermal niche width
  occ_sub[,temp.width := apply(.SD, 1, function(x) mean(trait_sub$temp_width[which(x==1)])), .SDcols = sp_pres]
  
  # Mean range size
  occ_sub[,range.size := apply(.SD, 1, function(x) mean(trait_sub$range_size[which(x==1)])), .SDcols = sp_pres]
  
  # Mean elevation range
  occ_sub[,elev.range := apply(.SD, 1, function(x) mean(trait_sub$elev_range[which(x==1)])), .SDcols = sp_pres]
  
  # Univariate functional diversity
  # verticality
  
  if (eltonbirds) {
    verttraits = trait_sub %>% dplyr::select(ForStrat.watbelowsurf:ForStrat.aerial)
    rownames(verttraits) = trait_sub$sciname
    
    # get functional entities for verticality
    fe.vert.fuzzy = sp.to.fe(sp_tr = verttraits,
                             tr_cat = data.frame(trait_name = colnames(verttraits), 
                                                 trait_type = rep("F",7), fuzzy_name = rep("vert",7)),
                             fe_nm_type = "fe_rank", check_input = T)
    
  } else {
    verttraits = trait_sub[,c("Arb", "Fos", "Ter", "Aqu")]
    rownames(verttraits) = trait_sub$sciname
    
    # get functional entities for verticality
    fe.vert.fuzzy = sp.to.fe(sp_tr = verttraits,
                             tr_cat = data.frame(trait_name = c("Arb", "Fos", "Ter", "Aqu"), 
                                                 trait_type = c("F", "F", "F", "F"), fuzzy_name = c("vert", "vert", "vert", "vert")),
                             fe_nm_type = "fe_rank", check_input = T)
  }
  
  # FUNCTIONAL ENTITIES FOR DIEL ACTIVITY
  dieltraits = trait_sub[, c("Noc", "Diu")] 
  rownames(dieltraits) = trait_sub$sciname
  
  fe.diel.fuzzy = sp.to.fe(sp_tr = dieltraits %>% mutate_all(.funs = function(x){as.numeric(x)}),
                           tr_cat = data.frame(trait_name = c("Noc", "Diu"), 
                                               trait_type = c("F", "F"), fuzzy_name = c("diel", "diel")),
                           fe_nm_type = "fe_rank", check_input = T)
  
  # make dataframe of only occurrence data without env data
  occ_sub.sp = occ_sub[, ..sp_pres]
  occ_sub.sp = as.data.frame(occ_sub.sp)
  rownames(occ_sub.sp) = 1:nrow(occ_sub.sp)
  
  alpha.fd.vert.fuzzy = alpha.fd.fe(asb_sp_occ = occ_sub.sp, sp_to_fe = fe.vert.fuzzy,
                                    ind_nm = c("fored", "fred"), check_input = T, details_returned = T)
  occ_sub[, for.vert := alpha.fd.vert.fuzzy[[1]][,"fored"]]
  
  alpha.fd.diel.fuzzy = alpha.fd.fe(asb_sp_occ = occ_sub.sp, sp_to_fe = fe.diel.fuzzy,
                                    ind_nm = c("fored", "fred"), check_input = T, details_returned = T)
  occ_sub[, for.diel := alpha.fd.diel.fuzzy[[1]][,"fored"]]
  
  
  # BODY SIZE DATA FRAME
  # calculate log body size
  traits_size = trait_sub %>% dplyr::select(body_size) %>% 
    mutate(log_size = log(body_size)) %>% 
    dplyr::select(!body_size)
  rownames(traits_size) = trait_sub$sciname
  
  # functional regularity of log(body size)
  fro_size = fro_calc(tr = traits_size, occ = occ_sub.sp)
  occ_sub[, fro.log.size := fro_size]
  
  # SES calculations
  # get species pools for different realms
  realms = unique(occ_sub$realm)
  realms = realms[which(!is.na(realms))]
  
  # make species pools for each realm
  pools = list()
  for (i in 1:length(realms)) {
    # get realm 
    realmi = realms[i]
    
    # which cells occur in realmi
    realmi.cells = which(occ_sub$realm == realmi)
    
    # subset occurrence dataframe to the given realm
    realmi.occ = occ_sub[realmi.cells,]
    
    # which species occur in realmi
    realmi.sp = apply(realmi.occ[,.SD, .SDcols = sp_pres], 2, sum)
    realmi.sp = which(realmi.sp > 0)
    
    pools[[realmi]] = realmi.sp
  }
  
  # start cluster
  # NEED TO CHECK THIS FUNCTION FOR FILTERING
  # MAKE SURE IT IS RELATIVE TO ENTIRE POOL NOT JUST FILTERED SPECIES
  cl = makeCluster(ncore)
  registerDoParallel(cl)
  
  sims = foreach(n = 1:nsim, .packages = c("mFD", "dplyr"), .inorder = TRUE, .combine = "cbind", .export = "fro_calc") %dopar% {
    
    # make blank occ matrix for adding simulated communities
    occ.sim = occ_sub.sp
    occ.sim[] <- NA
    
    for(r in 1:nrow(occ.sim)) {
      # identify realm and richness for a given cell
      realmi = occ_sub[[r, "realm"]]
      richi = occ_sub[[r, "richness"]]
      
      # sample species from the pool
      sim = sample(pools[[realmi]], richi, replace = F)
      
      # add sampled species to the simulated occurrence dataframe
      occ.sim[r, names(sim)] <- 1
      
    }
    
    occ.sim[is.na(occ.sim)] <- 0
    
    # calculate functional diversity indices
    
    # check if all species are present in at least one assemblage
    check = sum(apply(occ.sim, 2, sum) == 0)
    if (check > 0){
      notpres = which(apply(occ.sim, 2, sum) == 0)
      occ.sim = occ.sim[,-notpres]
      
      tr2 = trait_sub[-notpres, ,drop = FALSE]
      
      tr_size2 = traits_size[-notpres, , drop = FALSE]
      
      if (eltonbirds) {
        vt2 = tr2 %>% dplyr::select(ForStrat.watbelowsurf:ForStrat.aerial)
        rownames(vt2) = tr2$sciname
        
        fe_vert2 = sp.to.fe(sp_tr = vt2,
                            tr_cat = data.frame(trait_name = colnames(vt2), 
                                                trait_type = rep("F", 7), fuzzy_name = rep("vert", 7)),
                            fe_nm_type = "fe_rank", check_input = T)
      } else {
         vt2 = data.frame(Arb = tr2$Arb, Fos = tr2$Fos, Ter = tr2$Ter, Aqu = tr2$Aqu)
          rownames(vt2) = tr2$sciname
      
          fe_vert2 = sp.to.fe(sp_tr = vt2,
                          tr_cat = data.frame(trait_name = c("Arb", "Fos", "Ter", "Aqu"), 
                                              trait_type = c("F", "F", "F", "F"), fuzzy_name = c("vert", "vert", "vert", "vert")),
                          fe_nm_type = "fe_rank", check_input = T)
      }
      
     
      
      dt2 = data.frame(Diu = tr2$Diu, Noc = tr2$Noc)
      rownames(dt2) = tr2$sciname
      
      fe_diel2 = sp.to.fe(sp_tr = dt2,
                          tr_cat = data.frame(trait_name = c("Diu", "Noc"), 
                                              trait_type = c("F", "F"), fuzzy_name = c("diel", "diel")),
                          fe_nm_type = "fe_rank", check_input = T)
      
      # functional regularity for log body size
      fro.size.sim = fro_calc(tr = tr_size2, occ = occ.sim)
      
      # function over-redundancy for verticality
      for.vert.sim = alpha.fd.fe(asb_sp_occ = occ.sim, sp_to_fe = fe_vert2,
                                 ind_nm = c("fored"), check_input = T, details_returned = T)

      # functional over redundancy for diel activity
      for.diel.sim = alpha.fd.fe(asb_sp_occ = occ.sim, sp_to_fe = fe_diel2,
                                 ind_nm = c("fored", "fred"), check_input = T, details_returned = T)
      
      vert.mean.sim = apply(occ.sim, 1, function(x) {mean(tr2[which(x==1), "Verticality"])})
      vert.sd.sim = apply(occ.sim, 1, function(x) {sd(tr2[which(x==1), "Verticality"])})
      
      size.mean.sim = apply(occ.sim, 1, function(x){mean(tr2[which(x==1), "body_size"])})
      size.sd.sim = apply(occ.sim, 1, function(x){sd(tr2[which(x==1), "body_size"])})
      
      log.size.mean.sim = apply(occ.sim, 1, function(x){mean(log(tr2[which(x==1), "body_size"]))})
      log.size.sd.sim = apply(occ.sim, 1, function(x){sd(log(tr2[which(x==1), "body_size"]))})
      
      nocturnality.mean.sim = apply(occ.sim, 1, function(x) {mean(tr2[which(x==1), "Nocturnality"])})
    
    } else {
      # functional regularity for log body size
      fro.size.sim = fro_calc(tr = traits_size, occ = occ.sim)  
      
      # functional over redundancy for verticality
      for.vert.sim = alpha.fd.fe(asb_sp_occ = occ.sim, sp_to_fe = fe.vert.fuzzy,
                                 ind_nm = c("fored"), check_input = T, details_returned = T)
      
      # functional over redundancy for diel activity
      for.diel.sim = alpha.fd.fe(asb_sp_occ = occ.sim, sp_to_fe = fe.diel.fuzzy,
                                 ind_nm = c("fored", "fred"), check_input = T, details_returned = T)
      
      vert.mean.sim = apply(occ.sim, 1, function(x) {mean(trait_sub[which(x==1), "Verticality"])})
      vert.sd.sim = apply(occ.sim, 1, function(x) {sd(trait_sub[which(x==1), "Verticality"])})
      
      size.mean.sim = apply(occ.sim, 1, function(x){mean(trait_sub[which(x==1), "body_size"])})
      size.sd.sim = apply(occ.sim, 1, function(x){sd(trait_sub[which(x==1), "body_size"])})
      
      log.size.mean.sim = apply(occ.sim, 1, function(x){mean(log(trait_sub[which(x==1), "body_size"]))})
      log.size.sd.sim = apply(occ.sim, 1, function(x){sd(log(trait_sub[which(x==1), "body_size"]))})
      
      nocturnality.mean.sim = apply(occ.sim, 1, function(x) {mean(trait_sub[which(x==1), "Nocturnality"])})
    }
    
    # concatenate results into a single vector
    simi = c(fro.size.sim, for.vert.sim[[1]][,"fored"], for.diel.sim[[1]][,"fored"],
             vert.mean.sim, vert.sd.sim,
             size.mean.sim, size.sd.sim,
             log.size.mean.sim, log.size.sd.sim,
             nocturnality.mean.sim)
    
    return(simi)
    
    # fro.size.null[,n] = fro.size.sim
    # for.vert.null[,n] = for.vert.sim[[1]][,"fored"]
    # for.diel.null[,n] = for.diel.sim[[1]][,"fored"]
    # toc()
    # print(n)
  }
  stopCluster(cl)
  
  ncells = nrow(occ_sub)
  # calculate SES
  fro.size.null = sims[1:ncells,]
  fro.size.mean = apply(fro.size.null, 1, mean)
  fro.size.sd = apply(fro.size.null, 1, sd)
  occ_sub[, fro.log.size.ses := (fro.log.size - fro.size.mean)/fro.size.sd]
  
  for.vert.null = sims[(ncells+1):(2*ncells),]
  for.vert.mean = apply(for.vert.null, 1, mean)
  for.vert.sd = apply(for.vert.null, 1, sd)
  occ_sub[, for.vert.ses := (for.vert - for.vert.mean)/for.vert.sd]
  
  for.diel.null = sims[(2*ncells+1):(3*ncells),]
  for.diel.mean = apply(for.diel.null, 1, mean)
  for.diel.sd = apply(for.diel.null, 1, sd)
  occ_sub[, for.diel.ses := (for.diel - for.diel.mean)/for.vert.sd]
  
  vert.mean.null = sims[(3*ncells+1):(4*ncells),]
  vert.mean.mean = apply(vert.mean.null, 1, mean)
  vert.mean.sd = apply(vert.mean.null, 1, sd)
  occ_sub[, vert.mean.ses := (vert.mean - vert.mean.mean)/vert.mean.sd]
  
  vert.sd.null = sims[(4*ncells+1):(5*ncells),]
  vert.sd.mean = apply(vert.sd.null, 1, mean)
  vert.sd.sd = apply(vert.sd.null, 1, sd)
  occ_sub[, vert.sd.ses := (vert.sd - vert.sd.mean)/vert.sd.sd]
  
  size.mean.null = sims[(5*ncells+1):(6*ncells),]
  size.mean.mean = apply(size.mean.null, 1, mean)
  size.mean.sd = apply(size.mean.null, 1, sd)
  occ_sub[, size.mean.ses := (size.mean - size.mean.mean)/size.mean.sd]
  
  size.sd.null = sims[(6*ncells+1):(7*ncells),]
  size.sd.mean = apply(size.sd.null, 1, mean)
  size.sd.sd = apply(size.sd.null, 1, sd)
  occ_sub[, size.sd.ses := (size.sd - size.sd.mean)/size.sd.sd]
  
  log.size.mean.null = sims[(7*ncells+1):(8*ncells),]
  log.size.mean.mean = apply(log.size.mean.null, 1, mean)
  log.size.mean.sd = apply(log.size.mean.null, 1, sd)
  occ_sub[, log.size.mean.ses := (log.size.mean - log.size.mean.mean)/log.size.mean.sd]
  
  log.size.sd.null = sims[(8*ncells+1):(9*ncells),]
  log.size.sd.mean = apply(log.size.sd.null, 1, mean)
  log.size.sd.sd = apply(log.size.sd.null, 1, sd)
  occ_sub[, log.size.sd.ses := (log.size.sd - log.size.sd.mean)/log.size.sd.sd]
  
  nocturnality.mean.null = sims[(9*ncells+1):(10*ncells),]
  nocturnality.mean.mean = apply(nocturnality.mean.null, 1, mean)
  nocturnality.mean.sd = apply(nocturnality.mean.null, 1, sd)
  occ_sub[, nocturnality.mean.ses := (nocturnality.mean - nocturnality.mean.mean)/nocturnality.mean.sd]
  
  return(occ_sub) 
  
  
  
}



# Get gridcell vert parallel ----------------------------------------------
get_gridcell_vert_parallel = function(trait_dat, filter_col=NA, filter_val=NA, occ, env, rich_min = 5, ncore, nsim, eltonbirds = FALSE) {
  # trait_dat: datatable or dataframe with taxonomy, vertical position, and traits
  ##     species names should be in column named "sciname"  
  ##     This should be the whole taxonomic group so that we can get taxonomic subset with respect to whole group
  ##     Must have columns diu and noc with 1 indicating yes for each
  # filter_col = character; column to filter data taxonomically
  # filter_val = character; value to filter to
  # occ: gridcell x species dataframe (presence = 1, absence = 0)
  ##       first two columns are x, y; subsequent columns are species names
  ##       species names must be in same order as scinames in trait_dat
  # env: dataframe with x,y coordinates that match occ and environmental data
  # rich_min: numeric; remove gridcells with fewer than this number of species - will not be included in analysis
  
  # test if species names in trait_dat are in same order as occ
  if(sum(colnames(occ)[3:ncol(occ)] != trait_dat$sciname) > 0) {
    stop("Species names in occurrence dataframe and trait dataframe do not match")
  }
  
  # if filter arguments have values, then filter the trait data to the subset
  if(!is.na(filter_col)) {
    sub = which(trait_dat[,filter_col] == filter_val)
    trait_sub = trait_dat[sub,]
  } else{
    trait_sub = trait_dat
  }
  
  # test if occ is a datatable and convert if not
  if (class(occ)[1] != "data.table"){
    occ = as.data.table(occ)
  }
  
  
  # add env data to occurrence data
  occ = left_join(occ, env, by = c("x", "y"))
  
  # filter out realm == NA
  occ = occ %>% filter(!is.na(realm))
  
  # Richness
  occ[,richness := apply(.SD, 1, sum), .SDcols = which(colnames(occ) %in% trait_dat$sciname)]
  # remove rows with richness less than rich min
  occ = occ[richness >= 5,]
  
  
  # subset occurrence data columns to only include rows with species of interest
  # community data will be calculated only including species of interest
  
  # which columns in occ contain species of interest
  # if filter_val is NA, then trait_sub will equal trait_dat
  sp_cols = which(colnames(occ) %in% trait_sub$sciname)
  
  # calculate richness for group of interest
  occ[,richness_sub := apply(.SD, 1, sum), .SDcols = sp_cols]
  
  # which columns in occ contain env data and richness data
  env_cols = which(colnames(occ) %in% c(colnames(env), "richness", "richness_sub"))
  
  # subset rows with richness of group of interest > 0 and columns with species in group of interest
  occ_sub = occ[richness_sub>0, .SD, .SDcols = c(env_cols, sp_cols)]
  
  # subset occurrence data to species of interest
  #occ = subset(occ, select = c("x","y", trait_dat$sciname))
  
  #xy = occ[,1:2]
  #occ = occ[,3:ncol(occ)]
  
  
  # Calculate community data for the subset of species requested
  # get column indices for species in occ_sub
  #sp_cols = which(colnames(occ_sub) %in% trait_sub$sciname)
  
  # get number of cells each species occurs in
  sp_occ = apply(occ[,..sp_cols], 2, sum)
  sp_abs = which(sp_occ == 0)
  sp_abs = names(sp_abs)
  
  # remove species that do not occur in any cell
  # this could happen for rare species that only occurred in cells with fewer than 5 species
  occ_sub = occ_sub[,(sp_abs):=NULL]
  
  # which species occur in at least one cell
  # this list of species will be used to calculate community data on
  sp_pres = which(sp_occ >0)
  sp_pres = names(sp_pres)
  
  # remove species from trait dataframe that are no longer in occ dataframe
  trait_sub = trait_sub %>% filter(sciname %in% sp_pres)
  
  # testing to make sure order of species is still the same
  # sum(trait_sub$sciname != sp_pres)
  
  
  # Calculate richness of arboreal, terrestrial, and fossorial species per cell
  # If species have more than one classification, count them in both
  # Number arboreal species
  # .SDcols selects column index of species in trait_dat
  
  occ_sub[,arb := apply(.SD, 1, function(x) length(which(trait_sub$Arb[which(x==1)] == 1))), .SDcols = sp_pres]
  
  # Number terrestrial species
  occ_sub[,ter := apply(.SD, 1, function(x) length(which(trait_sub$Ter[which(x==1)] == 1))), .SDcols = sp_pres]
  
  # Number fossorial species
  occ_sub[,fos := apply(.SD, 1, function(x) length(which(trait_sub$Fos[which(x==1)] == 1))), .SDcols = sp_pres]
  
  # Number diurnal species
  occ_sub[,diu := apply(.SD, 1, function(x) length(which(trait_sub$Diu[which(x==1)] == 1))), .SDcols = sp_pres]
  
  # Number nocturnal species
  occ_sub[,noc := apply(.SD, 1, function(x) length(which(trait_sub$Noc[which(x==1)] == 1))), .SDcols = sp_pres]
  
  # Calculate proportion of arboreal, terrestrial, and fossorial species
  # sum of the three proportions can be greater than 1 because species can occupy more than one vertical stratum
  # Proportion arboreal species
  occ_sub[,p.arb := arb/richness_sub]
  
  # Proportion terrestrial species
  occ_sub[,p.ter := ter/richness_sub]
  
  # Proportion fossorial species
  occ_sub[,p.fos := fos/richness_sub]
  
  # Proportion diurnal species
  occ_sub[,p.diu := diu/richness_sub]
  
  # Proportion nocturnal species
  occ_sub[,p.noc := noc/richness_sub]
  
  # Proportion of species in group of interest out of species in whole group
  occ_sub[,p.richsub := richness_sub/richness]
  
  
  
  # Mean community verticality
  occ_sub[,vert.mean := apply(.SD, 1, function(x) mean(trait_sub$Verticality[which(x==1)])), .SDcols = sp_pres]
  
  # SD community verticality
  occ_sub[,vert.sd := apply(.SD, 1, function(x) sd(trait_sub$Verticality[which(x==1)])), .SDcols = sp_pres]
  
  # skewness community verticality
  occ_sub[,vert.skew := apply(.SD, 1, function(x) skewness(trait_sub$Verticality[which(x==1)])), .SDcols = sp_pres]
  
  # kurtosis community verticality
  occ_sub[,vert.kurtosis := apply(.SD, 1, function(x) kurtosis(trait_sub$Verticality[which(x==1)])), .SDcols = sp_pres]
  
  # Mean community nocturnality
  occ_sub[,nocturnality.mean := apply(.SD, 1, function(x) mean(trait_sub$Nocturnality[which(x==1)])), .SDcols = sp_pres]
  
  # Mean community body size
  occ_sub[,size.mean := apply(.SD, 1, function(x) mean(trait_sub$body_size[which(x==1)])), .SDcols = sp_pres]
  occ_sub[,log.size.mean := apply(.SD, 1, function(x) mean(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_pres]
  
  # SD body size
  occ_sub[,size.sd := apply(.SD, 1, function(x) sd(trait_sub$body_size[which(x==1)])), .SDcols = sp_pres]
  occ_sub[,log.size.sd := apply(.SD, 1, function(x) sd(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_pres]
  
  # kurtosis body size
  occ_sub[,log.size.kurtosis := apply(.SD, 1, function(x) kurtosis(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_pres]
  
  # kurtosis body size
  occ_sub[,log.size.skew := apply(.SD, 1, function(x) skewness(log(trait_sub$body_size[which(x==1)]))), .SDcols = sp_pres]
  
  # Mean community range size
  occ_sub[,range.size := apply(.SD, 1, function(x) mean(trait_sub$range_size[which(x==1)])), .SDcols = sp_pres]
  
  # Mean community food plasticity
  if(sum(colnames(trait_sub) %in% "food.plast")==1) {
    occ_sub[,food.plast := apply(.SD, 1, function(x) mean(trait_sub$food.plast[which(x==1)])), .SDcols = sp_pres]
  }
  
  # Mean vertical plasticity
  occ_sub[,vert.plast := apply(.SD, 1, function(x) mean(trait_sub$vert.plast[which(x==1)])), .SDcols = sp_pres]
  
  # Mean community thermal niche width
  occ_sub[,temp.width := apply(.SD, 1, function(x) mean(trait_sub$temp_width[which(x==1)])), .SDcols = sp_pres]
  
  # Mean range size
  occ_sub[,range.size := apply(.SD, 1, function(x) mean(trait_sub$range_size[which(x==1)])), .SDcols = sp_pres]
  
  # Mean elevation range
  occ_sub[,elev.range := apply(.SD, 1, function(x) mean(trait_sub$elev_range[which(x==1)])), .SDcols = sp_pres]
  
  # SES calculations
  # get species pools for different realms
  realms = unique(occ_sub$realm)
  realms = realms[which(!is.na(realms))]
  
  # make species pools for each realm
  pools = list()
  for (i in 1:length(realms)) {
    # get realm 
    realmi = realms[i]
    
    # which cells occur in realmi
    realmi.cells = which(occ_sub$realm == realmi)
    
    # subset occurrence dataframe to the given realm
    realmi.occ = occ_sub[realmi.cells,]
    
    # which species occur in realmi
    realmi.sp = apply(realmi.occ[,.SD, .SDcols = sp_pres], 2, sum)
    realmi.sp = which(realmi.sp > 0)
    
    pools[[realmi]] = realmi.sp
  }
  
  
  # make dataframe of only occurrence data without env data
  occ_sub.sp = occ_sub[, ..sp_pres]
  occ_sub.sp = as.data.frame(occ_sub.sp)
  rownames(occ_sub.sp) = 1:nrow(occ_sub.sp)
  
  # start cluster
  # NEED TO CHECK THIS FUNCTION FOR FILTERING
  # MAKE SURE IT IS RELATIVE TO ENTIRE POOL NOT JUST FILTERED SPECIES
  cl = makeCluster(ncore)
  registerDoParallel(cl)
  
  sims = foreach(n = 1:nsim, .packages = c("dplyr"), .inorder = TRUE, .combine = "cbind") %dopar% {
    
    # make blank occ matrix for adding simulated communities
    occ.sim = occ_sub.sp
    occ.sim[] <- NA
    
    for(r in 1:nrow(occ.sim)) {
      # identify realm and richness for a given cell
      realmi = occ_sub[[r, "realm"]]
      richi = occ_sub[[r, "richness_sub"]] # richness for group of interest
      
      # sample species from the pool
      sim = sample(pools[[realmi]], richi, replace = F)
      
      # add sampled species to the simulated occurrence dataframe
      occ.sim[r, names(sim)] <- 1
    }
    
    occ.sim[is.na(occ.sim)] <- 0
    
    # calculate functional diversity indices
    # check if all species are present in at least one assemblage
    check = sum(apply(occ.sim, 2, sum) == 0)
    if (check > 0){
      notpres = which(apply(occ.sim, 2, sum) == 0)
      occ.sim = occ.sim[,-notpres]
      
      tr2 = trait_sub[-notpres, ,drop = FALSE]
      
      vert.mean.sim = apply(occ.sim, 1, function(x) {mean(tr2[which(x==1), "Verticality"])})
      vert.sd.sim = apply(occ.sim, 1, function(x) {sd(tr2[which(x==1), "Verticality"])})
      
      size.mean.sim = apply(occ.sim, 1, function(x){mean(tr2[which(x==1), "body_size"])})
      size.sd.sim = apply(occ.sim, 1, function(x){sd(tr2[which(x==1), "body_size"])})
      
      log.size.mean.sim = apply(occ.sim, 1, function(x){mean(log(tr2[which(x==1), "body_size"]))})
      log.size.sd.sim = apply(occ.sim, 1, function(x){sd(log(tr2[which(x==1), "body_size"]))})
      
      nocturnality.mean.sim = apply(occ.sim, 1, function(x) {mean(tr2[which(x==1), "Nocturnality"])})
      
    } else {
      vert.mean.sim = apply(occ.sim, 1, function(x) {mean(trait_sub[which(x==1), "Verticality"])})
      vert.sd.sim = apply(occ.sim, 1, function(x) {sd(trait_sub[which(x==1), "Verticality"])})
      
      size.mean.sim = apply(occ.sim, 1, function(x){mean(trait_sub[which(x==1), "body_size"])})
      size.sd.sim = apply(occ.sim, 1, function(x){sd(trait_sub[which(x==1), "body_size"])})
      
      log.size.mean.sim = apply(occ.sim, 1, function(x){mean(log(trait_sub[which(x==1), "body_size"]))})
      log.size.sd.sim = apply(occ.sim, 1, function(x){sd(log(trait_sub[which(x==1), "body_size"]))})
      
      nocturnality.mean.sim = apply(occ.sim, 1, function(x) {mean(trait_sub[which(x==1), "Nocturnality"])})
    }
    
    # concatenate results into a single vector
    simi = c(vert.mean.sim, vert.sd.sim,
             size.mean.sim, size.sd.sim,
             log.size.mean.sim, log.size.sd.sim,
             nocturnality.mean.sim)
    
    return(simi)
    
  }
  stopCluster(cl)
  
  ncells = nrow(occ_sub)
  # calculate SES
  vert.mean.null = sims[1:ncells,]
  vert.mean.mean = apply(vert.mean.null, 1, mean)
  vert.mean.sd = apply(vert.mean.null, 1, sd)
  occ_sub[, vert.mean.ses := (vert.mean - vert.mean.mean)/vert.mean.sd]
  
  vert.sd.null = sims[(1*ncells+1):(2*ncells),]
  vert.sd.mean = apply(vert.sd.null, 1, mean)
  vert.sd.sd = apply(vert.sd.null, 1, sd)
  occ_sub[, vert.sd.ses := (vert.sd - vert.sd.mean)/vert.sd.sd]
  
  size.mean.null = sims[(2*ncells+1):(3*ncells),]
  size.mean.mean = apply(size.mean.null, 1, mean)
  size.mean.sd = apply(size.mean.null, 1, sd)
  occ_sub[, size.mean.ses := (size.mean - size.mean.mean)/size.mean.sd]
  
  size.sd.null = sims[(3*ncells+1):(4*ncells),]
  size.sd.mean = apply(size.sd.null, 1, mean)
  size.sd.sd = apply(size.sd.null, 1, sd)
  occ_sub[, size.sd.ses := (size.sd - size.sd.mean)/size.sd.sd]
  
  log.size.mean.null = sims[(4*ncells+1):(5*ncells),]
  log.size.mean.mean = apply(log.size.mean.null, 1, mean)
  log.size.mean.sd = apply(log.size.mean.null, 1, sd)
  occ_sub[, log.size.mean.ses := (log.size.mean - log.size.mean.mean)/log.size.mean.sd]
  
  log.size.sd.null = sims[(5*ncells+1):(6*ncells),]
  log.size.sd.mean = apply(log.size.sd.null, 1, mean)
  log.size.sd.sd = apply(log.size.sd.null, 1, sd)
  occ_sub[, log.size.sd.ses := (log.size.sd - log.size.sd.mean)/log.size.sd.sd]
  
  nocturnality.mean.null = sims[(6*ncells+1):(7*ncells),]
  nocturnality.mean.mean = apply(nocturnality.mean.null, 1, mean)
  nocturnality.mean.sd = apply(nocturnality.mean.null, 1, sd)
  occ_sub[, nocturnality.mean.ses := (nocturnality.mean - nocturnality.mean.mean)/nocturnality.mean.sd]
  
  return(occ_sub) 

}

