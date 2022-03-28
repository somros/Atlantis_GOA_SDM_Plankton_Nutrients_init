# 02/07/2022 this code calculates values of initial biomass of plankton groups and nutrients from ROMS-NPZ for virgin biomass calculations (init.nc).
# as of 02/07/2022 tjhis is based on Jan climatologies of ROMS NEP 10K data for 2017. We will need to use Jan data closest to the model initialzation date
# This code also breaks down detritus from NPZ to labile, refractory, and carrion
# It also converts NO3 and NH4+ concentrations (in millimoles typically) to their equivalents in mg N that Atlantis needs

library(tidyverse)
library(data.table)
library(sf)
library(rbgm)

select <- dplyr::select

#read in model geometry
atlantis_bgm <- read_bgm('../data/GOA_WGS84_V4_final.bgm')
atlantis_box <- atlantis_bgm %>% box_sf()

# reconstruct depth layers
z <- c(30,100,200,500,1000,4000) # for GOA

# list files
file_list <- list.files('../outputs/january_init/',full.names = TRUE)

# read files and bind them by row
npz_vars_list <- lapply(file_list, read.csv)

# write function to: 1. fill boxes with no NPZ data based on nearest neighbor; 2. set island boxes to 0
# I am wary of setting boundary boxes to 0, as advection from those boxes may be important in dynamic boxes

fill_boxes <- function(varframe){
  # need to set conc in island boxes as 0, and conc in empty boxes as the nearest conc. This takes a bit of reconstruction
  # repeat this for each variable to be sure although in theory it should be the same for all NPZ vars
  missing_boxes <- setdiff(0:108, varframe %>% pull(.bx0) %>% unique())
  island_boxes <- atlantis_box %>% filter(botz == 0) %>% pull(.bx0) %>% unique()
  no_data_boxes <- atlantis_box %>% filter(.bx0 %in% missing_boxes, botz != 0) %>% pull(.bx0) %>% unique()
  
  # calculate the nearest neighbor for each Atlantis box
  neighbor <- atlantis_box %>%
    st_set_geometry(NULL) %>%
    st_as_sf(coords = c(x='insideX',y='insideY')) 
  
  non_empty_neighbor <- neighbor %>% filter(.bx0 %in% setdiff(atlantis_box$.bx0,no_data_boxes))
  
  neighbor <- neighbor %>%
    mutate(nearest_empty = st_nearest_feature(geometry)-1) %>%
    rowwise() %>%
    mutate(index_nearest_not_empty = st_nearest_feature(geometry, non_empty_neighbor$geometry),
           nearest_not_empty = non_empty_neighbor[index_nearest_not_empty,]$.bx0,
           nearest = ifelse(nearest_empty %in% no_data_boxes, nearest_not_empty, nearest_empty))
  
  # make an empty list
  list_empty_boxes <- vector(mode = 'list', length = length(missing_boxes))
  
  for(i in 1:length(missing_boxes)){

    tb <- missing_boxes[i] 
    df <- data.frame('.bx0'= tb, 'botz' = -(atlantis_box %>% filter(.bx0==tb) %>% pull(botz)))
    
    if(df$botz==0){
      df1 <- data.frame(df, 'maxz' = 0)  
    } else {
      df1 <- data.frame(df, 'maxz' = c(z[z<=df$botz],df$botz))  
    }
    
    df1 <- df1 %>% select(-botz)
    
    df1 <- df1 %>% 
      distinct() %>%
      mutate(var_mean = NA,
             long_name = varframe$long_name[1],
             unit = varframe$unit[1])
    
    this_neighbor <- neighbor %>% filter(.bx0==tb) %>% pull(nearest)
    
    neighbor_conc_df <- varframe %>% filter(.bx0 == this_neighbor)
    
    # handle depths
    if(nrow(neighbor_conc_df)>nrow(df1)){ # cut to the deepest depth if the box we are filling is shallower than  its neighbor
      neighbor_conc_df <- neighbor_conc_df[1:nrow(df1),]
    } else if (nrow(neighbor_conc_df)<nrow(df1)) { # repeat the last row until the bottom if the box to fill is deeper than its neighbor
      diff_row <- nrow(df1)-nrow(neighbor_conc_df)
      neighbor_conc_df <- rbind(neighbor_conc_df, 
                                do.call('rbind', replicate(diff_row,neighbor_conc_df[nrow(neighbor_conc_df),],simplify = FALSE)))
    }
    
    df2 <- data.frame(df1 %>% select(.bx0,maxz), 
                      neighbor_conc_df %>% select(var_mean:unit)) 
    
    if(df2$maxz[1]==0){
      df2 <- df2 %>% mutate(var_mean=0) # if island box set concentration to 0
    }
    
    list_empty_boxes[[i]] <- df2
    
  }
  
  empty_box_reconstructed <- rbindlist(list_empty_boxes)
  
  
  # now bind the dataframes
  varframe <- rbind(varframe,empty_box_reconstructed)
  # arrange
  varframe <- varframe %>% arrange(.bx0,maxz)
  varframe
}

npz_vars_list_reconstructed <- lapply(npz_vars_list,fill_boxes)

# unlist
npz_vars <- rbindlist(npz_vars_list_reconstructed)

# drop Fe, we do not need it for now
npz_vars <- npz_vars %>% filter(long_name!='time-averaged iron concentration')

# print the units
print(npz_vars %>% select(long_name,unit) %>% distinct())

# map variables to Atlantis GOA groups

# Mapping NPZ to Atlantis groups

# Zooplankton in the NPZ model is organized in 5 groups: euphausiids, two copepod groups, large, and small microzooplankton
# Zoopl in Atlantis is organized in 3 grouptypes: large, medium, and small zoopl, and then euphausiids
# Mapping these two together is problematic. The two copepod groups have different vertical movements; the two microzoo groups are potentially cannibalistic
# Any grouping of the NPZ plankton is going to violate some restrictions, in addition there is no large zooplankton in NPZ and other sources are sparse
# Option 1: Eup as Eup; large Neocalanus as large zooplankton; small calanoids and large microzooplankton as meso; small micro as micro
# This causes cannibalism within the mesozooplankton (copepods eat microzoo)
# Option 2: Eup as Eup, large and small copepods as meso, large and small micro as micro. This means cannibalism in the micro, and
# probably a mix of vertical movements for meso (e.g. only some proportion moves). Also there is no large zooplankton, which will
# need to be approximated with something else. The advantage here is that I am not making too many assumptions of how to aggregate
# species on different tropchi levels (e.g. copepods and microzoo), and I am avoiding single-species groups. Approximate large zoopl
# distributions with Eup distributions and biomasses, aware that it will all shift pretty fast

# time-averaged small coastal copepod concentration <- Mesozooplankton
# time-averaged detritus concentration <- Detritus (split lab, ref, car as 40:40:20)
# time-averaged euphausiid concentration <- Euphausiids (copy it to macrozooplankton as well)
# time-averaged large microzooplankton concentration <- Microzooplankton
# time-averaged small microzooplankton concentration <- Microzooplankton
# time-averaged neocalanus spp. concentration <- Mesozooplankton
# time-averaged ammonia concentration <- NH3
# time-averaged nitrate concentration <- NO3
# time-averaged large phytoplankton concentration <- Diatoms
# time-averaged small phytoplankton concentration <- Picophytoplankton

npz_rows <- nrow(npz_vars)/length(unique(npz_vars$long_name))
npz_vars <- data.frame(npz_vars, 'Name'=c(rep('Mesozooplankton_N',npz_rows),
                                          rep('Detritus_N',npz_rows),
                                          rep('Euphausiids_N',npz_rows),
                                          rep('Microzooplankton_N',npz_rows),
                                          rep('Microzooplankton_N',npz_rows),
                                          rep('Mesozooplankton_N',npz_rows),
                                          rep('NH3_N',npz_rows),
                                          rep('NO3_N',npz_rows),
                                          rep('Diatoms_N',npz_rows),
                                          rep('Picophytoplankton_N',npz_rows)))


npz_atlantis <- npz_vars %>% 
  group_by(Name,.bx0,maxz,unit) %>%
  summarize(Value = sum(var_mean)) %>%
  ungroup() %>%
  arrange(Name,.bx0,desc(maxz)) %>%
  group_by(Name,unit,.bx0) %>%
  mutate(lyr = 0:(length(maxz)-1)) %>%
  ungroup() %>%
  complete(Name,.bx0,lyr) %>%
  select(-maxz,-unit)


# Calculate virgin biomass (t) for init calculations ------------------------------------------------

npz_biomass <- npz_vars %>% 
  group_by(Name,.bx0,maxz,unit) %>%
  summarize(Value = sum(var_mean)) %>%
  ungroup() %>%
  arrange(Name,.bx0,maxz) %>%
  group_by(Name,unit,.bx0) %>%
  mutate(minz=lag(maxz,default=0), # add depth layer thickness
         dz=maxz-minz) %>%
  ungroup() %>%
  left_join((atlantis_box %>% st_set_geometry(NULL) %>% select(.bx0,area))) %>%
  mutate(layer_volume_m3=dz*area)

# check the units
unique(npz_biomass$unit) # "milligram carbon meter-3"   "millimole nitrogen meter-3"

# for "milligram carbon meter-3" 
mgCtoTon <- 20/1e9 # assuming C is AFDW, which is 1/20 of WW (for consistency with other Atlantis calculations)
# for "millimole nitrogen meter-3"
mmolNtoTon <- 14.01*5.7*20/1e9 

npz_biomass1 <- npz_biomass %>%
  rowwise() %>%
  mutate(Biomass_t=ifelse(unit=="milligram carbon meter-3",Value*mgCtoTon*layer_volume_m3,Value*mmolNtoTon*layer_volume_m3))

# add it all up and write it out as .csv to cut-paste in the parameter spreadsheets
# for simplicity here we ignore the detritus and nutrient in the sediments - the sediment layer has 1 m thickness and is thus negligible here
npz_biomass1 %>%
  group_by(Name) %>%
  summarise(Biomass_goa_t=sum(Biomass_t)) %>%
  write.csv('npz_biomass_tot_goa.csv',row.names = F)

# this above should be 8vars*109boxes*6lyrs=5232 rows

############################### write to a text file to modify in a text editor and append to the cdf generated with make.init.nc()

# end result needs to be, for NO3 in GOA layers:
# NO3_N
# wc0,wc1,wc2,wc3,wc4,wc5,s1, # 6 layers and sed
# wc0,wc1,wc2,_,_,s1, # 3 layers and sed
# _,_,_,_,_,_,_, # what about island boxes with 0 wc layers? do they have sediment layers?
# ...
# where each row has 7 entries (6 wc and 1 sed layer), starting from bottom wc and going up and sed is the last
# for each variable, need as many rows as there are boxes in the model, 109 for GOA

# NO3 and NH3 are in millimoles meter-3. Everything else is in mg C meter-3. 
# Use Redfield Ratio - 5.7 of C:N for consistency with the rest of the calculations 

# Plankton
# convert mg C m-3 to mg N m-3 (/5.7)
plankton <- c('Euphausiids_N','Microzooplankton_N','Mesozooplankton_N','Diatoms_N','Picophytoplankton_N')

# assume that plankton concentrations in the sediment is 0

for (i in 1:length(plankton)){
  
  dat <- npz_atlantis %>%
    filter(Name==plankton[i]) %>%
    mutate(Value_N = signif(Value/5.7,digits = 5)) %>% # turn to mg N m-3
    replace_na(list(Value_N=0)) %>%
    select(Name,.bx0,Value_N) %>%
    nest(data = Value_N) %>%
    mutate(data = purrr::map(data,function(x) data.frame(matrix(c(t(x),0),nrow = 1)))) 
  
  write.table(rbindlist(dat$data),paste0('../outputs/init/',dat$Name[1],'.txt'), row.names = FALSE, col.names = FALSE, sep = ', ', eol = ',\n')
  
  # write out Silica for diatoms. Diatoms in the GOA should not be Si-limited (Hinckley et al. 2009), so we use a ratio of Si:N 3:1 for consistency with other Atlantis calculations and to err on the side of Si not being limiting
  if(plankton[i] == 'Diatoms_N'){
    write.table(rbindlist(dat$data)*3,paste0('../outputs/init/','Diatoms_S','.txt'), row.names = FALSE, col.names = FALSE, sep = ', ', eol = ',\n')
  }
  
}

# NO3 and NH3
# Convert millimol N m-3 to mg N m-3 (*14.01)
nuts <- c('NO3_N','NH3_N')

# assume that nutrient concentration in the sediment is same as the bottom layer at initial conditions

for (i in 1:length(nuts)){
  
  dat <- npz_atlantis %>%
    filter(Name==nuts[i]) %>%
    mutate(Value_N = signif(Value*14.01,digits = 5)) %>% # turn to mg N m-3
    replace_na(list(Value_N=0)) %>%
    select(Name,.bx0,Value_N) %>%
    nest(data = Value_N) %>%
    mutate(data = purrr::map(data,function(x) data.frame(matrix(c(t(x),t(x)[1]),nrow = 1)))) 
  
  write.table(rbindlist(dat$data),paste0('../outputs/init/',dat$Name[1],'.txt'), row.names = FALSE, col.names = FALSE, sep = ', ', eol = ',\n')
  
}

# Detritus
# Convert detritus concentrations to Detritus_labile, Detritus_refractory, and Carrion in mg N m-3 (*0.4,*0.4,*0.2)

# assume that detritus concentration in the sediment is same as the bottom layer at initial conditions

detritus <- npz_atlantis %>%
  filter(Name=='Detritus_N')

det_ref <- det_lab <- detritus %>% # ref and lab are the same to start
  mutate(Value_N = signif(Value*0.4/5.7,digits = 5)) %>%
  replace_na(list(Value_N=0)) %>%
  select(Name,.bx0,Value_N) %>%
  nest(data = Value_N) %>%
  mutate(data = purrr::map(data,function(x) data.frame(matrix(c(t(x),t(x)[1]),nrow = 1))))

carrion <-  detritus %>%
  mutate(Value_N = signif(Value*0.2/5.7,digits = 5)) %>%
  replace_na(list(Value_N=0)) %>%
  select(Name,.bx0,Value_N) %>%
  nest(data = Value_N) %>%
  mutate(data = purrr::map(data,function(x) data.frame(matrix(c(t(x),t(x)[1]),nrow = 1))))

write.table(rbindlist(det_ref$data),'../outputs/init/Detritus_refractory_N.txt', row.names = FALSE, col.names = FALSE, sep = ', ', eol = ',\n')
write.table(rbindlist(det_lab$data),'../outputs/init/Detritus_labile_N.txt', row.names = FALSE, col.names = FALSE, sep = ', ', eol = ',\n')
write.table(rbindlist(carrion$data),'../outputs/init/Carrion_N.txt', row.names = FALSE, col.names = FALSE, sep = ', ', eol = ',\n')

#####################

