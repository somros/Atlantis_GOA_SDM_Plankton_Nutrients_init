# this code takes csv files with concentrations by layer by box and transforms into biomasses and then in s1-s4

library(rbgm)
library(tidyverse)
library(sf)
library(viridis)
library(maps)
library(mapdata)
library(data.table)

# read in Atlantis boxes

atlantis_bgm <- read_bgm('../data/GOA_WGS84_V4_final.bgm')
atlantis_box <- atlantis_bgm %>% box_sf()
atlantis_crs <- atlantis_bgm$extra$projection
atlantis_bbox <- st_bbox(atlantis_box)

# data 
variable_files <- list.files('../outputs/s1s4/', pattern = '.csv', full.names = TRUE)
variables <- unique(gsub(".*[_]([^.]+)[.].*", "\\1", variable_files))

key <- data.frame('npz'=variables,
                  'fg'=c('Mesozooplankton',
                  'Detritus',
                  'Euphausiids',
                  'Iron',
                  'Microzooplankton',
                  'Microzooplankton',
                  'Mesozooplankton',
                  'NH3',
                  'NO3',
                  'Diatoms',
                  'Picophytoplankton'))

# make an empty list
biomass_list <- list()

for(v in 1:length(variables)) {
  
  this_variable <- variables[v]
  
  these_files <- variable_files[grep(this_variable,variable_files)]
  
  # read these four
  jfm <- read.csv(these_files[1])
  amj <- read.csv(these_files[2])
  jas <- read.csv(these_files[3])
  ond <- read.csv(these_files[4])
  
  # merge them into one dataframe with concentration climatologies for all 4 quarters
  all_seasons <- Reduce(function(x, y) merge(x, y, by = c('.bx0','maxz','long_name','unit')), list(jfm,amj,jas,ond))
  
  # rename 
  all_seasons <- all_seasons %>% 
    set_names(c('.bx0','maxz','long_name','unit','conc_1','conc_2','conc_3','conc_4')) %>%
    arrange(.bx0,maxz)
  
  #######################################################################################
  # need to set conc in island boxes as 0, and conc in empty boxes as the nearest conc. This takes a bit of reconstruction
  missing_boxes <- setdiff(0:108, all_seasons %>% pull(.bx0) %>% unique())
  island_boxes <- atlantis_box %>% filter(.bx0 %in% missing_boxes, botz == 0) %>% pull(.bx0) %>% unique()
  no_data_boxes <- atlantis_box %>% filter(.bx0 %in% missing_boxes, botz != 0) %>% pull(.bx0) %>% unique()
  
  # reconstruct depth layers
  z <- c(30,100,200,500,1000,4000)
  
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
    
    df1 <- df1 %>% 
      distinct() %>%
      mutate(long_name = all_seasons$long_name[1],
             unit = all_seasons$unit[1])
    
    this_neighbor <- neighbor %>% filter(.bx0==tb) %>% pull(nearest)
    
    neighbor_conc_df <- all_seasons %>% filter(.bx0 == this_neighbor)
    
    if(nrow(neighbor_conc_df)>nrow(df1)){ # cut to the deepest depth
      neighbor_conc_df <- neighbor_conc_df[1:nrow(df1),]
    } else if (nrow(neighbor_conc_df)<nrow(df1)) { # repeat the last row until the bottom
      diff_row <- nrow(df1)-nrow(neighbor_conc_df)
      neighbor_conc_df <- rbind(neighbor_conc_df, 
                                do.call('rbind', replicate(diff_row,neighbor_conc_df[nrow(neighbor_conc_df),],simplify = FALSE)))
    }
    
    df2 <- data.frame(df1, neighbor_conc_df %>% select(conc_1:conc_4)) 
    
    if(df2$botz[1]==0){
      df2 <- df2 %>% mutate(across(starts_with("conc"), ~.*0))
    }
    
    df2 <- df2 %>% select(-botz)
    
    list_empty_boxes[[i]] <- df2
    
  }
  
  empty_box_reconstructed <- rbindlist(list_empty_boxes)
  
  
  # now bind the dataframes
  all_seasons <- rbind(all_seasons,empty_box_reconstructed)
  # arrange
  all_seasons <- all_seasons %>% arrange(.bx0,maxz)
  
  # add dz, calculate volume for each cell
  all_seasons <- all_seasons %>% 
    full_join(atlantis_box, by='.bx0') %>%
    group_by(box_id) %>%
    mutate(minz = lag(maxz, default = 0)) %>%
    ungroup() %>%
    mutate(dz = maxz-minz,
           vol_m3 = area*dz)
  
  # turn the concentrations in boundary boxes to 0
  all_seasons <- all_seasons %>%
    rowwise() %>%
    mutate(conc_1 = ifelse(isTRUE(boundary),0,conc_1),
           conc_2 = ifelse(isTRUE(boundary),0,conc_2),
           conc_3 = ifelse(isTRUE(boundary),0,conc_3),
           conc_4 = ifelse(isTRUE(boundary),0,conc_4))
  
  # calculate weight (biomass or nutrient weight) per layer per box from concentrations
  biomass_layers <- all_seasons %>% 
    mutate(biomass_s1 = conc_1*vol_m3,
           biomass_s2 = conc_2*vol_m3,
           biomass_s3 = conc_3*vol_m3,
           biomass_s4 = conc_4*vol_m3) %>%
    select(box_id,maxz,biomass_s1:biomass_s4) 
  
  
  # sum over layers for the total biomass per box, then calclate proportions
  season_biomass <- biomass_layers %>%
    group_by(box_id) %>%
    summarise(across(starts_with("biomass"), ~ sum(.x, na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(fg = key[v,2])
  
  biomass_list[[v]] <- season_biomass
  
}

biomass_by_group <- rbindlist(biomass_list) %>%
  pivot_longer(cols = biomass_s1:biomass_s4, names_to = 'season', values_to = 'biomass') %>%
  group_by(fg, season, box_id) %>%
  summarise(biomass = sum(biomass,na.rm = T)) %>%
  ungroup() %>%
  group_by(fg, season) %>%
  mutate(prop = biomass/sum(biomass)) %>%
  ungroup()

# reshape this so it is wide again
s1_s4_all <- biomass_by_group %>%
  select(-biomass) %>%
  pivot_wider(names_from = season, values_from = prop) %>%
  set_names(c('fg','box_id','s1','s2','s3','s4'))

# loop over FGs, write out a file and make a plot

FGs <- key$fg

for(i in 1:length(FGs)){
  
  this_fg <- FGs[i]
  this_s <- s1_s4_all %>% filter(fg==this_fg)
  
  # write out file
  write.csv(this_s, paste0('../outputs/s1s4/final_for_parameters/',this_fg,'.csv'), row.names = FALSE)
  
  # make a plot, with s1-s4
  coast <- maps::map('worldHires', regions = c('USA','Canada'), plot = FALSE, fill = TRUE)
  coast_sf <- coast %>% st_as_sf() %>% st_transform(crs=atlantis_crs)
  
  p <- atlantis_box %>%
    left_join(this_s %>% pivot_longer(cols = s1:s4, names_to = 'Quarter', values_to = 'Proportion'),
              by = 'box_id') %>%
    ggplot()+
    geom_sf(aes(fill=Proportion))+
    scale_fill_viridis()+
    geom_sf(data = coast_sf)+
    coord_sf(xlim=c(atlantis_bbox$xmin,atlantis_bbox$xmax),ylim=c(atlantis_bbox$ymin,atlantis_bbox$ymax))+
    theme_bw()+
    labs(title=paste0('S1-S4 for ', this_fg))+
    facet_wrap(~Quarter)
  
  ggsave(paste0('../outputs/s1s4/final_for_parameters/',this_fg,'.png'), p, width = 12, height = 8)
  
}
