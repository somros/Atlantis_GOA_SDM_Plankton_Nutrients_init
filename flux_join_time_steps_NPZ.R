# 1/14/2022
# Alberto Rovellini
# This code takes ROMS output (NEP 10K for now) and maps concentrations of plankton and nutrients to the Atlantis GOA boxes
# Two steps:
# 1. take the January climatology from ROMS files. Do this with nco (or cdo?) and run it from this code on the appropriate folder. This step results in a 1-ts netcdf
# 2. Take the 1-ts netcdf and apply the code to map data to Atlantis
# 3. Be careful with concentrations (which these all will be), the integration process will be different from state variables. For plankton we need a biomass value per cell (sum in a volume), whereas for nutrients it may be a concentration (mean per cell)

# Packages and settings
library(raster)
library(tidyverse)
library(sf)
library(here)
library(tidync)
library(rbgm)
library(viridis)
library(angstroms)
library(tabularaster)

select <- dplyr::select
here <- here::here
map <- purrr::map
options(dplyr.summarise.inform=FALSE)

# STEP 1: use NCO or CDO to go from ROMS data to obtain a Jan climatology from the ROMS-NPZ data
# If you want to do this you need to have NCO or CDO installed (http://nco.sourceforge.net/, https://code.mpimet.mpg.de/projects/cdo)
# For now I am doing this from Cygwin command line (I am on Windows) - there is a package to do this in R
# Here are the commands from Cygwin, after navigating to the folder with the ROMS-NPZ data of interest:

# I cannot do this in CDO because the first AND second steps of merging over time and then averaging drop 1-D variables (including Vtransform) from my ROMS data, see thread at https://code.mpimet.mpg.de/boards/1/topics/1177
# Leaving here the CDO commands I tried for future reference

# cdo mergetime nep5_avg_08*.nc one_file.nc    # NO THIS KILLS 1-d VARS
# cdo monavg one_file.nc monthly_means.nc  #gets monthly means for all vars 
# cdo selmon,1 monthly_means.nc jan.nc  #keeps jan only
# cdo selyear,2017 jan.nc jan2017.nc  #only need this last one because one file had a mislabeled time step (1900-01-01), keep Jan 2017 only

# NCO seems to work better, in that all variables are retained after concatenating and taking the average

# ncrcat nep5_avg_0804.nc nep5_avg_0805.nc nep5_avg_0806.nc nep5_avg_0807.nc nep5_avg_0808.nc one_file.nc # (not very practical for larger sets of files)
# ncra -d ocean_time,6,35 one_file.nc jan_mean.nc

# UPDATE 17/01/2022
# We can use an extension of these NCO commands to calculate seasonal climatologies follow these steps

# ncrcat nep5_avg_0804.nc nep5_avg_0805.nc nep5_avg_0806.nc nep5_avg_0807.nc nep5_avg_0808.nc nep5_avg_0809.nc nep5_avg_0810.nc nep5_avg_0811.nc nep5_avg_0812.nc nep5_avg_0813.nc /cygdrive/c/Users/Alberto\ Rovellini/Documents/GOA/SDM/Plankton_and_nutrients/data/ROMS_climatologies/s2.nc
# 
# ncra -d ocean_time,6,94 s1.nc s1_clim.nc 
# 
# ncrcat nep5_avg_0813.nc nep5_avg_0814.nc nep5_avg_0815.nc nep5_avg_0816.nc nep5_avg_0817.nc nep5_avg_0818.nc nep5_avg_0819.nc nep5_avg_0820.nc nep5_avg_0821.nc nep5_avg_0822.nc /cygdrive/c/Users/Alberto\ Rovellini/Documents/GOA/SDM/Plankton_and_nutrients/data/ROMS_climatologies/s2.nc
# 
# ncra -d ocean_time,6,96 s2.nc s2_clim.nc
# 
# ncrcat nep5_avg_0822.nc nep5_avg_0823.nc nep5_avg_0824.nc nep5_avg_0825.nc nep5_avg_0826.nc nep5_avg_0827.nc nep5_avg_0828.nc nep5_avg_0829.nc nep5_avg_0830.nc nep5_avg_0831.nc /cygdrive/c/Users/Alberto\ Rovellini/Documents/GOA/SDM/Plankton_and_nutrients/data/ROMS_climatologies/s3.nc
# 
# ncra -d ocean_time,7,98 s3.nc s3_clim.nc
# 
# ncrcat nep5_avg_0831.nc nep5_avg_0832.nc nep5_avg_0833.nc nep5_avg_0834.nc nep5_avg_0835.nc nep5_avg_0836.nc nep5_avg_0837.nc nep5_avg_0838.nc nep5_avg_0839.nc nep5_avg_0840.nc /cygdrive/c/Users/Alberto\ Rovellini/Documents/GOA/SDM/Plankton_and_nutrients/data/ROMS_climatologies/s4.nc
# 
# ncra -d ocean_time,9,100 s4.nc s4_clim.nc

# STEP 2: run code on the jan climatology

# read ROMS data
romsfile <- 'C:/Users/Alberto Rovellini/Documents/GOA/SDM/Plankton_and_nutrients/data/ROMS_climatologies/s4_clim.nc'
roms <- tidync(romsfile)
# read GOA ROMS grid
romsfile2 <- 'C:/Users/Alberto Rovellini/Documents/GOA/ROMS/data/roms/NEP_grid_5a.nc'
roms2 <- tidync(romsfile2)
# read Atlantis BGM
atlantis_bgm <- read_bgm('C:/Users/Alberto Rovellini/Documents/GOA/ROMS/data/atlantis/GOA_WGS84_V4_final.bgm')
#Atlantis geometry as an sf shapefile
atlantis_sf <- atlantis_bgm %>% box_sf()

# get list of ROMS variables
roms_vars <- hyper_grids(roms) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables associated with that grid and make a reference table
  purrr::map_df(function(x){
    roms %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })
# ... and grid variables
roms2_vars <- hyper_grids(roms2) %>% # all available grids in the ROMS ncdf
  pluck("grid") %>% # for each grid, pull out all the variables asssociated with that grid and make a reference table
  purrr::map_df(function(x){
    roms2 %>% activate(x) %>% hyper_vars() %>% 
      mutate(grd=x)
  })

# extract time
time_grd <- roms_vars %>% filter(name=="ocean_time") %>% pluck('grd')
roms_time <- roms %>% activate(time_grd) %>% hyper_tibble() %>% pull()

# ROMS grids
# Rho grid
# Horizontal coordinates: xi, eta, lat and lon of rho points
latlon_rhogrd <- roms2_vars %>% filter(name=="lat_rho") %>% pluck('grd')
roms_rho <- roms2 %>% activate(latlon_rhogrd) %>% hyper_tibble() %>%
  select(lon_rho,lat_rho,xi_rho,eta_rho) %>% 
  mutate(rhoidx=row_number()) # add index

# Vertical coordinates: h, s_rho, Cs_r, Cs_w
h <- roms2 %>% activate(latlon_rhogrd) %>% hyper_tibble() %>% select(xi_rho,eta_rho,h)
s_rho_grd <- roms_vars %>% filter(name=="s_rho") %>% pluck('grd')
s_rho <- roms %>% activate(s_rho_grd) %>% hyper_tibble()
## Cs_r is the S-coord stretching (length is num layers from -1 to 0 describing what portion of the w.c. each layer spans)
## We pull the Cs_r values from the roms ncdf
## one Cs_r value per s-coordinate
Cs_r <- s_rho %>% pluck('Cs_r')
s_w_grd <- roms_vars %>% filter(name=="s_w") %>% pluck('grd')
s_w <- roms %>% activate(s_w_grd) %>% hyper_tibble()
Cs_w <- s_w %>% pluck('Cs_w')

# u and v grids
# u and v grids from ROMS data
latlon_ugrd <-roms2_vars %>% filter(name=="lat_u") %>% pluck('grd')
latlon_vgrd <-roms2_vars %>% filter(name=="lat_v") %>% pluck('grd')

# pull the lon/lats
roms_u <- roms2 %>% activate(latlon_ugrd) %>% hyper_tibble() %>% select(lon_u,lat_u,xi_u,eta_u) %>% mutate(uidx=row_number())
roms_v <- roms2 %>% activate(latlon_vgrd) %>% hyper_tibble() %>% select(lon_v,lat_v,xi_v,eta_v) %>% mutate(vidx=row_number())

# append coordinates to the native lat lon from ROMS
append_xy_coords <- function(lonlatdat,xyproj=atlantis_bgm$extra$projection,lon_col="lon_rho",lat_col="lat_rho"){
  lonlatdat %>% 
    st_as_sf(coords=c(lon_col,lat_col),crs=4326,remove=F) %>%  # convert to spatial object
    st_transform(xyproj) %>%  # convert to Atlantis coords
    mutate(x = st_coordinates(.)[,1],
           y = st_coordinates(.)[,2]) # grab x and y coordinates and add them as attributes
}

rhoxy<- append_xy_coords(roms_rho,lon_col="lon_rho",lat_col="lat_rho") %>% mutate(rhoidx=row_number())
uxy <- append_xy_coords(roms_u,lon_col="lon_u",lat_col="lat_u")%>% mutate(uidx=row_number())
vxy <- append_xy_coords(roms_v,lon_col="lon_v",lat_col="lat_v")%>% mutate(vidx=row_number())

# join ROMS grids with Atlantis boxes (rho grid) and faces (u and v grids)
# Boxes
boxes_rho_join <- atlantis_sf %>% st_join(rhoxy)

# Faces
faces <- atlantis_bgm$faces %>% select(-label)

faces_sf <- atlantis_bgm %>% face_sf() %>% 
  mutate(label = 0:(length(label)-1)) %>% # creates a new index 'face_id' starting from 0 and increasing, as the 'label' column produced by rbgm::face_sf() is incorrect (tested in R 4.0.4)
  # join attribute data
  left_join(faces,by=c('label'='.fx0')) %>% 
  rename(.fx0=label)

# construct a buffer around each face
faces_buffer <- st_buffer(faces_sf,dist=10000) #TODO: make this an argument
# join u points
faces_u_join <- faces_buffer %>% st_join(uxy)
# join v points
faces_v_join <- faces_buffer %>% st_join(vxy)

# just a couple of checks
# insert a check for which boxes have no overlapping rho points
empty_boxes<- boxes_rho_join %>% 
  st_set_geometry(NULL) %>% 
  filter(is.na(rhoidx)) %>% 
  select(.bx0) %>% 
  distinct() %>% pull()
paste0("Atlantis boxes with no ROMS points are boxes ",paste(empty_boxes,collapse = ","))

# ... and one for which faces do not intercept u and v points
empty_faces<- faces_u_join %>% 
  st_set_geometry(NULL) %>% 
  filter(is.na(uidx)) %>% 
  select(.fx0) %>% 
  distinct() %>% pull()
paste0("Atlantis faces with no u points are faces ",paste(empty_faces,collapse = ","))

# get indeces of rho, u, and v points that overlap with Atlantis geometry, to subset large ROMS files and reduce memory chokes
min_xi_rho <- min(boxes_rho_join$xi_rho, na.rm = TRUE)
max_xi_rho <- max(boxes_rho_join$xi_rho, na.rm = TRUE)
min_eta_rho <- min(boxes_rho_join$eta_rho, na.rm = TRUE)
max_eta_rho <- max(boxes_rho_join$eta_rho, na.rm = TRUE)

min_xi_u <- min(faces_u_join$xi_u, na.rm = TRUE)
max_xi_u <- max(faces_u_join$xi_u, na.rm = TRUE)
min_eta_u <- min(faces_u_join$eta_u, na.rm = TRUE)
max_eta_u <- max(faces_u_join$eta_u, na.rm = TRUE)

min_xi_v <- min(faces_v_join$xi_v, na.rm = TRUE)
max_xi_v <- max(faces_v_join$xi_v, na.rm = TRUE)
min_eta_v <- min(faces_v_join$eta_v, na.rm = TRUE)
max_eta_v <- max(faces_v_join$eta_v, na.rm = TRUE)

# Atlantis depth
# enter the depth breaks in the model
atlantis_z <- c(-30,-100,-200,-500,-1000,-4000) #TODO: make this an argument

# small function to build enough depth layers for each box, starting with the defined layers above
# in this function, botz is the bottom depth given for each box, which is available in the .bgm file 
build_Atlantis_depths <- function(botz,lyrs){
  # bottom of each layer, starting from shallowest
  lyrbot<-lyrs
  # layers to use are all those that are shallower than the given botz
  lyr_vec <- lyrbot[lyrbot>botz]
  # the depth of the deepest layer is equal to botz
  lyr_vec <- c(lyr_vec,botz)
  # in Atlantis, each box has the same number of depth layers, but some layers have zero thickness
  # so we have to pad with zeroes to make all boxes have the same number of layers
  nzeroes <- length(lyrs)-length(lyr_vec)
  lyr_vec <- c(lyr_vec,rep(0,nzeroes))
  return(lyr_vec)
}
# construct the depth profiles of each Atlantis box

atlantis_depths <- atlantis_bgm$boxes %>% select(.bx0,botz) %>% 
  # apply the function above to create the layers
  mutate(maxz=purrr::map(botz,~build_Atlantis_depths(.,lyrs=atlantis_z))) %>% 
  unnest(cols=c(maxz)) %>% 
  # add a minimum depth for each layer
  group_by(.bx0) %>% 
  mutate(minz=lag(maxz,1,default = 0),atlantis_layer=1:length(atlantis_z)) %>% 
  # add a layer thickness calculation
  mutate(dz=minz-maxz) %>% 
  # "dummy" layers (layers too deep for a given Atlantis box) should have minz and dz=0
  mutate(minz=ifelse(maxz==0,0,minz),dz=ifelse(maxz==0,0,dz)) %>% 
  ungroup() %>% 
  select(.bx0,atlantis_layer,minz,maxz,dz)

# face depths (for later flux calcs)
face_depths <- faces %>% 
  left_join(atlantis_depths %>% select(.bx0,atlantis_layer,dz),by=c("left"=".bx0")) %>% 
  mutate(left_area=length*dz) %>% 
  rename(dz_left=dz) %>% 
  left_join(atlantis_depths %>% select(.bx0,atlantis_layer,dz),by=c("right"=".bx0","atlantis_layer")) %>% 
  mutate(right_area=length*dz) %>% 
  rename(dz_right=dz) %>% 
  # area of the face is the smallest area, maintaining NAs (if one box is deeper than its neighbor, no flux)
  rowwise() %>%
  mutate(dz_max = max(dz_left,dz_right),
         face_area=ifelse((left_area>0 & right_area >0), min(left_area, right_area), NA)) 

# For GOA: custom version of romshcoords()
# TODO: put this in a separate file when packaging

ncget <- function(x, varname) {
  nc <- ncdf4::nc_open(x)
  on.exit(ncdf4::nc_close(nc))
  ncdf4::ncvar_get(nc, varname)
}

set_indextent <- function(x) {
  setExtent(x, extent(0, ncol(x), 0, nrow(x)))
}

romshcoords_goa <- function(x, y, grid_type = "rho", slice, ..., S = "Cs_r", depth = "h", simple = FALSE){
  h <- romsdata(x, varname = depth)
  Cs_r <- ncget(y, S)
  v <- values(h)
  if (simple) {
    ## simplistic, early version - probably should be defunct
    out <- set_indextent(brick(array(rep(rev(Cs_r), each = length(v)) * v, 
                                     c(ncol(h), nrow(h), length(Cs_r))), transpose = TRUE))
  } else {
    grid_type <- match.arg(tolower(grid_type),c("rho","psi","u","v","w"))
    
    Vtransform <- as.integer(ncget(y,"Vtransform"))
    if (!Vtransform %in% c(1,2)) stop("Vtransform must be 1 or 2")
    
    hc <- ncget(y,"hc")
    
    depth_grid <- if (grid_type=="w") "w" else "rho"
    
    zeta <- if (missing(slice)) 0 else stop("slice not supported yet")##angstroms::romsdata2d(x,"zeta",slice=slice,transpose=FALSE)
    N <- length(ncget(y,"Cs_r"))
    Np <- N+1
    
    h <- ncget(x,"h")
    hmin <- min(h)
    hmax <- max(h)
    
    Lp <- dim(h)[1]
    Mp <- dim(h)[2]
    L <- Lp-1
    M <- Mp-1
    
    z <- array(NA,dim=c(Lp,Mp,if (grid_type=="w") Np else N))
    
    ## Compute vertical stretching function, C(k):
    ##stretch <- stretching(x,depth_grid)
    if (depth_grid=="w") {
      stretch <- list(C=ncget(y,"Cs_w"),s=ncget(y,"s_w"))
    } else {
      stretch <- list(C=ncget(y,"Cs_r"),s=ncget(y,"s_rho"))
    }
    
    ## Average bathymetry and free-surface at requested C-grid type.
    if (grid_type=="rho") {
      hr <- h
      zetar <- zeta
    } else if (grid_type=="psi") {
      hp <- 0.25*(h[1:L,1:M]+h[2:Lp,1:M]+h[1:L,2:Mp]+h[2:Lp,2:Mp])
      zetap <- 0.25*(zeta[1:L,1:M]+zeta[2:Lp,1:M]+zeta[1:L,2:Mp]+zeta[2:Lp,2:Mp])
    } else if (grid_type=="u") {
      hu <- 0.5*(h[1:L,1:Mp]+h[2:Lp,1:Mp])
      zetau <- 0.5*(zeta[1:L,1:Mp]+zeta[2:Lp,1:Mp])
    } else if (grid_type=="v") {
      hv <- 0.5*(h[1:Lp,1:M]+h[1:Lp,2:Mp])
      zetav <- 0.5*(zeta[1:Lp,1:M]+zeta[1:Lp,2:Mp])
    } else if (grid_type=="w") {
      hr <- h
      zetar <- zeta
    } else {
      stop("unsupported grid_type: ",grid_type)
    }
    
    ## Compute depths (m) at requested C-grid location.
    
    if (Vtransform == 1) {
      if (grid_type=="rho") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hr
          z[,,k] <- z0 + zetar*(1.0 + z0/hr)
        }
      } else if (grid_type=="psi") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hp
          z[,,k] <- z0 + zetap*(1.0 + z0/hp)
        }
      } else if (grid_type=="u") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hu
          z[,,k] <- z0 + zetau*(1.0 + z0/hu)
        }
      } else if (grid_type=="v") {
        for (k in seq_len(N)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hv;
          z[,,k] <- z0 + zetav*(1.0 + z0/hv)
        }
      } else if (grid_type=="w") {
        z[,,1] <- -hr
        for (k in seq(from=2,to=Np,by=1)) {
          z0 <- (stretch$s[k]-stretch$C[k])*hc + stretch$C[k]*hr
          z[,,k] <- z0 + zetar*(1.0 + z0/hr)
        }
      } else {
        stop("unsupported grid_type: ",grid_type)
      }
    } else if (Vtransform == 2) {
      if (grid_type=="rho") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hr)/(hc+hr)
          z[,,k] <- zetar+(zeta+hr)*z0
        }
      } else if (grid_type=="psi") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hp)/(hc+hp)
          z[,,k] <- zetap+(zetap+hp)*z0
        }
      } else if (grid_type=="u") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hu)/(hc+hu)
          z[,,k] <- zetau+(zetau+hu)*z0
        }
      } else if (grid_type=="v") {
        for (k in seq_len(N)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hv)/(hc+hv)
          z[,,k] <- zetav+(zetav+hv)*z0
        }
      } else if (grid_type=="w") {
        for (k in seq_len(Np)) {
          z0 <- (hc*stretch$s[k]+stretch$C[k]*hr)/(hc+hr)
          z[,,k] <- zetar+(zetar+hr)*z0
        }
      } else {
        stop("unsupported grid_type: ",grid_type)
      }
    } else {
      stop("Vtransform must be 1 or 2")
    }
    ## FIXME all these flips and twirls can be applied more efficiently (or avoided)
    ## though should layers start at the surface and go down or ...
    
    out <- raster::flip(set_indextent(raster::brick(z, transpose = TRUE)), "y")
    ## NO - we want to start at the bottom, so we match romsdata3d
    #out <- raster::subset(out, rev(seq_len(raster::nlayers(out))))
    
  } 
  
  out
}

# convert ROMS s-coordinates to depth with Mike Sumner's angstroms package

romsdepths <- romshcoords_goa(x = romsfile2, y = romsfile, S = "Cs_r", depth = "h")

# using tabularaster to convert to tibble
# and a indexing template with "by_column" filling
romsi <- crossing(xi_rho=1:dim(romsdepths)[2],eta_rho=1:dim(romsdepths)[1]) %>% arrange(-eta_rho) %>% mutate(cellindex=row_number()) # making sure that the join by cellindex below is correct - doing this for consistency with the way tabularaster::as_tibble() unpacks the raster cells 
romsdepthsdf <- tabularaster::as_tibble(romsdepths,dim=F) %>% 
  arrange(cellindex) %>% 
  left_join(romsi,by='cellindex') %>% 
  set_names(c("romsdepth","cellindex","xi_rho","eta_rho")) %>% 
  group_by(cellindex,xi_rho,eta_rho) %>% 
  nest(romsdepth=c(romsdepth)) %>% ungroup() %>% 
  mutate(romsdepth=purrr::map(romsdepth,function(x)x[['romsdepth']])) %>%
  filter(between(xi_rho, min_xi_rho, max_xi_rho) & between(eta_rho, min_eta_rho, max_eta_rho))

# do the same for depth at w
romsdepths_w <- romshcoords_goa(x = romsfile2, y = romsfile, grid_type = "w", S = "Cs_w", depth = "h")
#romsdepths_w <- subset(romsdepths_w, dim(romsdepths_w)[3]:1) # angstroms returns depths from shallowest to deepest, while below we extract ROMS variables from deepest to shallowest. flipping here or else it maps ROMS variables upside-down

# using tabularaster to convert to tibble
romsdepthsdf_w <- tabularaster::as_tibble(romsdepths_w,dim=F) %>% 
  arrange(cellindex) %>% 
  left_join(romsi,by='cellindex') %>% 
  set_names(c("romsdepth","cellindex","xi_rho","eta_rho")) %>% 
  group_by(cellindex,xi_rho,eta_rho) %>% 
  nest(romsdepth=c(romsdepth)) %>% ungroup() %>% 
  mutate(romsdepth=purrr::map(romsdepth,function(x)x[['romsdepth']])) %>%
  filter(between(xi_rho, min_xi_rho, max_xi_rho) & between(eta_rho, min_eta_rho, max_eta_rho))

# Matching keys for boxes and faces
# TODO: shift this above - this will not need repeating at each time step
# build a matching key for Atlantis boxes...
boxes_rho_thin <- boxes_rho_join %>% 
  st_set_geometry(NULL) %>% 
  dplyr::select(.bx0,xi_rho,eta_rho,rhoidx,area) %>% 
  drop_na()
boxes_rho_join_with_depth <- atlantis_depths %>% 
  left_join(boxes_rho_thin,by=c(".bx0")) %>% 
  ungroup() %>% 
  drop_na()

# make vectors of rho point indices
boxes_depths_rhopts <- boxes_rho_join_with_depth %>% 
  select(-dz) %>% 
  group_by(.bx0,atlantis_layer,minz,maxz) %>% 
  nest(rhopts=c(xi_rho,eta_rho,rhoidx)) %>% 
  mutate(rhovec=purrr::map(rhopts,~pluck(.,'rhoidx')))%>%
  ungroup()

# ... and for Atlantis faces
faces_u_thin <- faces_u_join %>%
  st_set_geometry(NULL) %>%
  select(.fx0,xi_u,eta_u,uidx) %>%
  drop_na()
faces_u_join_with_depth <- face_depths %>%
  left_join(faces_u_thin, by = c(".fx0")) %>%
  ungroup() %>%
  #drop_na() %>% # turn this on or off depending on whether we want to have NA fluxes in non-existing layers or not, cannot recall what HydroCOnstruct wants
  select(.fx0,atlantis_layer,dz_max,xi_u,eta_u,uidx,face_area) %>%
  group_by(uidx,.fx0) %>%
  mutate(maxz=-cumsum(dz_max),minz=-lag(-maxz,default=0)) %>%
  ungroup()%>%
  select(-dz_max)

# again for v
faces_v_thin <- faces_v_join %>%
  st_set_geometry(NULL) %>%
  select(.fx0,xi_v,eta_v,vidx) %>%
  drop_na()
faces_v_join_with_depth <- face_depths %>%
  left_join(faces_v_thin, by = c(".fx0")) %>%
  ungroup() %>%
  #drop_na() %>%
  select(.fx0,atlantis_layer,dz_max,xi_v,eta_v,vidx,face_area) %>%
  group_by(vidx,.fx0) %>%
  mutate(maxz=-cumsum(dz_max),minz=-lag(-maxz,default=0)) %>%
  ungroup()%>%
  select(-dz_max)

# make vectors of u point indices
faces_depths_upts <- faces_u_join_with_depth %>% 
  group_by(.fx0,atlantis_layer,minz,maxz) %>% 
  nest(upts=c(xi_u,eta_u,uidx)) %>% 
  mutate(uvec=map(upts,~pluck(.,'uidx'))) %>%
  ungroup()
# and v
faces_depths_vpts <- faces_v_join_with_depth %>% 
  group_by(.fx0,atlantis_layer,minz,maxz) %>% 
  nest(vpts=c(xi_v,eta_v,vidx)) %>% 
  mutate(vvec=map(vpts,~pluck(.,'vidx'))) %>%
  ungroup()
# now join
faces_depths_uvpts <- faces_depths_upts %>% full_join(faces_depths_vpts, by = c('.fx0','atlantis_layer','maxz','minz','face_area'))

# pull ROMS static variables (salt, temp, w, etc.) from Vertical interpolation
# function to wrap 'spline' and apply it to ROMS data and depth vector (1 m intervals)
interp_foo <- function(romsdepths,romsvar) {
  depths_out <- seq(round(min(romsdepths)),0,by=1) # 1m interpolation, starting from deepest
  interp <- spline(romsdepths,romsvar,xout=depths_out) %>% pluck('y')
  return(tibble(depth=depths_out,val=interp))
}

# function that pulls variables from ROMS at rho points, and interpolates the vertical values at each 1 m
interpolate_var <- function(variable, time_step){
  grd <- roms_vars %>% filter(name==variable) %>% pluck('grd')
  # pull the env data
  # interpolate the env data
  # do this step conditional to join with the appropriate depth data frame depending on the variable
  # if variable is horizontal velocity
  # if(variable == "u") {
  #   dat <- roms %>% activate(grd) %>%
  #     hyper_tibble(select_var=variable, 
  #                  xi_u = between(xi_u, min_xi_u, max_xi_u), 
  #                  eta_u = between(eta_u, min_eta_u, max_eta_u),
  #                  ocean_time = ocean_time == time_step)
  #   
  #   interp_dat <- dat %>% 
  #     dplyr::select(xi_u,eta_u,!!variable) %>% 
  #     nest(data=c(!!variable))%>% 
  #     mutate(evar=purrr::map(data,~.[[1]])) %>% 
  #     inner_join(romsdepthsdf,by=c('xi_u'='xi_rho','eta_u'='eta_rho')) %>% 
  #     inner_join(roms_u,by=c('xi_u','eta_u')) 
  # } else if (variable == "v") {
  #   dat <- roms %>% activate(grd) %>%
  #     hyper_tibble(select_var=variable, 
  #                  xi_v = between(xi_v, min_xi_v, max_xi_v), 
  #                  eta_v = between(eta_v, min_eta_v, max_eta_v),
  #                  ocean_time = ocean_time == time_step)
  #   
  #   interp_dat <- dat %>% 
  #     dplyr::select(xi_v,eta_v,!!variable) %>% 
  #     nest(data=c(!!variable))%>% 
  #     mutate(evar=purrr::map(data,~.[[1]])) %>% 
  #     inner_join(romsdepthsdf,by=c('xi_v'='xi_rho','eta_v'='eta_rho')) %>% 
  #     inner_join(roms_v,by=c('xi_v','eta_v')) 
  # } else { # if variable is a state variable or vertical velocity
    dat <- roms %>% activate(grd) %>%
      hyper_tibble(select_var=variable, 
                   xi_rho = between(xi_rho, min_xi_rho, max_xi_rho), 
                   eta_rho = between(eta_rho, min_eta_rho, max_eta_rho),
                   ocean_time = ocean_time == time_step)
    
    interp_dat <- dat %>% 
      dplyr::select(xi_rho,eta_rho,!!variable) %>% 
      nest(data=c(!!variable))%>% 
      mutate(evar=purrr::map(data,~.[[1]])) %>%
      inner_join(romsdepthsdf,by=c('xi_rho','eta_rho')) %>% 
      inner_join(roms_rho,by=c('xi_rho','eta_rho')) 
  #   if (variable == 'w') { # if the variable is w we need a different vertical mapping
  #     interp_dat <- interp_dat %>%
  #       inner_join(romsdepthsdf_w,by=c('xi_rho','eta_rho')) %>% 
  #       inner_join(roms_rho,by=c('xi_rho','eta_rho'))
  #     # drop NAs - there are a lot in w - might want to check why
  #     # idx <- unlist(lapply(interp_dat$evar, function(x)length(which(is.na(x)))/length(x)))
  #     # interp_dat <- interp_dat[-which(idx==1),]
  #   } else { # for all other state variables
  #     interp_dat <- interp_dat %>%
  #       inner_join(romsdepthsdf,by=c('xi_rho','eta_rho')) %>% 
  #       inner_join(roms_rho,by=c('xi_rho','eta_rho')) 
  #   }
  # }
  interp_dat <- interp_dat %>% 
    mutate(interp = purrr::map2(romsdepth,evar,interp_foo)) %>% 
    dplyr::select(-data,-evar,-romsdepth)
  return(interp_dat)
}

# From here down: repeat at each time step
# TODO: put this in a separate file / code block. 
# Also the functions that we call in here will need to go elsewhere and only be called once

init_vars <- c('NO3','NH4','PhS','PhL','MZS','MZL','Cop','NCa','Eup','Det','Iron')
this_time <- roms_time[1] # do not actually need time here
  
# apply interpolate_var
# TODO: operationalize this so that it works with whatever state variables people need to pull
# salt_interp <- interpolate_var('salt', time_step)
# temperature_interp <- interpolate_var('temp', time_step) 
# w_interp <- interpolate_var('w', time_step) 

plankton <- function(this_variable){
  
  var_interp <- interpolate_var(this_variable, this_time)
  
  this_unit <- ncmeta::nc_atts(romsfile,this_variable) %>% 
    filter(name %in% c('long_name','units')) %>% 
    tidyr::unnest(cols = c(value)) %>% 
    pull(value)
  
  # match interpolated data to Atlantis boxes, faces, and layers
  match_interpolated_data <- function(interp_dat,idx_vec,minz,maxz){
    #idx <- names(interp_dat[grep('idx',names(interp_dat))])
    interp_var<- interp_dat
    if('uidx'%in%names(interp_var)){
      interp_var<-interp_var%>%filter(uidx%in%idx_vec)
    }else if ('vidx'%in%names(interp_var)){
      interp_var<-interp_var%>%filter(vidx%in%idx_vec)
    }else{
      interp_var<-interp_var%>%filter(rhoidx%in%idx_vec)
    }
    interp_var<-interp_var %>%    # filter depths
      mutate(interp_within_depths=purrr::map(interp, ifelse(maxz==0,function(df) df %>% mutate(val = NA), function(df) df %>% filter(between(depth,maxz,minz))))) %>% # write out NAs if the thickness of the depth layer is 0 m, i.e. for non-existing Atlantis layers
      # unpack the interpolated data and calculate a mean across relevant rho points and depths
      pluck('interp_within_depths') %>% 
      bind_rows() %>% 
      pluck('val') %>% 
      mean(na.rm=T)
    
    return(interp_var)
  }
  
  # state variables, but not w
  boxes_statevars_interp <- boxes_depths_rhopts %>% 
    # apply the function to the interpolated salinity and temperature data
    mutate(var_mean=purrr::pmap_dbl(list(idx_vec=rhovec,minz=minz,maxz=maxz),match_interpolated_data,interp_dat=var_interp))
  
  var_out <- boxes_statevars_interp %>%
    drop_na() %>%
    #mutate(time_step = time_step, maxz = -maxz) %>% # how to turn this into iteration over time stesp? Either ts from NetCDF, or just iteration # of the function
    mutate(maxz = -maxz) %>%
    #ungroup() %>%
    #select(time_step,.bx0,maxz,var_mean) %>%
    select(.bx0,maxz,var_mean) %>%
    mutate(long_name = this_unit[1], unit = this_unit[2])
  
  write.csv(var_out, paste0('../outputs/s1s4/s4_',this_variable,'.csv'), row.names = FALSE)
}

purrr::map(init_vars,plankton)
