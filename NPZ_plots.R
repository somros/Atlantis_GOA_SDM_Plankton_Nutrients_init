---
title: "Plots for NPZ variables"
author: "Alberto Rovellini"
date: "8/26/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
# read NPZ variables, visualise them as plot

library(rbgm)
library(tidyverse)
library(viridis)
library(maps)
library(mapdata)
```

```{r}
variable_files <- list.files('outputs/', pattern = '.csv', full.names = TRUE)

test <- read.csv(variable_files[1])

test <- test %>% mutate(date = as.POSIXct(time_step,origin='1900-01-01',tz='UTC'))

# read in Atlantis boxes

atlantis_bgm <- read_bgm('data/atlantis/GOA_WGS84_V4_final.bgm')
atlantis_box <- atlantis_bgm %>% box_sf()
atlantis_crs <- atlantis_bgm$extra$projection
atlantis_bbox <- st_bbox(atlantis_sf)
```

```{r}
# join
var_join <- atlantis_box %>% left_join(test, by='.bx0')

# map depths to layers (where 1 is the surface and it increases with depth)
z <- c(0,30,100,200,500,1000,4000)
var_join <- var_join %>% mutate(atlantis_layer = findInterval(maxz,z,left.open = TRUE))

# TODO: pad this with zeroes - or think if it should be done in the main code instead.
```

View. Get the coastline.
```{r}
# view. get coastline
coast <- maps::map('worldHires', regions = c('USA','Canada'), plot = FALSE, fill = TRUE)
coast_sf <- coast %>% st_as_sf() %>% st_transform(crs=atlantis_crs)

var_join %>% slice_min(atlantis_layer) %>%
  ggplot()+
  geom_sf(aes(fill=var_mean))+
  scale_fill_viridis()+
  geom_sf(data = coast_sf)+
  coord_sf(xlim=c(atlantis_bbox$xmin,atlantis_bbox$xmax),ylim=c(atlantis_bbox$ymin,atlantis_bbox$ymax))+
  theme_minimal()+
  labs(fill = var_join$unit, title = var_join$long_name)
```

```{r, fig.width=10, fig.height=18}
# view vertical distribution

var_join %>% ggplot()+
  geom_point(aes(x=var_mean,y=atlantis_layer))+
  geom_line(aes(x=var_mean,y=atlantis_layer))+
  scale_y_reverse()+
  theme_minimal()+
  facet_wrap(~.bx0, ncol = 6)
```

