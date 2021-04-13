library(measurements)
library(gganimate)
library(gifski)
library(png)
library(rgdal)
library(raster)
library(ggplot2)
library(rgeos)
library(mapview)
library(broom)
library(dplyr)
library(sp)
library(sf)
library(magick)
library(ggrepel)
#### Upload dataset containing weekly case info and lat long 
setwd("/Users/Harley/Documents/DRC Research")

#Import the gadm DRC map data, other necessary data
drc_shp_2 <- st_read("rdc_zone_de_sante_09092019")
weekly_cases <- read.csv("weekly_cases.csv")
names(weekly_cases) <- c("ID", "Date", "weekly_confirm", "SCOnset")

### Number Weeks
Date <- unique(weekly_cases$Date)
Weeks <- data.frame(Date)
Weeks$Week <- seq.int(nrow(Weeks))
weeks_comb <- left_join(weekly_cases, Weeks, by=c("Date"))
weeks_comb <- weeks_comb %>%
  group_by(SCOnset) %>%
  mutate(cumulative = cumsum(weekly_confirm)) %>%
  filter(!SCOnset %in% c("mwenga"))
weeks_comb$SCOnset <- gsub("mangurujipa", "manguredjipa", weeks_comb$SCOnset)

## Get the cumulative number of cases, aggregated weekly
weeks_comb$total <- cumsum(weeks_comb$weekly_confirm)

## Change the information on coordinates to the correct latitude longitude form for mapping
drc_shp_2$Nom <- tolower(drc_shp_2$Nom)
drc_joined <- left_join(weeks_comb, drc_shp_2[,c("Area", "Nom", "geometry", "coordx", "coordy")],
                        by=c("SCOnset"="Nom"))

drc_joined <- st_as_sf(drc_joined)
drc_joined$coordx <- gsub("°", "d", drc_joined$coordx)
drc_joined$coordx <- gsub(",", ".", drc_joined$coordx)
drc_joined$coordy <- gsub("°", "d", drc_joined$coordy)
drc_joined$coordy <- gsub(",", ".", drc_joined$coordy)

drc_joined$long <- as.numeric(char2dms(drc_joined$coordx, chd="d", 
                                             chm = "'", chs = '\"'))
drc_joined$lat <- as.numeric(char2dms(drc_joined$coordy, chd="d", 
                                             chm = "'", chs = '\"'))


## Pull map data from GADM
drc2 <- getData("GADM", country="Democratic Republic of the Congo", level=2)

## Create the ggplot which will be animated, details include zooming into map, 
## size of dots, size of labels
drc_plot <- ggplot() + geom_polygon(data = drc2, aes(long, lat, group=group), fill="whitesmoke") +
  coord_sf(xlim = c(28, 31), ylim = c(-0.7,2.1), expand=FALSE) +
  geom_path(data = drc2, aes(long, lat, group=group), color="darkgrey") +
  geom_point(data = drc_joined, aes(x=long, y=lat, size = weekly_confirm), color="red",  
             alpha=0.4) +
  geom_text_repel(data=drc_joined, aes(label=SCOnset, x=long, y=lat, size=2.2), 
                                 show.legend = FALSE, point.size=NA) +
  labs(title = "Week of Epidemic: {frame_time}", size = 'Weekly Cases') +
  scale_size(range = c(.1,25)) +
  transition_time(Week) 

## Animate the ggplot, including the speed and size of map
gganimate::animate(drc_plot, renderer=magick_renderer(), duration=20, width=800, height=800)

## Save the animation
anim_save(filename = "drc_animate_final.gif", animation = last_animation())



