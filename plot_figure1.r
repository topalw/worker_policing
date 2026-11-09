setwd("C:/Users/topalw/Desktop/ants/Behavioural/genetic_analyses")

library(tidyverse)
library(rworldxtra)
library(sf)
library(tmap)
library(dplyr)
library(tmaptools)
library(terra)
library(maptiles)

# colours 
colSupercolony <- '#ffba49'
colCleo <- '#20a39e'

# admixture proportions 
admix <- read.table('1.data/All_pruned_para_K2.admix.2.Q')
species <- read.table('1.data/species_assignment_mtDNA.txt',h=T)
# samples & species
s <- read.table('1.data/samples_species_bamorder.txt')
names(s) <- c('id','sp')
s %>% filter(sp == 'paralugubris') -> s

# order clusters 
as1 <- as.numeric(which.max(admix[54,])) # cluster 1 
as2 <- as.numeric(which.max(admix[84,])) # cluster 2
admix <- admix[,c(as1,as2)]
dom <- apply(admix,1,which.max)
# coordinates 
coord <- read.table('1.data/coordinates_dataset.txt',sep=' ',h=T)
names(coord) <- c('id','date','azimuth','size','alates','altitude','habitat','living',
                  'comments','e','n','lon','lat','sampled')
# sampled nests for behavioural proj 
b <- read.table('1.data/full_observed_nests.txt',sep='\t',h=T)
names(b) <- c('id','size','geometry','collection','sexual')
# filter out lugubris 
# manually add 22st304b which was missing from mtdna seqs
b %>% filter(id %in% c(species$id[species$sp == 'paralugubris'],'22st304b')) -> b


# change coordinates to N/E
geom_clean <- gsub("c\\(|\\)", "", b$geometry)
parts <- lapply(strsplit(geom_clean, ","), trimws)
b$e <- as.numeric(sapply(parts, `[`, 1))
b$n <- as.numeric(sapply(parts, `[`, 2))
pts <- st_as_sf(b, coords = c("e", "n"), crs = 2056, remove = FALSE)
pts <- st_transform(pts, crs = 4326)
crds <- st_coordinates(pts)
b$lon <- crds[, "X"]
b$lat <- crds[, "Y"]
b$sampled <- 1 
# merge 
coord <- bind_rows(coord,b)


# make europe map for top right 
data(countriesHigh)
world <- countriesHigh %>% st_as_sf()
world <- st_transform(world, crs = st_crs(3035))

europe_bbox <- st_bbox(c(xmin = -10, xmax = 40, ymax = 65, ymin = 30), crs = st_crs(4326))
europe_bbox <- st_bbox(st_transform(st_as_sfc(europe_bbox), 3035))

europe_map_cropped <- world %>% st_as_sf() %>% st_crop(europe_bbox)

supercol_sf <- st_sf(id='Supercolony',
                     geometry = st_sfc(st_point(c(6.188346, 46.54316)), crs = 4326))

t1 <- tm_shape(europe_map_cropped) +
  tm_polygons() +
  tm_graticules() +
  tm_shape(supercol_sf) + 
  tm_dots(fill='black',size=.7) + 
  tm_layout(bg.color = "white",
            space.color = "white",
            earth.boundary.lwd = 1)
pdf('global_map.pdf',8,6)
print(t1)
dev.off()

# study area 
pts_sf <- st_as_sf(coord, coords = c("lon", "lat"), crs = 4326)


# make distance matrix 
dist_mat <- st_distance(pts_sf)
dist_mat <- dist_mat[match(b$id,coord$id),match(b$id,coord$id)]
colnames(dist_mat) <- b$id
dim(dist_mat)
write.table(dist_mat,'1.data/distance_matrix.txt',row.names = F,quote=F)

# make box buffer 
BUFFER_M <- 1500   # try 2000-4000 m, tune to taste
# buffer in a metric CRS (LV95), then convert the box back to WGS84 for read_osm
bb <- st_bbox(st_transform(pts_sf, 2056))
bb["xmin"] <- bb["xmin"] - BUFFER_M
bb["xmax"] <- bb["xmax"] + BUFFER_M
bb["ymin"] <- bb["ymin"] - BUFFER_M
bb["ymax"] <- bb["ymax"] + BUFFER_M
bb <- st_bbox(bb, crs = 2056)
bbox_zoom <- st_bbox(st_transform(st_as_sfc(bb), 4326))
bbox_zoom_vect <- vect(st_as_sfc(bbox_zoom))      # bbox_zoom from before, in 4326
bg <- get_tiles(bbox_zoom_vect, provider = "Esri.WorldImagery", zoom = 15, crop = TRUE)
crds <- st_coordinates(pts_sf)
# the single outlier point, as its own sf object
outlier_sf <- st_sf(id = "Cleo",
                    geometry = st_sfc(st_point(c(6.240267, 46.56208)), crs = 4326))
supercol_sf <- st_sf(id='Supercolony',
                     geometry = st_sfc(st_point(c(6.188346, 46.54316)), crs = 4326))
t2 <- tm_shape(bg) +
  tm_rgb() +
  tm_shape(pts_sf) +
  #tm_shape(hull) +
  #tm_borders(col = colSupercolony, lwd = 2) +
  tm_shape(supercol_sf) +
  tm_dots(fill = colSupercolony, size = 2) +
  tm_shape(outlier_sf) +
  tm_dots(fill = colCleo, size = 2) +
   tm_graticules() 
# plot
pdf('local_map.pdf',8,6)
print(t2)
dev.off()

# supercolony plot 
coord2 <- coord[-grep('Cleo',coord$id),]
# add diagnostics for colouring
coord2$sampled <- ifelse(coord2$id %in% b$id,1,0)
coord2$genetic <- ifelse(coord2$id %in% s$id,1,0)
#admixture component 
coord2$admix <- NA
coord2$admix[match(s$id,coord2$id)] <- dom # 

### make genetic distance
# identify neighbours 
dist_mat <- st_distance(pts_sf)
colnames(dist_mat) <- coord$id
genetic_ids <- unique(coord2$id[coord2$genetic == 1 & coord2$id %in% colnames(dist_mat)])

idxs <- match(genetic_ids,coord2$id)
# for each sample find closest sequenced sample, get its name 
# its distance and its maximum Admixture cluster assignment 
coord2$nearest_genetic_id <- NA
coord2$nearest_physical_dist <-NA
coord2$nearest_max_cluster <- NA
coord2$nearest_cluster_id <- NA
for(i in 1:nrow(coord2)){
  min <- which.min(dist_mat[i,idxs])
  coord2$nearest_genetic_id[i] <- names(min)
  coord2$nearest_physical_dist[i] <- dist_mat[i,idxs][min]
  coord2$nearest_max_cluster[i] <- max(admix[match(coord2$nearest_genetic_id[i],s$id),])
  coord2$nearest_cluster_id[i] <-as.numeric(which.max(admix[match(coord2$nearest_genetic_id[i],s$id),]))
}

coord2[match(b$id[1:36],coord2$id),] %>% select(id,azimuth,size,lon,lat,
                                                contains('nearest')) -> tmp

write.table(tmp,'1.data/samples_with_cluster_assignment.txt',quote=F,row.names=F)
# get genetic distance of samples we want 
# read covariance 
cov <- as.matrix(read.table('1.data/All_pruned_para.cov'))
D <- as.matrix(dist(cov))

# extract genetic distance for sampled nests 
surrogates <- coord2$nearest_genetic_id[match(b$id,coord2$id)]
surrogates <- surrogates[1:36] # remove Cleo
gen_dist <- D[match(surrogates,s$id),match(surrogates,s$id)]
colnames(gen_dist) <- b$id[1:36]
write.table(gen_dist,'1.data/genetic_distance_matrix.txt',quote=F,row.names=F)

pts_sf <- st_as_sf(coord2, coords = c("lon", "lat"), crs = 4326)
BUFFER_M <- 500   # try 2000-4000 m, tune to taste
# buffer in a metric CRS (LV95), then convert the box back to WGS84 for read_osm
bb <- st_bbox(st_transform(pts_sf, 2056))
bb["xmin"] <- bb["xmin"] - BUFFER_M
bb["xmax"] <- bb["xmax"] + BUFFER_M
bb["ymin"] <- bb["ymin"] - BUFFER_M
bb["ymax"] <- bb["ymax"] + BUFFER_M
bb <- st_bbox(bb, crs = 2056)
bbox_zoom <- st_bbox(st_transform(st_as_sfc(bb), 4326))
bbox_zoom_vect <- vect(st_as_sfc(bbox_zoom))    
bg <- get_tiles(bbox_zoom_vect, provider = "Esri.WorldImagery", zoom = 17, crop = TRUE)
crds <- st_coordinates(pts_sf)
# point sets 
base_pts    <- pts_sf[coord2$sampled == 0 & coord2$genetic == 0, ]
admix1_pts  <- pts_sf[coord2$genetic == 1 & coord2$admix == 1, ]
admix2_pts  <- pts_sf[coord2$genetic == 1 & coord2$admix == 2, ]
sampled_pts <- pts_sf[coord2$sampled == 1, ]


t2 <- tm_shape(bg) +
  tm_rgb() +
  tm_shape(base_pts) +
  tm_dots(fill = "grey", size = 0.3) +
  tm_shape(admix1_pts) +
  tm_dots(fill = "#E46B6B", fill_alpha=.8, size = 1) +
  tm_shape(admix2_pts) +
  tm_dots(fill = "#527FBD", fill_alpha = .8, size = 1) +
  tm_shape(sampled_pts) +
  tm_dots(fill=colSupercolony, size = .6) +  
  tm_scalebar(text.size=.6, text.color='white') + tm_graticules()
t2

# plot
pdf('supercolony_map.pdf',8,6)
print(t2)
dev.off()



