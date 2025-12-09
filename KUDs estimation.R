###################################################################################################################
# Miguel Gandra || m3gandra@gmail.com || November 2023 ############################################################
###################################################################################################################

# Acoustic Telemetry unified script - multispecies analyses

#  • computation of spatiotemporal overlap index
#  • estimation of kernel utilization distributions (KUDs)
#  • estimation of movement metrics (distance, ROM, etc.)
#  • diel and seasonal analysis


#######################################################################################################
# Automatically install required libraries  ###########################################################
#######################################################################################################

if(!require(readxl)){install.packages("readxl"); library(readxl)}
if(!require(adehabitatHR)){install.packages("adehabitatHR"); library(adehabitatHR)}
if(!require(raster)){install.packages("raster"); library(raster)}
if(!require(sf)){install.packages("rgdal"); library(sf)}
if(!require(maptools)){install.packages("maptools"); library(maptools)}
if(!require(plotrix)){install.packages("plotrix"); library(plotrix)}
if(!require(plyr)){install.packages("plyr"); library(plyr)}
if(!require(gdistance)){install.packages("gdistance"); library(gdistance)}
if(!require(zoo)){install.packages("zoo"); library(zoo)}
if(!require(proxy)){install.packages("proxy"); library(proxy)}
if(!require(mapdata)){install.packages("mapdata"); library(mapdata)}
if(!require(geosphere)){install.packages("geosphere"); library(geosphere)}
if(!require(glmmTMB)){install.packages("glmmTMB"); library(glmmTMB)}
if(!require(car)){install.packages("car"); library(car)}
if(!require(DHARMa)){install.packages("DHARMa"); library(DHARMa)}
if(!require(effects)){install.packages("effects"); library(effects)}
if(!require(fitdistrplus)){install.packages("fitdistrplus"); library(fitdistrplus)}
#library(poseidon)
source("http://highstat.com/Books/BGS/GAMM/RCodeP2/HighstatLibV6.R")
source("./Code/mcp_function.R")


#######################################################################################################
# Load Data  ##########################################################################################
#######################################################################################################

# import data
data_files <- list.files("./COAs", full.names=T)
data <- lapply(data_files, function(x) read.csv2(x, header=T))
names(data) <- sub("\\..*", "", list.files("./COAs"))
data <- mapply(function(data, filename){data$file<-filename; return(data)}, data=data, filename=names(data), SIMPLIFY=F)

# import info of tagged fish 
fish_info <- read_excel("./Data/Fish_info_v5.xlsx", sheet=1, trim_ws=T)
fish_info <- data.frame(fish_info)
fish_info$tagging_date <- as.POSIXct(fish_info$tagging_date, format="%Y-%m-%d", tz="UTC")
fish_info$length_cm <- as.numeric(fish_info$length_cm)

# import methodological info
methodological_info <- read_excel("./Data/Methodological_info_v3.xlsx", sheet=1, trim_ws=T)
methodological_info <- data.frame(methodological_info)
files_names <- sub("\\..*", "", list.files("./Data/COAs"))
methodological_info <- methodological_info[order(match(methodological_info$file, files_names)),]

# check if there is any missing entries
file_refs <- unlist(lapply(data, function(x) unique(x$file)))
names(file_refs) <- NULL
if(any(!file_refs %in% methodological_info$file)){
  missing_files <- file_refs[!file_refs %in% methodological_info$file]
  missing_files <- paste0(paste("•", missing_files), collapse="\n")
  stop(paste0("Missing data in methodological info file\n", missing_files))
  missing_files
}

#import receiver info
receivers_list <- read.csv2("./Receivers_info.csv", header=T)
receiver_info <- aggregate(receivers_list$receiver, by=list(receivers_list$file), function(x) length(unique(x)))
colnames(receiver_info) <- c("RefID", "N_receivers")

# import coastline layer (https://www.eea.europa.eu/data-and-maps/data/eea-coastline-for-analysis-2)
coastline <- read_sf(dsn="./EEA_Coastline_Polygon_Shape", layer="EEA_Coastline_20170228")
coastline <- st_transform(coastline, st_crs(4326))

# format data
data <- lapply(data, function(x) {x$transmitter<-as.factor(x$transmitter); return(x)})
data <- lapply(data, function(x) {x$timebin<-as.character(x$timebin); return(x)})
data <- lapply(data, function(x) {x$timebin<-as.POSIXct(x$timebin, format="%Y-%m-%d %H:%M:%S", tz="UTC"); return(x)})
data <- lapply(data, function(x) {x$longitude<-as.numeric(x$longitude); return(x)})
data <- lapply(data, function(x) {x$latitude<-as.numeric(x$latitude); return(x)})

# order by ID and timebin
data <- lapply(data, function(x) {x[order(x$transmitter, x$timebin),]; return(x)})

# reorder columns
data <- lapply(data, function(x) x[,c("file", "timebin", "transmitter", "longitude", "latitude", "detections", "species")])


#######################################################################################################
# Function to assign reproductive season (spawning vs resting)  #######################################
#######################################################################################################

# assignReprodPeriod <- function(data) {
#   
#   # grab spawning periods defined on the methodological file
#   print(unique(data$file))
#   reprod_periods <- methodological_info[methodological_info$file==unique(data$file),]
#   if(is.na(reprod_periods$spawning_start)){
#     data$reprod_season <- NA; return(data)
#   }
#   spawning_start <- reprod_periods$spawning_start
#   spawning_end <- reprod_periods$spawning_end
#   
#   # create vector with all dates and assign reprod season
#   start <- lubridate::floor_date(min(data$timebin), unit="day")
#   end <- lubridate::ceiling_date(max(data$timebin), unit="day")
#   date_vec <- data.frame("date"=seq.POSIXt(start, end, by="day"))
#   date_vec$reprod_season <- getReprodPeriod(date_vec$date, spawning.start=spawning_start,
#                                             spawning.end=spawning_end, format="%d/%b")
#   
#   # distinguish betweeen non-consecutive seasons
#   consec_seasons <- rle(as.character(date_vec$reprod_season))
#   season_counts <- table(consec_seasons$values)
#   season_counts <- unlist(sapply(season_counts, function(x) 1:x))
#   season_counts <- season_counts[order(season_counts)]
#   date_vec$season_id  <- paste0(date_vec$reprod_season, "/", rep(sprintf("%02d", season_counts), consec_seasons$lengths))
#   
#   # assign numbered seasons to each data record
#   data$date <- strftime(data$timebin, "%Y-%m-%d", tz="UTC")
#   data <- plyr::join(data, date_vec[,c("date","season_id")], by="date", type="left")
#   s_order <- order(as.numeric(gsub("\\D", "", unique(data$season_id))))
#   data$season_id <- factor(data$season_id, levels=unique(data$season_id)[s_order])
#   colnames(data)[which(colnames(data)=="season_id")] <- "reprod_season"
#   return(data)
# }
# 
# data <- lapply(data, assignReprodPeriod)



#######################################################################################################
# Function to assign diel phase (day vs night)  #######################################################
#######################################################################################################

# calculateDielPhase <- function(data) {
#   
#   print(unique(data$file))
#   coordinates <- SpatialPoints(data[,c("longitude","latitude")], proj4string=CRS("+proj=longlat +datum=WGS84"))
#   coordinates <- cbind(data$longitude, data$latitude)
#   data$timeofday <- getDielPhase(data$timebin, phases=4, coordinates=coordinates)
#   return(data)
# }
# 
# data <- lapply(data, calculateDielPhase)



#######################################################################################################
# Remove detection prior to tagging  ##################################################################
#######################################################################################################

removeFalseDetections <- function(data) {
  cat(unique(data$file))
  raw_total <- nrow(data)
  tagging_info <- data.frame("file"=unique(data$file), "transmitter"=levels(data$transmitter))
  tagging_info <- plyr::join(tagging_info, fish_info[,c("file", "transmitter", "tagging_date")], by=c("file", "transmitter"), type="left", match="first")
  tagging_dates <- tagging_info$tagging_date
  data_transmitter <- split(data, f=data$transmitter)
  data_transmitter <- mapply(function(data, tagdate) subset(data, data$timebin>=tagdate), 
                             data=data_transmitter, tagdate=tagging_dates, SIMPLIFY=F) 
  data_filtered <- do.call("rbind", data_transmitter)
  data_filtered$transmitter <- droplevels(data_filtered$transmitter)
  deleted <- raw_total - nrow(data_filtered)
  cat(paste0(" - ", deleted, " detections removed\n"))
  return(data_filtered)
}

data <- lapply(data, removeFalseDetections)


# discard spurious detections in Epinephelus_marginatus1 dataset
original_rows <- nrow(data$Epinephelus_marginatus1 )
data$Epinephelus_marginatus1 <- data$Epinephelus_marginatus1[data$Epinephelus_marginatus1$latitude>38.1,]
paste(original_rows-nrow(data$Epinephelus_marginatus1), "detection(s) discarded\n")
  

#######################################################################################################
# Summary Function  ###################################################################################
#######################################################################################################

acousticSummary <- function(data) {
  
  print(unique(data$file))
  tagging_info <- data.frame("file"=unique(data$file), "transmitter"=levels(data$transmitter))
  tagging_info <- join(tagging_info, fish_info[,c("file", "transmitter", "tagging_date")], by=c("file", "transmitter"), type="left", match = "first")
  tagging_dates <- tagging_info$tagging_date
  tagdates <- strftime(tagging_dates, format="%d/%m/%Y", tz="UTC")
  n_receivers <- receiver_info$N_receivers[receiver_info$RefID==unique(data$file)]
  
  if(any(is.na(tagging_dates))){return(list("complete_table"=NA, "summary_table"=NA))}
  
  # retrieve last detections dates
  last_dates <- tapply(X=data$timebin, INDEX=data$transmitter, FUN=max)
  last_dates <- as.POSIXct(last_dates, origin='1970-01-01', tz="UTC")
  last_dates <- strftime(last_dates, format="%Y-%m-%d", tz="UTC")
  last_dates <- as.POSIXct(last_dates, format="%Y-%m-%d", tz="UTC")
  last_dates_print <- strftime(last_dates, format="%d/%m/%Y", tz="UTC")
  
  # days between 1st and last detection (Di)
  getTimeSeqs <- function(start, end) {if(is.na(end)) {return(NA)}else{seq.POSIXt(start, end, by="day")}}
  timeseqs <- mapply(getTimeSeqs, start=tagging_dates, end=last_dates, SIMPLIFY=F)
  Di <- as.numeric(lapply(timeseqs, function(x) length(x[!is.na(x)])))
  Di[Di==0] <- NA
  
  
  # calculate nº monitored hours
  monitored_hours <- mapply(function(start, end){difftime(end, start, units="hours")},
                            start=tagging_dates, end=last_dates, SIMPLIFY=F)
  monitored_hours <- as.numeric(monitored_hours)
  monitored_hours[monitored_hours==0] <- NA
  
  
  # days with detections (Dd)
  data$date <- strftime(data$timebin, format="%d-%m-%Y", tz="UTC")
  Dd <- by(data, data$transmitter, function (x) length(unique(x$date)))
  Dd <- as.integer(Dd)
  data <- data[,-ncol(data)]
  
  # complete IR
  Ir <- round(Dd/Di,2)
  
  # individual detections
  individual_detections <- aggregate(data$detections, by=list(data$transmitter), sum)$x
  
  
  # create results table
  complete_table <- data.frame("File"=unique(data$file), "Species"=unique(data$species), "Transmitter"=levels(data$transmitter), "Tagging_date"=tagdates, 
                               "Last_detection"=last_dates_print, "Detections"=individual_detections, "Monitored_hours"=monitored_hours,
                               "Di"=Di, "Dd"=Dd, "IR"=Ir, check.names=F, row.names=NULL)
  
  # calculate summary stats
  timespan <- paste(strftime(min(tagging_dates), format="%m/%Y", tz="UTC"), "-", strftime(max(last_dates), format="%m/%Y", tz="UTC"))
  average_Di <- paste(round(mean(Di, na.rm=T)), "±", round(std.error(Di)))
  average_Dd <- paste(round(mean(Dd, na.rm=T)), "±", round(std.error(Dd)))
  average_Ir <- paste(sprintf("%.2f", mean(Ir, na.rm=T)), "±", sprintf("%.2f", std.error(Ir)))
  
  summary_table <- data.frame("Species"=unique(data$species), "Timespan"=timespan, "Sample Size"=nlevels(data$transmitter),
                              "Nº detections"=sum(data$detections), "Nº receivers"=n_receivers, "Di"=average_Di, "Dd"=average_Dd, "IR"=average_Ir)
  
  # retrieve mean stats (non-formatted)
  mean_values <- data.frame("File"=unique(data$file), "Species"=unique(data$species), "SampleSize"=nlevels(data$transmitter), "NReceivers"=n_receivers, 
                            "NDetections"=sum(data$detections), "Di"=mean(Di, na.rm=T), "Dd"=mean(Dd, na.rm=T), "IR"=mean(Ir, na.rm=T))
  
  # return results
  return(list("complete_table"=complete_table, "summary_table"=summary_table, "mean_values"=mean_values))
}



#######################################################################################################
## Compute basic stats ################################################################################
#######################################################################################################

# calculate simple space-use metrics
results <- lapply(data, acousticSummary)


# aggregate summarized results
summary_table <- do.call("rbind", lapply(results, function(x) x$summary_table))
summary_table <- summary_table[order(summary_table$Species),]

# aggregate mean values by dataset
species_data <- do.call("rbind", lapply(results, function(x) x$mean_values))

# aggregate complete data (all individuals)
complete_data <- do.call("rbind", lapply(results, function(x) x$complete_table))
rownames(complete_data) <- NULL
colnames(complete_data)[which(colnames(complete_data)=="Di")] <- "Tp"

# remove individuals with < 5 detections
complete_data <- complete_data[complete_data$Detections>=5,]


#######################################################################################################
## Discard datasets ###################################################################################
#######################################################################################################

discarded_datasets <- c("Diplodus_sargus3", "Diplodus_vulgaris", "Epinephelus_marginatus3", 
                        "Scorpaena_porcus", "Scorpaena_scrofa", "Serranus_cabrilla", 
                        "Serranus_scriba", "Xyrichtys_novacula")

data <- data[!names(data) %in% discarded_datasets]
complete_data <- complete_data[!complete_data$File %in% discarded_datasets,]
receivers_list <- receivers_list[!receivers_list$file %in% discarded_datasets,]
methodological_info <- methodological_info[!methodological_info$file %in% discarded_datasets,]
summary_table <- summary_table[!rownames(summary_table) %in% discarded_datasets,]

#write.csv2(summary_table, file="./summary_table.csv",  row.names=F, fileEncoding="Windows-1252")



#######################################################################################################
## Grab fish info #####################################################################################
#######################################################################################################

colnames(fish_info) <- tools::toTitleCase(colnames(fish_info))
complete_data <- join(complete_data, fish_info[,c("File","Transmitter", "Length_cm", "Tag")], 
                      by=c("File","Transmitter"), type="left")
ordered_cols <- c("File", "Species", "Transmitter", "Length_cm", "Tag", "Tagging_date",
                  "Last_detection", "Detections", "Monitored_hours", "Tp", "Dd", "IR")
complete_data <- complete_data[,ordered_cols]



#######################################################################################################
## Grab species and methodological variables ##########################################################
#######################################################################################################

# Assign unique array ID based on receiver coordinates
receivers_list$longitude <- as.numeric(receivers_list$longitude)
receivers_list$latitude <- as.numeric(receivers_list$latitude)
arrays_list <- aggregate(receivers_list[,c("longitude","latitude")], by=list(receivers_list$file), mean)
arrays_list$longitude <- round(arrays_list$longitude)
arrays_list$latitude <- round(arrays_list$latitude)
colnames(arrays_list)[1] <- "File"
unique_arrays <- unique(arrays_list[,c("longitude","latitude")])
unique_arrays$ArrayID <- 1:nrow(unique_arrays) 
arrays_list <- plyr::join(arrays_list, unique_arrays, by=c("longitude","latitude"), type="left")
arrays_list <- arrays_list[order(arrays_list$ArrayID),]
complete_data <- plyr::join(complete_data, arrays_list[,c("File","ArrayID")], by="File", type="left")

# Calculate nº receivers
array_pts <- split(receivers_list, f=list(receivers_list$file))
n_receivers <- lapply(array_pts, function(x) data.frame("File"=unique(x$file), "NReceivers"=nrow(unique(x[,c("longitude","latitude")]))))
n_receivers <- do.call("rbind", n_receivers)
complete_data <- plyr::join(complete_data, n_receivers, by="File", type="left")

# Calculate monitored areas (in km2)
array_pts <- split(receivers_list, f=list(receivers_list$file))
array_pts <- lapply(array_pts, function(x) st_multipoint(as.matrix(x[,c("longitude","latitude")])))
array_pts <- lapply(array_pts, function(x) st_sfc(x, crs=st_crs(4326)))
epsg_codes <- methodological_info
epsg_codes <- epsg_codes[order(match(epsg_codes$file, names(array_pts))),]
array_pts <- mapply(function(pts, code){st_transform(pts, st_crs(code))}, pts=array_pts, code=epsg_codes$epsg_code)
#detec_ranges <- methodological_info$detection_range
detec_ranges <- 250
array_areas <- mapply(function(pts, range){st_buffer(pts, dist=range, byid=T)}, pts=array_pts, range=detec_ranges, SIMPLIFY=F)
array_areas <- lapply(array_areas, function(x) round(st_area(x)/1000000, 2))
array_areas <- reshape2::melt(array_areas)
colnames(array_areas) <- c("MonitoredArea_km2", "File")
complete_data <- join(complete_data, array_areas, by="File", type="left")

# Estimate spatial scale
array_length <- lapply(array_pts, function(x) round(max(dist(st_coordinates(x)))/1000,2))
array_length <- reshape2::melt(array_length)
colnames(array_length) <- c("SpatialScale_km", "File")
complete_data <- join(complete_data, array_length, by="File", type="left")

# Estimate receiver density (nº receivers per square km of MCP) #'
array_pts_sp <- lapply(array_pts, function(x) as(x, "Spatial"))
array_pts_sp <- lapply(array_pts_sp, function(x) sp::SpatialPoints(x@coords))
array_mcp <- lapply(array_pts_sp, function(x) mcp_receivers(x, percent=100, unin="m", unout="km2"))
# remove land overlaps
pb <- txtProgressBar(min=1, max=length(array_mcp), initial=0, style=3)
for(i in 1:length(array_mcp)){
  setTxtProgressBar(pb, i)
  epsg_code <- epsg_codes$epsg_code[i]
  coastline_projected <- st_transform(coastline, st_crs(epsg_code))
  coastline_study <- st_crop(st_geometry(coastline_projected), extent(array_mcp[[i]]))
  array_mcp[[i]] <- st_as_sf(array_mcp[[i]])
  st_crs(array_mcp[[i]]) <- st_crs(epsg_code)
  if(length(coastline_study)>0){
    coastline_study <- st_union(st_geometry(coastline_study))
    array_mcp[[i]] <- st_difference(st_geometry(array_mcp[[i]]), st_geometry(coastline_study))
  }
}
close(pb)
final_mcp <- reshape2::melt(lapply(array_mcp, function(x) st_area(x)/1000000))
colnames(final_mcp) <- c("MCP_km2", "File")
final_mcp$MCP_km2 <- as.numeric(final_mcp$MCP_km2)
final_mcp$MCP_km2 <- round(final_mcp$MCP_km2, 4)
complete_data <- join(complete_data, final_mcp, by="File", type="left")
complete_data$ReceiverDensity <- round(complete_data$NReceivers/complete_data$MCP_km2, 2)

# Get species traits
species_data <- read_excel("./Data/Species_traits_revised.xlsx", sheet=1)
species_data <- species_data[!is.na(species_data$FBname),]
sampling_cols <- which(grepl(".", colnames(species_data), fixed=T))
species_data <- species_data[,-sampling_cols]
excluded_cols <- c("FeedingType", "DietTroph", "MaxLengthSL", "RepGuild1", "DepthMin", "ComDepthMin", "ComDepthMax")
species_data <- species_data[,-which(colnames(species_data) %in% excluded_cols)]
complete_data <- join(complete_data, species_data, by="Species", type="left")



#######################################################################################################
# Function to remove individuals with less than 5 records #############################################
# required to avoid errors in the kernel_UD function

cleanData <- function(data) {
  data$transmitter <-droplevels(data$transmitter)
  for (i in levels(data$transmitter)) {
    coords <- st_coordinates(data[data$transmitter==i,])
    if (nrow(coords) < 5 || is.null(nrow(coords))) {
      data <- subset(data, data$transmitter!=i)
    }
  }
  data$transmitter <- droplevels(data$transmitter)
  return(data)
}


#######################################################################################################
# Function to substract area on land ##################################################################

subtractLand <- function(land_layer) {
  
  for (i in 1:length(data)){
    if(is.null(intersect(data[i,], land_layer))) {next}
    data$area[i] <-  data$area[i] - sum(area(intersect(data[i,], land_layer)))/1000000
  }  
  return(data)
}



#######################################################################################################
# Calculate home-ranges ###############################################################################
#######################################################################################################

calculateUDs <- function(data, bandwidth, grid_res=100) {
  
  cat(unique(data$file))
  epsg_code <- methodological_info$epsg_code[methodological_info$file==unique(data$file)]
  if(is.na(epsg_code)){return(list("complete_KUDs"=NA, "summary_table"=NA, "mean_values"=NA))}
  
  # convert data to point feature and project coordinates
  coord_cols <- c("longitude", "latitude")
  data_cols <- colnames(data)[!colnames(data) %in% coord_cols]
  spatial_data <- st_as_sf(data, coords=coord_cols, crs=st_crs(4326))
  spatial_data <- st_transform(spatial_data, st_crs(epsg_code))
  
  # clean spatial data
  spatial_data <- cleanData(spatial_data)
  if(nrow(spatial_data)==0){cat("\n"); return(list("complete_KUDs"=NA, "summary_table"=NA, "mean_values"=NA))}
  
  # estimate bounding box
  if(nrow(unique(st_coordinates(spatial_data)))<=2){
    bbox <- st_bbox(st_buffer(spatial_data, dist=500)) 
  }else{
    bbox <- st_bbox(spatial_data)
  }
  
  # check bounding box dimensions
  if(abs(bbox[3]-bbox[1])<grid_res){
    bbox[1] <- bbox[1] - grid_res/2
    bbox[3] <- bbox[3] + grid_res/2
  }
  if(abs(bbox[4]-bbox[2])<grid_res){
    bbox[2] <- bbox[2] - grid_res/2
    bbox[4] <- bbox[4] + grid_res/2
  }
  
  # create spatial grid
  my_grid <- raster(extent(bbox)*2, res=grid_res)
  my_grid <- as(my_grid, 'SpatialPixels')
  
  # set kernel smoothing factor
  selected_bandwidth <- bandwidth
  
  # compute KUDs (increase grid extent if necessary)
  spatial_data <- as(spatial_data, "Spatial")
  kud_complete <- kernelUD(spatial_data[,"transmitter"], h=selected_bandwidth, grid=my_grid)
  k95_complete <- try(getverticeshr(kud_complete, percent=95, unin="m", unout="km2"), silent=T)
  counter <- 0
  multiplier <- seq(2, 5, by=0.5)
  while(class(k95_complete)=="try-error") {
    counter <- counter+1
    my_grid <- raster(extent(raster(my_grid))*multiplier[counter], res=grid_res)
    my_grid <- as(my_grid, 'SpatialPixels')
    kud_complete <- kernelUD(spatial_data[,"transmitter"], h=selected_bandwidth, grid=my_grid)
    k95_complete <- try(getverticeshr(kud_complete, percent=95, unin="m", unout="km2"), silent=T)
    if(class(k95_complete)!="try-error"){cat(paste0(" - Grid extended (x", multiplier[counter], ")"))}
  }
  k50_complete <- getverticeshr(kud_complete, percent=50, unin="m", unout="km2")
  #complete_results <- list("kernelUDs"=kud_complete, "IDs"=k50_complete$id, "K50"=k50_complete$area, "K95"=k95_complete$area)
  complete_results <- data.frame("File"=unique(data$file), "ID"=k50_complete$id, "KUD50"=k50_complete$area, "KUD95"=k95_complete$area,
                                 "Species"=unique(data$species), "Week"=unique(data$week))
  gc()
  
  # summary table
  average_k50 <-  paste(sprintf("%.2f", mean(k50_complete$area, na.rm=T)), "±", sprintf("%.2f", std.error(k50_complete$area)))
  average_k95 <- paste(sprintf("%.2f", mean(k95_complete$area, na.rm=T)), "±", sprintf("%.2f", std.error(k95_complete$area)))
  summary_table <- data.frame("Species"=unique(data$species), "Nº fish"=nrow(k50_complete), "K50_mean"=average_k50, "k95_mean"=average_k95, "bandwidth"=selected_bandwidth) 
  
  # mean stats (non-formatted)
  mean_values <- data.frame("Species"=unique(data$species), "k50_mean"=mean(k50_complete$area, na.rm=T), "k95_mean"=mean(k95_complete$area, na.rm=T))
  
  # return results
  cat("\n")
  return(list("complete_KUDs"=complete_results, "summary_table"=summary_table, "mean_values"=mean_values))
  
}


#######################################################################################################
# Calculate home-ranges (complete periods) ############################################################
#######################################################################################################

start.time <- Sys.time()
kud_results <- lapply(data, calculateUDs, bandwidth=200, grid_res=100)
end.time <- Sys.time()
end.time - start.time
#save(kud_results, file=c("./kud_results_h200.RData"))
#load("./kud_results_h200.RData")


# save summary table
kud_summary_table <- do.call("rbind", lapply(kud_results, function(x) x$summary_table))
kud_summary_table <- kud_summary_table[order(kud_summary_table$Species),]
write.csv2(kud_summary_table, "./kud_summary_table_h200.csv", row.names=F)

# retrieve results by ID
kud_values <- lapply(kud_results, function(x) x$complete_KUDs)
kud_values <- lapply(kud_values, function(x) data.frame("Transmitter"=x$IDs, "KUD50"=x$K50, "KUD95"=x$K95))
kud_values <- mapply(function(table, dataset){table$File<-dataset; return(table)}, table=kud_values, dataset=names(kud_results), SIMPLIFY=F)
kud_values <- do.call("rbind", kud_values)
rownames(kud_values) <- NULL
kud_values$KUD50 <- round(kud_values$KUD50, 3)
kud_values$KUD95 <- round(kud_values$KUD95, 3)
complete_data <- join(complete_data, kud_values, by=c("Transmitter","File"), type="left")


#######################################################################################################
# Calculate home-ranges (weekly) ######################################################################
#######################################################################################################

# assign week (of the year) and split data
data <- lapply(data, function(x) {x$week <- strftime(x$timebin, "%W/%Y", tz="UTC"); return(x)})
data_weekly <- lapply(data, function(x) split(x, f=x$week))

# calculate KUD areas for each individual and week
kud_week_results <- vector("list", length(data_weekly))
names(kud_week_results) <- names(data_weekly)
start.time <- Sys.time()
for(i in 1:length(data_weekly)){
  kud_week_results[[i]] <- lapply(data_weekly[[i]], calculateUDs, bandwidth=200, grid_res=100)
}
end.time <- Sys.time()
end.time - start.time

# aggregate results
kud_weekly_table <- list()
for(i in 1:length(data_weekly)){
  kud_week <- lapply(kud_week_results[[i]], function(x) x$complete_KUDs)
  empty_weeks <- which(unlist(lapply(kud_week, function(x) all(is.na(x)))))
  if(length(empty_weeks)>0){kud_week <- kud_week[-empty_weeks]}
  kud_week <- do.call("rbind", kud_week)
  rownames(kud_week) <- NULL
  kud_week$File <- names(data_weekly)[i]
  kud_weekly_table[[i]] <- kud_week
}
kud_weekly_table <- do.call("rbind", kud_weekly_table)

# join information
colnames(kud_weekly_table)[2] <- "Transmitter"
kud_weekly_table <- plyr::join(kud_weekly_table, complete_data[,c("File", "Transmitter", "Length_cm")],
                               by=c("File", "Transmitter"), type="left")
kud_weekly_table <- kud_weekly_table[,c("File", "Transmitter", "KUD50", "KUD95", "Week", "Length_cm")]


# format
kud_weekly_table$KUD50 <- round(kud_weekly_table$KUD50, 3)
kud_weekly_table$KUD95 <- round(kud_weekly_table$KUD95, 3)
kud_weekly_table <- kud_weekly_table[order(kud_weekly_table$File, kud_weekly_table$Transmitter, kud_weekly_table$Week),]


# save results
write.csv2(kud_weekly_table, "./Data/KUDs by individual-week.csv", row.names=F)



#######################################################################################################
# Calculate home-ranges (spawning vs non-spawning) ####################################################
#######################################################################################################

# data_reprod <- lapply(data, function(x) split(x, f=x$reprod_season))
# 
# # keep only datasets that include at least one resting and one spawing season
# selected_data <-  unlist(lapply(data_reprod, function(x) any(grepl("resting", names(x), fixed=T)) & any(grepl("spawning", names(x), fixed=T))))
# data_reprod <- data_reprod[selected_data]
# 
# # calculate KUD areas for each individual and reprod season
# kud_reprod_results <- vector("list", length(data_reprod))
# names(kud_reprod_results) <- names(data_reprod)
# start.time <- Sys.time()
# for(i in 1:length(data_reprod)){
#   kud_reprod_results[[i]] <- lapply(data_reprod[[i]], calculateUDs, bandwidth=200, grid_res=100)
# }
# end.time <- Sys.time()
# end.time - start.time
# 
# # aggregate results
# kud_reprod_table <- list()
# for(i in 1:length(data_reprod)){
#   kud_reprod <- lapply(kud_reprod_results[[i]], function(x) x$complete_KUDs)
#   empty_seasons <- which(unlist(lapply(kud_reprod, function(x) all(is.na(x)))))
#   if(length(empty_seasons)>0){kud_reprod <- kud_reprod[-empty_seasons]}
#   kud_reprod <- lapply(kud_reprod, function(x) data.frame("Transmitter"=x[[2]], "KUD50"=x[[3]], "KUD95"=x[[4]]))
#   kud_reprod <- mapply(function(table, season){table$Reprod_season <- season; return(table)}, table=kud_reprod, season=names(kud_reprod), SIMPLIFY=F)
#   kud_reprod <- do.call("rbind", kud_reprod)
#   rownames(kud_reprod) <- NULL
#   kud_reprod$File <- names(data_reprod)[i]
#   kud_reprod_table[[i]] <- kud_reprod
# }
# kud_reprod_table <- do.call("rbind", kud_reprod_table)
# 
# 
# # join information
# kud_reprod_table <- join(kud_reprod_table, complete_data[,c("File", "Transmitter", "Length_cm")],
#                          by=c("File", "Transmitter"), type="left")
# kud_reprod_table <- kud_reprod_table[,c("File", "Transmitter", "KUD50", "KUD95", "Reprod_season", "Length_cm")]
# 
# 
# # format
# kud_reprod_table$KUD50 <- round(kud_reprod_table$KUD50, 3)
# kud_reprod_table$KUD95 <- round(kud_reprod_table$KUD95, 3)
# kud_reprod_table <- kud_reprod_table[order(kud_reprod_table$File, kud_reprod_table$Transmitter, kud_reprod_table$Reprod_season),]
# 
# 
# # save results
# write.csv2(kud_reprod_table, "./Data/KUDs by individual-season.csv", row.names=F)



#######################################################################################################
# Calculate distances traveled #######################################################################
#######################################################################################################

calculateDistances <- function(data, grid.resolution) {
  
  # show current dataset in console
  cat(paste0(unique(data$file), "\n"))
  
  # retrieve epsg code
  epsg_code <- methodological_info$epsg_code[methodological_info$file==unique(data$file)]
  
  coords <- data[,c("longitude", "latitude")]
  sp_points <- st_multipoint(as.matrix(coords))
  sp_points <- st_sfc(sp_points,crs=st_crs(4326))
  bbox <- extent(st_bbox(sp_points))*1.5
  coastline_study <- st_crop(st_geometry(coastline), bbox)
  coastline_projected <- st_transform(coastline_study, st_crs(epsg_code))
  coastline_projected <- as(coastline_projected, "Spatial")
  
  discard_tags <- names(table(data$transmitter))[which(table(data$transmitter)<3)]
  data <- data[!data$transmitter %in% discard_tags,]
  data$transmitter <- droplevels(data$transmitter)
  data_transmitter <- split(data, f=data$transmitter)
  transmitter_results <- list()
  
  # calculate linear distances (if no land intersections)
  if(length(coastline_study)==0){
    dists<-c()
    for (t in 1:length(data_transmitter)) {
      # print progress to console
      cat(paste0("  • id = ", names(data_transmitter)[[t]], "\n"))
      # if individual doesn't have any detections, jump to next
      if(nrow(data_transmitter[[t]])<=0){next}
      # else, if individual has only a single detection, return point geometry
      if(nrow(data_transmitter[[t]])==1){transmitter_results[[t]]<-NA; next}
      # else, calculate great circle distances
      coords <- as.matrix(data_transmitter[[t]][,c("longitude","latitude")])
      transmitter_results[[t]] <- c(geosphere::distVincentyEllipsoid(coords), NA)
    }
    # else calculate shortest in-water paths (whenever land is intersected)
  }else{
    # compute transition matrix if not already available 
    transition_file <- paste0("./Data/trCost/trCost", grid.resolution, "m_", unique(data$file), ".RData")
    if(!file.exists(transition_file)){
      cat(paste0(" - creating transition object\n"))
      projected_pts <- st_transform(sp_points, st_crs(epsg_code))
      bbox <- extent(st_bbox(projected_pts))*1.5
      template_raster <- raster(bbox, res=grid.resolution)
      template_raster[] <- 0
      land_raster <- rasterize(coastline_projected, template_raster, update=T)
      values(land_raster)[values(land_raster)>0] <- 10000
      values(land_raster)[values(land_raster)==0] <- 1
      trCost <- transition(1/land_raster, transitionFunction=mean, directions=16)
      trCost <- geoCorrection(trCost, type="c")
      save(trCost, file=transition_file )
    }else{
      cat(paste0(" - transition object loaded\n"))
      load(transition_file)
    }
    # iterate over each individual
    for (t in 1:length(data_transmitter)) {
      # print progress to console
      cat(paste0("  • id = ", names(data_transmitter)[[t]], "\n"))
      coords <- as.matrix(data_transmitter[[t]][,c("longitude","latitude")])
      coords_proj <- st_sfc(st_multipoint(coords), crs=st_crs(4326))
      coords_proj <- st_transform(coords_proj, crs=st_crs(epsg_code))
      coords_proj <- st_coordinates(coords_proj)[,1:2]
      dists <- geosphere::distVincentyEllipsoid(coords)
      for(r in 1:nrow(coords)){
        if(r==nrow(coords)){break}
        segment <- st_linestring(rbind(coords[r,], coords[r+1,]))
        segment <- st_sfc(segment, crs=st_crs(4326))
        is_pt <- all(coords[r,] == coords[r+1,])
        in_land <- lengths(sf::st_intersects(segment, sf::st_as_sf(coastline_study), sparse=T))>0
        # if segment overlaps land and is not a point, calculate shortest in-water path
        if(in_land==T & is_pt==F & floor(dists[r])>grid.resolution){
          segment <- gdistance::shortestPath(trCost, coords_proj[r,], coords_proj[r+1,], output="SpatialLines")
          segment <- st_as_sf(segment)
          sf::st_crs(segment) <- st_crs(epsg_code)
          segment_wgs84 <- st_transform(segment, crs=st_crs(4326))
          dists[r] <- sum(geosphere::distVincentyEllipsoid(sf::st_coordinates(segment_wgs84)[,1:2]))
        }
      }
      transmitter_results[[t]] <- c(dists,NA)
    }
  }
  
  # assign transmitter IDs to complete results  
  names(transmitter_results) <- names(data_transmitter)
  
  # summary table
  max_hourly_dist <- unlist(lapply(transmitter_results, function(x) max(x, na.rm=T)))
  total_dist <- unlist(lapply(transmitter_results, function(x) sum(x, na.rm=T)))/1000
  
  max_ROM <- paste(round(mean(max_hourly_dist, na.rm=T)), "±", round(std.error(max_hourly_dist)))
  mean_total_dist <- paste(round(mean(total_dist, na.rm=T)), "±", round(std.error(total_dist)))
  
  number_ids <- length(which(!is.na(total_dist)))
  
  summary_results <- data.frame("species"=unique(data$species), "file"=unique(data$file), "analysed individuals"=number_ids,
                                "max ROM (m/h)"=max_ROM, "mean total distance (km)"=mean_total_dist, check.names=F)
  
  # mean stats (non-formatted)
  mean_values <- data.frame("Species"=unique(data$species), "max_ROM"=mean(max_hourly_dist, na.rm=T), "total_dist"=mean(total_dist, na.rm=T))
  
  
  return(list("summary_table"=summary_results, "mean_values"=mean_values, "complete_results"=transmitter_results))
}


#######################################################################################################
## Run function ###########################################################################################
#######################################################################################################

# calculate distances (km) and max rate of movement (mh-1)
start.time <- Sys.time()
distance_results <- lapply(data, calculateDistances, grid.resolution=50)
end.time <- Sys.time()
end.time - start.time
#save(distance_results, file="./Data/distance_results.RData")

# load results
#load("./Data/distance_results.RData")

# summary table (with formatted values ± standard error)
distance_summary <- do.call("rbind", lapply(distance_results, function(x) x$summary_table))
distance_summary <- distance_summary[order(distance_summary$species),]
rownames(distance_summary) <- NULL
#write.csv2(distance_summary, "./distance_summary_table.csv", row.names=F)

# mean values (by dataset)
dist_values <- do.call("rbind", lapply(distance_results, function(x) x$mean_values))
dist_values$File <- rownames(dist_values)
rownames(dist_values) <- NULL
complete_data <- join(complete_data, dist_values, by="File", type="left")

# complete results (by individual)
dist_values <- lapply(distance_results, function(x) x$complete_results)
total_dists <- lapply(dist_values, function(x) lapply(x, function(y) sum(y, na.rm=T)/1000))
total_dists <- lapply(total_dists, function(x) reshape2::melt(x))
total_dists <- mapply(function(table, dataset){table$File<-dataset; return(table)}, table=total_dists, dataset=names(distance_results), SIMPLIFY=F)
total_dists <- do.call("rbind", total_dists)
rownames(total_dists) <- NULL
colnames(total_dists)[1:2] <- c("Total_distance", "Transmitter")
complete_data <- join(complete_data, total_dists, by=c("Transmitter" ,"File"), type="left")
complete_data$Total_distance <- round(complete_data$Total_distance, 1)
complete_data$ROM_mh <- round(complete_data$Total_distance/complete_data$Monitored_hours, 1)
complete_data <- complete_data[,-which(colnames(complete_data)=="Total_distance")]
write.csv2(complete_data, "~/Desktop/metrics_by_individual.csv", row.names=F, fileEncoding="Windows-1252")


#######################################################################################################
## Prepare final data frame ###########################################################################
#######################################################################################################

response_cols <- c("IR", "k50_mean", "k95_mean", "mean_ROM", "max_ROM")
taxonomic_cols <- c("Species", "Order", "Family", "Genus")
methological_cols <- c("ArrayID", "NDetections", "SampleSize", "NReceivers", "ReceiverDensity", "SpatialScale_km")
trait_cols <- c("DemersPelag", "LongevityWild", "ReproMode", "Troph", "DepthMax", "MaxLengthTL", "Vulnerability")

final_data <- complete_data[, c("File", response_cols, taxonomic_cols, methological_cols, trait_cols)]
final_data$Species <- as.factor(final_data$Species)
final_data$Order <- as.factor(final_data$Order)
final_data$Family <- as.factor(final_data$Family)
final_data$Genus <- as.factor(final_data$Genus)
final_data$ArrayID <- as.factor(final_data$ArrayID)
final_data$DemersPelag <- as.factor(final_data$DemersPelag)
final_data$ReproMode <- as.factor(final_data$ReproMode)

levels(final_data$DemersPelag)[levels(final_data$DemersPelag)=="pelagic-neritic"] <- "pelagic"
final_data$DemersPelag <- factor(final_data$DemersPelag, levels=levels(final_data$DemersPelag)[c(4,2,1,3)])


#######################################################################################################
## Check experimental design biases ###################################################################
#######################################################################################################

final_data <- complete_data[, c("File", response_cols, taxonomic_cols, methological_cols, trait_cols)]
final_data <- final_data[order(final_data$File),]
final_data$ID <- 1:nrow(final_data)
final_data$NDetections <- final_data$NDetections/1000

xvars <- c("NReceivers", "NReceivers", "NReceivers", "SpatialScale_km")
yvars <- c("NDetections", "SpatialScale_km", "ReceiverDensity", "ReceiverDensity")
xlabs <- c(rep("Nº Receivers",3), "Spatial Scale (km)")
ylabs <- c("Nº Detections (x1000)", "Spatial Scale (km)", "Receiv Density (receiv/km2)", "Receiv Density (receiv/km2)")


pdf("~/Desktop/Experimental_variables.pdf", height=7, width=7, useDingbats=F)
m <- matrix(c(1,2,3,4,5,5), ncol=2, byrow=T)
layout(m)
par(mar=c(4,5,2,2))
for(i in 1:length(xvars)){
  xvar <- final_data[, xvars[i]]
  yvar <- final_data[, yvars[i]]
  plot(x=xvar, y=yvar, type="n", xlab=xlabs[i], ylab=ylabs[i], axes=F)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="gray96")
  mod <- lm(yvar ~ xvar)
  newx <- seq(par("usr")[1], par("usr")[2], length.out=100)
  preds <- predict(mod, newdata=data.frame(xvar=newx), interval='confidence', level=0.95)
  polygon(c(rev(newx), newx), c(rev(preds[,3]), preds[,2]), col='grey90', border = NA)
  abline(mod, col="blue", lty=2)
  cor <- cor.test(xvar, yvar)
  cor_string <- paste0("cor=", round(cor$estimate,2), " | p=", round(cor$p.value, 3))
  legend("topright", legend=cor_string, bty="n", cex=0.8)
  points(x=xvar, y=yvar, pch=21, col="blue", bg="white", cex=1.7, lwd=0.8)
  axis(1, at=pretty(xvar), labels=pretty(xvar), cex.axis=0.9)
  axis(2, at=pretty(yvar), labels=pretty(yvar), cex.axis=0.9, las=1)
  text(x=xvar, y=yvar, labels=final_data$ID, cex=0.6)
  if(cor$p.value<0.05){box(col="red")}else{box()}
}
plot.new()
ids <- paste(final_data$ID, "-", final_data$File)
legend("left", legend=ids, bty="n", ncol=3, cex=0.8)
dev.off()

#fit <- loess(yvar~xvar, span=1)
#xvals <- seq(min(xvar), max(xvar), (max(xvar) - min(xvar))/1000)
#lines(xvals, predict(fit, xvals), col='blue', lty=2)


#######################################################################################################
## Trophic level Plots ################################################################################
#######################################################################################################

xvars <- c("Troph", "Troph", "Troph", "Troph")
yvars <- c("k95_mean", "IR", "mean_ROM", "max_ROM")
xlabs <- c(rep("Trophic level",4))
ylabs <- c("KUD 95% (km2)", "Residency", "Average ROM (m/h)", "Max ROM (m/h)")

pdf("~/Desktop/Trophic_level.pdf", height=7, width=7, useDingbats=F)
m <- matrix(c(1,2,3,4,5,5), ncol=2, byrow=T)
layout(m)
par(mar=c(4,5,2,2))
for(i in 1:length(xvars)){
  xvar <- final_data[, xvars[i]]
  yvar <- final_data[, yvars[i]]
  plot(x=xvar, y=yvar, type="n", xlab=xlabs[i], ylab=ylabs[i], axes=F)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="gray96")
  mod <- lm(yvar ~ xvar)
  newx <- seq(par("usr")[1], par("usr")[2], length.out=100)
  preds <- predict(mod, newdata=data.frame(xvar=newx), interval='confidence', level=0.95)
  polygon(c(rev(newx), newx), c(rev(preds[,3]), preds[,2]), col='grey90', border = NA)
  abline(mod, col="blue", lty=2)
  cor <- cor.test(xvar, yvar)
  cor_string <- paste0("cor=", round(cor$estimate,2), " | p=", round(cor$p.value, 3))
  legend("topright", legend=cor_string, bty="n", cex=0.8)
  points(x=xvar, y=yvar, pch=21, col="blue", bg="white", cex=1.7, lwd=0.8)
  axis(1, at=pretty(xvar), labels=pretty(xvar), cex.axis=0.9)
  axis(2, at=pretty(yvar), labels=pretty(yvar), cex.axis=0.9, las=1)
  text(x=xvar, y=yvar, labels=final_data$ID, cex=0.6)
  if(cor$p.value<0.05){box(col="red")}else{box()}
}
plot.new()
ids <- paste(final_data$ID, "-", final_data$File)
legend("left", legend=ids, bty="n", ncol=3, cex=0.8)
dev.off()



#######################################################################################################
## Habitat ############################################################################################
#######################################################################################################

xvars <- c("DemersPelag", "DemersPelag", "DemersPelag", "DemersPelag")
yvars <- c("k95_mean", "IR", "mean_ROM", "max_ROM")
tlabs <- c("KUD 95% (km2)", "Residency", "Average ROM (m/h)", "Max ROM (m/h)")

pdf("~/Desktop/Habitat.pdf", height=7, width=9.5, useDingbats=F)
par(mar=c(4,5,2,2), mfrow=c(2,2))
for(i in 1:length(xvars)){
  xvar <- final_data[, xvars[i]]
  yvar <- final_data[, yvars[i]]
  b<- boxplot(yvar~xvar, axes=F, type="n", xlab="", ylab="", ylim=c(min(yvar)*0.9, max(yvar)*1.15))
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="gray96")
  boxplot(yvar~xvar, axes=F, add=T, col="gray80", ylim=c(min(yvar)*0.9, max(yvar)*1.15))
  title(main=tlabs[i], cex.main=0.9)
  axis(1, at=1:nlevels(xvar), labels=rep("", 4), cex.axis=0.75)
  mtext(side=1, at=1:nlevels(xvar), text=levels(xvar), line=0.8, cex=0.75 )
  axis(2, at=pretty(yvar), labels=pretty(yvar), cex.axis=0.9, las=1)
  stats <- table(xvar)
  offline <- (max(yvar)-min(yvar))/15
  text(x=1:nlevels(xvar), y=b$stats[5,]+offline, labels=paste0("n=",stats), cex=0.7)
  ktest <- kruskal.test(yvar~xvar)
  ktest_label <- paste0("X2=", sprintf("%.2f", ktest$statistic), "; p=", sprintf("%.3f", ktest$p.value))
  mtext(ktest_label, side=1, line=2.2, cex=0.65)
  if(ktest$p.value<0.05){box(col="red")
  }else{box()}
}
dev.off()



#######################################################################################################
## Run models #########################################################################################
#######################################################################################################

# scale and center numeric columns
final_data$NDetections  <- as.numeric(scale(final_data$NDetections))
final_data$NReceivers  <- as.numeric(scale(final_data$NReceivers))
final_data$SpatialScale_km <- as.numeric(scale(final_data$SpatialScale_km))
final_data$ReceiverDensity <- as.numeric(scale(final_data$ReceiverDensity))
final_data$LongevityWild <- as.numeric(scale(final_data$LongevityWild))
final_data$Troph <- as.numeric(scale(final_data$Troph))
final_data$DepthMax <- as.numeric(scale(final_data$DepthMax))
final_data$MaxLengthTL <- as.numeric(scale(final_data$MaxLengthTL))
final_data$Vulnerability <- as.numeric(scale(final_data$Vulnerability))

# check multicollinearity between predictors
corvif(final_data[,-c(1:10)])
numeric_cols <- which(unlist(lapply(final_data[,-c(1:8)], is.numeric)))
predictor_cors <- cor(final_data[,names(numeric_cols)], use="complete.obs")
predictor_cors[upper.tri(predictor_cors, diag=T)] <- NA
cor_matrix <- round(predictor_cors, 2)
cor_matrix

# check responses distributions
k95 <- final_data$k95_mean
hist(k95, 50)
hist(log(k95), 50)
qqnorm(log(k95))
qqline(log(k95))
descdist(k95, boot=1000)
descdist(log(k95), boot=1000)

ir <- final_data$IR
hist(ir,50)
qqnorm(ir)
qqline(ir)
descdist(ir, discrete=F)


# fit gamm




# fit models

model1 <- glmmTMB(log(k95_mean) ~ NReceivers +SampleSize +DemersPelag +LongevityWild +ReproMode +Troph +DepthMax 
                  +MaxLengthTL +Vulnerability + (1|ArrayID) + (1|Family/Genus/Species), data=final_data)

model2 <- glmmTMB(k95_mean ~ NReceivers +SampleSize +DemersPelag +LongevityWild +ReproMode +Troph +DepthMax 
                  +MaxLengthTL +Vulnerability + (1|ArrayID) + (1|Family/Genus/Species), data=final_data, family=Gamma)

model3 <- glmmTMB(k95_mean ~ NReceivers +SampleSize +DemersPelag +LongevityWild +ReproMode +Troph +DepthMax 
                  +MaxLengthTL +Vulnerability + (1|Family/Genus/Species), data=final_data, family=Gamma)

plot(fitted(selected_model), residuals(selected_model), type="n", main="Residuals vs Fitted Values", 
     xlab="Fitted", ylab="Residuals", las=1)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="grey97")
points(fitted(selected_model), residuals(selected_model))
abline(h=0, col="red")

selected_model <- model2
Anova(selected_model)
summary(selected_model)
plot(residuals(selected_model), fitted(selected_model))
dharma_res <- suppressMessages(simulateResiduals(selected_model))
plotQQunif(dharma_res)
plotResiduals(dharma_res)  
testDispersion(dharma_res)
testZeroInflation(dharma_res)



# effects plots
selected_model <- model2
pdf("~/Desktop/fixed_effects.pdf", height=14.5, width=14, useDingbats=F)
eff_result <- allEffects(selected_model)
for(e in 1:length(eff_result)){eff_result[[e]]$response <- "KUD 95%"}
plot(eff_result, type="response", colors=c("#25a0bb", "#ee6c11"), rows=4, cols=3)
dev.off()




#######################################################################################################
## Plot map with receivers and overall KUDs ###########################################################################################
#######################################################################################################


# import coastline layer (https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html)
coastline <- readOGR("~/Desktop/Marine Connectivity/Coastline/gshhg-shp-2.3.7/GSHHS_shp/i", layer="GSHHS_i_L1")
coastline <- spTransform(coastline, CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
coords <- do.call("rbind", lapply(data, function(x) x[,c("longitude","latitude")]))
coords <-  SpatialPoints(coords, proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
coastline <- crop(coastline, extent(coords)*1.1)

pdf("./Map.pdf", width=12, height=12, useDingbats=F )
plot(coastline, col="gray60", border=NA)
points(coords, pch=16, cex=0.6)
box()
dev.off()


#######################################################################################################
## Test influence of individual size in computed metrics ##############################################
#######################################################################################################

full_data <- do.call("rbind", lapply(results, function(x) x$complete_table))
rownames(full_data) <- NULL
colnames(full_data) <- tolower(colnames(full_data))
full_data <- join(full_data, fish_info[,c("file", "transmitter", "length_cm", "sex")], by=c("file", "transmitter"), type="left")

complete_k95s <- do.call("rbind", lapply(kud_results, function(x) data.frame("transmitter"=x$complete_KUDs$IDs,"K95"=x$complete_KUDs$K95)))
complete_k50s <- do.call("rbind", lapply(kud_results, function(x) data.frame("transmitter"=x$complete_KUDs$IDs,"K50"=x$complete_KUDs$K50)))
complete_k95s$file <- sub("\\..*", "", rownames(complete_k95s))
complete_k50s$file <- sub("\\..*", "", rownames(complete_k50s))
rownames(complete_k95s) <- NULL
rownames(complete_k50s) <- NULL
full_data <- join(full_data, complete_k95s, by=c("file", "transmitter"), type="left")
full_data <- join(full_data, complete_k50s, by=c("file", "transmitter"), type="left")


# standardize fish length within each species (normalized values between 0 and 1)
groupped_data <- split(full_data, f=full_data$file)
groupped_data <- lapply(groupped_data, function(x) {x$transmitter <- as.factor(x$transmitter); return(x)})
groupped_data <- lapply(groupped_data, function(x){x$transmitter<-droplevels(x$transmitter); return(x)})
groupped_data <- lapply(groupped_data, function(x) {x$length_scaled <- round(rescale(x$length_cm, c(0,1)),2); return(x)})
full_data <- do.call("rbind", groupped_data)

# only include datasets with > 10 fish 
groupped_data <- groupped_data[lapply(groupped_data, function(x) nlevels(x$transmitter))>10]
groupped_data <- groupped_data[lapply(groupped_data, function(x) length(which(!is.na(x$length_scaled))))>10]


pdf("~/Desktop/Ir_by_size.pdf", height=12, width=7, useDingbats=F)
par(mfrow=c(5,2), mar=c(4,4,2,2))
for(i in 1:length(groupped_data)) {
  cat(paste(names(groupped_data)[i], "\n"))
  data <- groupped_data[[i]]
  plot(y=data$ir, x=data$length_scaled, pch=16, cex=0.8, col="blue", cex.lab=1, ylim=c(-0.05,1.05), xlim=c(-0.05,1.05),
       xlab="", ylab="",main=names(groupped_data)[i], axes=F)
  title(xlab="Fish length (normalized)", line=2.2)
  title(ylab="IR", line=2.8)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="gray96")
  mod <- lm(ir~length_scaled, data=data)
  newx <- seq(par("usr")[1], par("usr")[2], length.out=100)
  preds <- predict(mod, newdata=data.frame(length_scaled=newx), interval='confidence', level=0.95)
  polygon(c(rev(newx), newx), c(rev(preds[,3]), preds[,2]), col='grey90', border=NA)
  abline(mod, col="blue", lty=2)
  cor <- cor.test(data$ir, data$length_scaled)
  cor_string <- paste0("cor=", round(cor$estimate,2), " | p=", round(cor$p.value, 3))
  legend("topright", legend=cor_string, bty="n", cex=0.8)
  points(x=data$length_scaled, y=data$ir, pch=21, col="blue", bg="white", cex=1, lwd=0.8)
  axis(1, at=seq(0,1,by=0.2), labels=seq(0,1,by=0.2), cex.axis=1)
  axis(2, at=seq(0,1,by=0.2), labels=seq(0,1,by=0.2), cex.axis=1, las=1)
  if(cor$p.value<0.05){box(col="red")}else{box()}
}
dev.off()



pdf("~/Desktop/KUD95_by_size.pdf", height=12, width=7, useDingbats=F)
par(mfrow=c(5,2), mar=c(4,4,2,2))
for(i in 1:length(groupped_data)) {
  cat(paste(names(groupped_data)[i], "\n"))
  data <- groupped_data[[i]]
  plot(y=data$K95, x=data$length_scaled, pch=16, cex=0.8, col="blue", cex.lab=1, xlim=c(-0.05,1.05),
       xlab="", ylab="", main=names(groupped_data)[i], axes=F)
  title(xlab="Fish length (normalized)", line=2.2)
  title(ylab="KUD 95%", line=3.1)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col="gray96")
  mod <- lm(K95~length_scaled, data=data)
  newx <- seq(par("usr")[1], par("usr")[2], length.out=100)
  preds <- predict(mod, newdata=data.frame(length_scaled=newx), interval='confidence', level=0.95)
  polygon(c(rev(newx), newx), c(rev(preds[,3]), preds[,2]), col='grey90', border=NA)
  abline(mod, col="blue", lty=2)
  cor <- cor.test(data$K95, data$length_scaled)
  cor_string <- paste0("cor=", round(cor$estimate,2), " | p=", round(cor$p.value, 3))
  legend("topright", legend=cor_string, bty="n", cex=0.8)
  points(x=data$length_scaled, y=data$K95, pch=21, col="blue", bg="white", cex=1, lwd=0.8)
  axis(1, at=seq(0,1,by=0.2), labels=seq(0,1,by=0.2), cex.axis=1)
  axis(2, at=pretty(data$K95), labels=sprintf("%.2f",pretty(data$K95)), cex.axis=1, las=1)
  if(cor$p.value<0.05){box(col="red")}else{box()}
}
dev.off()


