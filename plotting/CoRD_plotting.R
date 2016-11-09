# ===============================================================================
# Copyright 2016 Daniel Cadol
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# ===============================================================================

#install.packages("rgdal") #needs GDAL installed on machine already
#install.packages("tiff")
#install.packages("sp")
#install.packages("raster")
#library(rgdal) #thought I would need this, but haven't yet
library(sp)
library(raster)


############################
####     Functions      ####
############################

### Function to create polygons that bound a user-defined number of segments
seg <- function(upstream_endpoint,downstream_endpoint,segments){
	valley_endpts <- cbind(c(upstream_endpoint[1],downstream_endpoint[1]),
		c(upstream_endpoint[2],downstream_endpoint[2]))
	valley_segpts <- cbind(seq(valley_endpts[1,1],valley_endpts[2,1],length.out=(segments+1)),
		seq(valley_endpts[1,2],valley_endpts[2,2],length.out=(segments+1)))

	### Create endpoints of segment slicing lines
	#find valley line slop
	valley_azmuth <- (valley_endpts[1,2]-valley_endpts[2,2])/(valley_endpts[1,1]-valley_endpts[2,1])
	#slope of slicing line must be the negative of the inverse of the valley line slope
	seg_azmuth <- -1/valley_azmuth
	#define how far to the side of the valley line to slice
	half_valley_width <- 500 #meters
	#Function to find slice endpoints given a center line point
	seg_slice_lines <- function (cen_x,cen_y){
		east_x <- cen_x+half_valley_width*cos(atan(seg_azmuth))
		east_y <- cen_y+half_valley_width*cos(atan(seg_azmuth))*seg_azmuth
		west_x <- cen_x-half_valley_width*cos(atan(seg_azmuth))
		west_y <- cen_y-half_valley_width*cos(atan(seg_azmuth))*seg_azmuth
		slice_endpts <- cbind(c(east_x, west_x),c(east_y, west_y))
		return(slice_endpts)
	}
	# apply the slice line generating function to each segmenting point
	slices_lines <- mapply(seg_slice_lines, valley_segpts[,1], valley_segpts[,2],SIMPLIFY=T)

	### Create polygons that outline each valley segment
	#function to sort out the points into the right order
	make_poly <- function(seg){
		poly_x <- c(slices_lines[1:2,seg],slices_lines[2:1,seg+1])
		poly_y <- c(slices_lines[3:4,seg],slices_lines[4:3,seg+1])
		#polygon(poly_x, poly_y) #just for plotting
		poly_pts <- cbind(poly_x, poly_y)
		return(poly_pts)
	}
	#apply sorting function to create list of sorted points for each segment
	seg_poly_pts <- lapply(1:segments, make_poly)
	# Create a list containing SpatialPolygon type data objects (from library(sp))
	seg_polygons <- lapply (seg_poly_pts,Polygon)
	return(seg_polygons)
}

### Function to break a single raster into the user-defined number of segments created in the seg function
one_ras_segmentation <- function(vegmap, seg_polygons)	{
	
	### Create a list of raster objects, each one a segment of the valley
	# function to mask out a raster segment
	clip_by_mask <- function(seg_poly){
		polygon_obj <- spPolygons(seg_poly)
		vegmap_seg <- mask(vegmap, polygon_obj)
		return(vegmap_seg)
	}
	# apply masking function to each segment polygon
	vegmap_segments <- lapply(seg_polygons, clip_by_mask)
	vegmap_segments_comb <- stack(vegmap_segments)
	return(vegmap_segments_comb)
}

### function to break a multi-raster stack into a user-defined number of segments
multi_ras_segmentation <- function(vegmap_stack,seg_polygons)	{
	### Create a list of RasterBrick objects, each one a segment of the valley
	# function to mask out a raster segment
	clip_by_mask <- function(seg_poly){
		polygon_obj <- spPolygons(seg_poly)
		vegmap_stack_seg <- mask(vegmap_stack, polygon_obj)
		return(vegmap_stack_seg)
	}
	# apply masking function to each segment polygon
	vegmap_stack_segments <- lapply(seg_polygons, clip_by_mask)
	return(vegmap_stack_segments)
}
############################
####  End of Functions  ####
############################

####################################
###  Multiple time step rasters  ###

### Inputs
	# raster for each year of the model run
		vegmap_yr1 <- raster("Documents/ResearchProjects/VirtualWatershed/CoRD_Plotting/vegout_clip_utm.tif")
		vegmap_yr2 <- raster("Documents/ResearchProjects/VirtualWatershed/CoRD_Plotting/vegout_clip_utm.tif")
		vegmap_yr3 <- raster("Documents/ResearchProjects/VirtualWatershed/CoRD_Plotting/vegout_clip_utm.tif")
		vegmap_yr4 <- raster("Documents/ResearchProjects/VirtualWatershed/CoRD_Plotting/vegout_clip_utm.tif")
		vegmap_yr5 <- raster("Documents/ResearchProjects/VirtualWatershed/CoRD_Plotting/vegout_clip_utm.tif")
			# etc (and with different maps for each year in a real analysis!)
	#peak flows for each year in cms
		flow_hist <- c(50,15,200,50,50)
	# user defined valley line endpoints
		upstream_endpoint <- c(344575,3954305)
		downstream_endpoint <- c(344070,3953500)
	# user defined number of segments
		segments <- 10
	#folder where plots should be sent
		plot_path <- "Documents/ResearchProjects/VirtualWatershed/CoRD_Plotting/Plots"

### Stack annual maps into a raster brick
	vegmaps_all <- list(vegmap_yr1, vegmap_yr2, vegmap_yr3, vegmap_yr4, vegmap_yr5)
	vegmap_stack <- brick(vegmaps_all)
	vegmap_stack <- trim(vegmap_stack)
	
### Run multi-raster segmentation !! ###
	seg_polygons <- seg(upstream_endpoint,downstream_endpoint,segments)
	vegmap_stack_segments <- multi_ras_segmentation(vegmap_stack, seg_polygons)

### Plot Maps
tiff(filename=paste(plot_path,"VegMapsByYear.tif",sep="/"), width=7, height=7, units="in", res=150, pointsize=10,compression="lzw")
	plot(vegmap_stack)
dev.off()
	#plot(vegmap_stack[[1]])
	# seg <- 6
	# year <- 1
	# plot(vegmap_stack_segments[[seg]][[year]])

### Find the number of cells in each unique vegetation class through time for the WHOLE STUDY REACH
	# find all unique classes (assume all classes are present in the first year)
	classes_present <- unique(vegmap_yr1)
	# convert raster to x,y,value data frame
	vegmap_pts_by_year <- lapply(as.list(vegmap_stack),rasterToPoints)
	for (year in 1:length(vegmaps_all)){
		#function to count number of cells (rows in the data frame) in a given class
		count_cells <- function(class) {
			sum(vegmap_pts_by_year[[year]][,3]==class)
		}
		#apply counting function to all unique classes
		cell_count <- sapply(classes_present, count_cells)
		if(year==1) {cell_counts_all <- data.frame(cell_count)}
		else {cell_counts_all <- cbind(cell_counts_all,cell_count)}
		names(cell_counts_all)[year] <- paste("year",year,sep="_")
	}

### Find the number of cells in each unique vegetation class through time for EACH SEGMENT
	# find all unique classes (assume all classes are present in the first year)
	classes_present <- unique(vegmap_yr1)
	# convert raster to x,y,value data frame
	vegmap_stack_pts_by_year <- lapply(as.list(vegmap_stack_segments),rasterToPoints)
	# have a list pre-set
	seg_cell_counts <- vector("list",segments)
	for (seg in 1:segments){
		for(year in 1:length(vegmaps_all)){
			#function to count number of cells (rows in the data frame) in a given class
			count_cells <- function(class) { 
				sum(vegmap_stack_pts_by_year[[seg]][,year+2]==class)
			}
			#apply counting function to all unique classes
			cell_count <- sapply(classes_present, count_cells)
			if(year==1) {cell_counts_by_yr <- data.frame(cell_count)}
			else {cell_counts_by_yr <- cbind(cell_counts_by_yr,cell_count)}
			names(cell_counts_by_yr)[year] <- paste("year",year,sep="_")
		}
		seg_cell_counts[[seg]] <-cell_counts_by_yr
	}

# Descriptive stats of the class counts (e.g., mean of all years, max of all years, etc)
seg_cell_stats <- vector("list",segments) 
for(i in 1:segments){
	seg_cell_stats[[i]] <- data.frame(cbind(apply(seg_cell_counts[[i]],1,mean),apply(seg_cell_counts[[i]],1,min),apply(seg_cell_counts[[i]],1,max),apply(seg_cell_counts[[i]],1,sd)))
	colnames(seg_cell_stats[[i]]) <- c("mean","min","max","st_dev")
}

# create a color ramp for use with the various classes
num_of_classes <- length(classes_present)
color_set<-colorRampPalette(c('blue','purple','red','orange','gold'))
classcols <- color_set(num_of_classes)[as.numeric(cut(1:num_of_classes,breaks = num_of_classes))]

# create a color ramp for use with the various segments
color_set<-colorRampPalette(c('gold','green','cyan','blue'))
segcols <- color_set(segments)[as.numeric(cut(1:segments,breaks = segments))]

### plot change in each class's cell count through time
tiff(filename=paste(plot_path,"ClassCellCountByYear.tif",sep="/"), width=7, height=7, units="in", res=150, pointsize=10,compression="lzw")
	par(fig=c(0,1,0.3,1), mar=c(2,4,1,1), new=FALSE)
	class <- 1
	plot(1:length(vegmaps_all),cell_counts_all[class,],ylim=c(0,max(cell_counts_all)),xlab="Year",ylab="Full Map Cell Count",col=classcols[class],type='o')
	for (class in 2:length(classes_present)){
		points(1:length(vegmaps_all),cell_counts_all[class,],col=classcols[class],type='o')
	}
	legend(1,150000, classes_present, pch=1, col=classcols, ncol=3)

	par(fig=c(0,1,0,0.3), mar=c(4,4,1,1), new=TRUE)
	plot(1:length(flow_hist),flow_hist,xlab="Year",ylab="Peak Flow (cms)",col="blue",type='o')
dev.off()

### plot each class's cell statistics as they vary down the reach
seg_cell_stats_array <- array(unlist(seg_cell_stats),c(length(classes_present),4,segments))

tiff(filename=paste(plot_path,"ClassStatsBySeg.tif",sep="/"), width=7, height=7, units="in", res=150, pointsize=10,compression="lzw")
	class<-1
	plot(1:segments,seg_cell_stats_array[class,1,],col=classcols[class],type='l',lwd=2,ylim=c(0,max(seg_cell_stats_array)),xlab="Segment Number",ylab="Number of Cells in Class")
	polygon(c(1:segments,rev(1:segments)),c(seg_cell_stats_array[class,2,],rev(seg_cell_stats_array[class,3,])), col=paste0(classcols[class],"40"),border="transparent" )
	for (class in 2:length(classes_present)){
		points(1:segments,seg_cell_stats_array[class,1,],col=classcols[class],type='l',lwd=2)
		polygon(c(1:segments,rev(1:segments)),c(seg_cell_stats_array[class,2,],rev(seg_cell_stats_array[class,3,])), 	col=paste0(classcols[class],"40"),border="transparent")
	}
	legend("topleft", legend=classes_present, lwd=2, col=classcols, ncol=3)
dev.off()


### Plot Variation in Moran's I through time
moran_by_year <- vector("numeric",length(vegmaps_all))
for (year in 1:length(vegmaps_all)) {
	moran_by_year[year] <- Moran(vegmap_stack[[year]])
}

seg_morans_by_year <- vector("list", segments)
for (seg in 1:segments){
	single_seg_morans <- vector("numeric",length(vegmaps_all))
	for (year in 1:length(vegmaps_all)) {
		single_seg_morans[year] <- Moran(vegmap_stack_segments[[seg]][[year]])
	}
	seg_morans_by_year[[seg]] <- single_seg_morans
}

tiff(filename=paste(plot_path,"MoranIByTime.tif",sep="/"), width=7, height=6, units="in", res=150, pointsize=10,compression="lzw")
	plot(1:length(vegmaps_all),moran_by_year,type='o',xlab="Year of simulation",ylab="Moran's I",ylim=c(min(unlist(seg_morans_by_year)),max(unlist(seg_morans_by_year))),lwd=2)
	for (seg in 1:segments){
		points(1:length(vegmaps_all),seg_morans_by_year[[seg]],col=segcols[seg],type='o')
	}
	legend("bottomright",legend=c("Full Map", paste("seg",1:segments)),col=c("#000000",segcols),pch=1,lwd=c(2,rep(1,segments)))
dev.off()

####################################
###    Analyze a single raster   ###
raster_path <- "Documents/ResearchProjects/VirtualWatershed/CoRD_Plotting/vegout_clip_utm.tif"
vegmap <- raster(raster_path)
vegmap <-trim(vegmap)

### user defined valley line endpoints
upstream_endpoint <- c(344575,3954305)
downstream_endpoint <- c(344070,3953500)
### user defined number of segments
segments <- 10

### Run single raster segmentation !! ###
seg_polygons <- seg(upstream_endpoint,downstream_endpoint,segments)
vegmap_segments_comb <- one_ras_segmentation(vegmap,seg_polygons)


### Plot raster, valley line, and segmenting points
plot(vegmap)
points(valley_endpts[,1],valley_endpts[,2],type='l')
points(valley_segpts[,1],valley_segpts[,2])

plot(vegmap_segments_comb)

### plot historgram of all unique values in the raster
barplot(vegmap,las=3)

### Calculate and plot Moran's I
	# Whole map
	map_moran <- Moran(vegmap)
	# By Segment
	for(i in 1:segments){
		ith_moran <- Moran(vegmap_segments[[i]])
		if (i==1){seg_morans <- ith_moran}
		else {seg_morans <- append(seg_morans, ith_moran)}
	}
	plot(1:segments,seg_morans,type='o',pch=16)
	abline(h=map_moran,lwd=3)

### Plot downstream change in # of cells in each class
	classes_present <- unique(vegmap)
	class_count_by_seg <- data.frame(array(rep(classes_present, segments+1),c(length(classes_present),segments+1)))
	#function to count number of cells (rows in the data frame) in a given class
	count_cells <- function(class) { sum(vegmap_seg_pts[,3]==class) }
	for(i in 1:segments){
		vegmap_seg_pts <- rasterToPoints(vegmap_segments[[i]])
		cell_count <- sapply(classes_present, count_cells)
		class_count_by_seg[,i+1] <- cell_count
		colnames(class_count_by_seg)[i+1] <- paste("seg",i,sep="_")
	}
# cumulative version of class count, for plotting purposes
class_cum_by_seg <- cumsum(class_count_by_seg[,1:segments+1])

# create a color ramp for use with the various classes
num_of_classes <- length(classes_present)
color_set<-colorRampPalette(c('blue','purple','red','orange','gold'))
classcols <- color_set(num_of_classes)[as.numeric(cut(1:num_of_classes,breaks = num_of_classes))]

i<-1
plot(1:segments, class_cum_by_seg[i,1:segments],type='o',col=classcols[i],ylim=c(0,max(class_cum_by_seg)),xlab="Segment",ylab="Cell Count")
for(i in 2:num_of_classes){
	points(1:segments, class_cum_by_seg[i,1:segments],type='o',col=classcols[i])
}





