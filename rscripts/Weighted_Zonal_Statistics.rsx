##Polygons=file
##Raster=file
##Names=string NULL
##Output=file
# ##Output=output vector
##Format=selection GeoPackage;ESRI_Shapefile
##Statistics=string NULL
##Method=selection exactextractr;terra;fast;rough

# #!/usr/bin/env Rscript
# Title         : weighted_zonal_statistics.R
# Description   : Perform zonal statistics of raster data over vector polygons, weighted by fractional pixel area
# Date          : Nov 2023
# Version       : 2.1
# Copyright name: CC BY-NC-SA
# License       : GPL v3
# Authors       : Federico Filipponi
# Maintainer    : Federico Filipponi <federico.filipponi@gmail.com>
# ########################################################

# # ### for debug
# Polygons <- "/media/DATA/ISPRA/new7/Incendi_e_Biodiversita_NBFC/sites/monte_Morrone/incendio_2017/Burned_area_2017.gpkg"
# Raster <- "/media/DATA/ISPRA/new7/CSA/SIMonA/pubblicazione_prodotti/ECM/ecm-f4-2020_v1.0/ECM-F4-2020_v1.0/ECM-F4_2020_v1.0.tif"
# Names <- "A,B,C,D,E"
# Names <- NULL
# Output <- "/media/DATA/download/test_zonal/test01.gpkg"
# Format <- 0
# Statistics <- "mean,sum,min,max"
# Method <- 1
# Method <- 0

# ###################################
# Initialize process
message("-------------------------------------")
message(paste(c("Started at: "), Sys.time(), sep=""))
message("Initializing process ...")

ptm <- proc.time()

# get options
# out_vfmt <- "GPKG"
input_vector <- normalizePath(path=as.character(Polygons), winslash = "/", mustWork = FALSE)
input_raster <- normalizePath(path=as.character(Raster), winslash = "/", mustWork = FALSE)
if(Names == "NULL"){
  band_names_in <- "NULL"
} else {
  band_names_in <- as.character(Names)
}
output_file <- normalizePath(path=Output, winslash = "/", mustWork = FALSE)
vformat <- c("GeoPackage", "ESRI_Shapefile")
out_vformat <- as.character(vformat[Format + 1])

method_list <- c("exactextractr","terra","fast","rough")
extract_method <- as.character(method_list[as.integer(Method)+1])
if(sum(as.integer(extract_method %in% c("terra","exactextractr","fast","rough"))) == 0){
  message("ERROR: argument 'Method' shoud be set to one of the followings: 'terra', 'exactextractr', 'fast', 'rough'.")
  q("no")
}

if(Statistics == "NULL"){
  message("ERROR: At least one statistic should be requested with parameter 'Statistics'.")
  q("no")
}
extract_stats <- unlist(strsplit(x = Statistics, split=","))
if("area" %in% extract_stats){
  categorical <- TRUE
} else {
  categorical <- FALSE
}
if(!categorical){
  if(sum(!(extract_stats %in% c("mean","sum","min","max","median","mode","variance","stdev"))) > 0){
    message("ERROR: At least one of the requested statistics is not allowed.")
    q("no")
  } else {
    if(sum(as.integer(extract_method %in% c("terra","fast"))) > 0){
      if(sum(!(extract_stats %in% c("mean","sum","min","max"))) > 0){
        message(paste("ERROR: At least one of the requested statistics is not supported for analysis with method '", extract_method, "'." , sep=""))
        message("ERROR: Run again using supported methods: 'extractextractr' (statistics weighted by fractional pixel area), 'rough'.")
        q("no")
      }
    }
  }
}

# check environment
invisible(tryCatch(find.package("terra"), error=function(e){
  message("ERROR: R package 'terra' not installed.")
  q("no")
  }))
if(extract_method %in% c("terra","fast","rough")){
  terra_version <- as.character(packageVersion("terra"))
  terra_version <- as.double(paste(unlist(strsplit(terra_version, split='\\.'))[1], ".", unlist(strsplit(terra_version, split='\\.'))[2], sep=""))
  if(terra_version < 1.7){
    message("ERROR: R package 'terra' version should be >= '1.7.0'.")
    q("no")
  }
}
if(extract_method == "exactextractr"){
  invisible(tryCatch(find.package("exactextractr"), error=function(e){
    message("ERROR: R package 'exactextractr' not installed. Consider to use 'terra' package as extraction 'Method'.")
    q("no")
  }))
  # check exactextractr version
  exactextractr_version <- as.character(packageVersion("exactextractr"))
  exactextractr_version <- as.double(paste(unlist(strsplit(exactextractr_version, split='\\.'))[1], ".", unlist(strsplit(exactextractr_version, split='\\.'))[2], sep=""))
  if(exactextractr_version < 0.9){
    stop("R package 'exactextractr' version should be >= '0.9.0'.")
  }
}

# check arguments
if(!file.exists(input_raster)){
  stop(paste("Input raster file '", input_raster, "' not found", sep=""))
}
if(!file.exists(input_vector)){
  stop(paste("Input vector file '", input_vector, "' not found", sep=""))
}
if(file.exists(output_file)){
  stop(paste("Output file '", output_file, "' already exists", sep=""))
}

# create output folder
if(!dir.exists(dirname(output_file))){
  dir.create(path=dirname(output_file), recursive = TRUE, showWarnings = TRUE)
}

# load required libraries
suppressWarnings(suppressPackageStartupMessages(require(sf, quietly=TRUE)))
suppressWarnings(suppressPackageStartupMessages(require(terra, quietly=TRUE)))
terraOptions(progress=0)
if(extract_method == "exactextractr"){
  suppressWarnings(suppressPackageStartupMessages(require(exactextractr, quietly=TRUE)))
}

# set environment
if(extract_method %in% c("terra")){
  tex_fast <- FALSE
  tex_exact <- TRUE
}
if(extract_method %in% c("fast")){
  tex_fast <- TRUE
  tex_exact <- FALSE
}
if(extract_method %in% c("rough")){
  tex_fast <- FALSE
  tex_exact <- FALSE
}

# ###################################
# Import data
message("-------------------------------------")
message("Importing data ...")

# read vector data
# deal with 'Z' dimension error (need to drop 'Z')
vecto <- sf::st_zm(sf::st_read(input_vector, stringsAsFactors=FALSE, quiet=TRUE))
# subset only polygons and multipolygons features
vecto <- vecto[sf::st_is(vecto, c("POLYGON", "MULTIPOLYGON")),]

if(nrow(vecto) == 0){
  stop("No polygons found in input vector layer")
}

# read raster data
r <- terra::rast(input_raster)
# set band names
rnl <- terra::nlyr(r)

band_names <- unlist(strsplit(x = band_names_in, split=","))
if(length(band_names) != rnl | band_names == "NULL"){
  if(!categorical){
    warning("Provided band names is not the same length of input raster band number. Continue with standard band names.")
  }
  names(r) <- paste("BAND_", sprintf(fmt="%003i", c(1:rnl)), sep="")
} else {
  if(out_vformat == "ESRI_Shapefile"){
    if(sum(nchar(band_names) > 10)){
      warning("At least one band name exceed 10 digits allowed by 'ESRI Shapefile' output file format.")
    }
  }
  names(r) <- band_names
}

if(categorical){
  # check if input file has 1 band
  if(terra::nlyr(r) > 1){
    stop("Supported raster file for 'categorical' processing should be a 1-band raster.")
    q("no")
  }
  # check if data is of supported datatype
  rdt <- terra::datatype(r)
  if(sum(as.numeric(rdt %in% c("INT1U", "INT2U"))) == 0){
    stop("Supported raster file for 'categorical' processing should be of type 'INT1U' or 'INT2U'.")
    q("no")
  }
}

if(sum(as.integer(!(sf::st_is_valid(vecto)))) > 0){
  message("Invalid geometries found in input vector data. Attempt to fix topological errors ...")
  vecto <- tryCatch(sf::st_make_valid(vecto), error=function(e){ 
    message("ERROR: Cannot fix input vecotr geometries.")
    q("no")
  })
}

# reproject vector data if necessary
vecto_epsg <- sf::st_crs(vecto)$epsg
raster_epsg <- as.integer(terra::crs(r, describe=TRUE)$code)
if(vecto_epsg != raster_epsg){
  vecto <- sf::st_transform(vecto, raster_epsg)
}

# ###################################
# query raster over spatial polygons
message("-------------------------------------")
message("Querying raster file ...")

if(exists("vecto")){
  
  if(categorical){
    
    # set terra RAM usage
    terra::terraOptions(datatype="INT1U")
    terra::terraOptions(todisk=TRUE)
    
    # generate set of values to query
    classV <- sort(as.integer(unlist(terra::unique(r, na.rm=TRUE))))
    ll <- length(classV)
    
    if(ll > 256){
      message("Warning: Number of classes to process for area calculation is greater than 256.")
    }
    
    # cycle on categorical raster classes
    for(u in 1:ll){
      
      message(paste("Processing class ", u, "/", ll, sep=""))
      
      v <- classV[u]
      rV <- terra::which.lyr(r == v)
      # rV <- terra::rast(r, vals=0)
      # rV[which(r[] == v)] <- 1
      invisible(gc())
      
      # query raster on spatial polygons
      if(extract_method == "exactextractr"){
        qut <- exactextractr::exact_extract(x = rV, y = vecto, weights="area", default_weight=0, fun="sum", progress=FALSE)
      }
      if(extract_method == "terra"){
        qut <- tryCatch(suppressWarnings(terra::extract(rV, vecto, weights=tex_fast, exact=tex_exact, cells=FALSE, ID=FALSE, fun="sum", na.rm=TRUE)), error=function(e){
          message("Error when using highest accuracy for raster zonal statistics. Going on with lower accuracy (force 'Fast' option).")
          suppressWarnings(terra::extract(rV, vecto, weights=tex_exact, exact=tex_fast, cells=FALSE, ID=FALSE, fun="sum", na.rm=TRUE))
        })
      }
      if(extract_method %in% c("fast","rough")){
        qut <- suppressWarnings(terra::extract(rV, vecto, weights=tex_fast, exact=tex_exact, cells=FALSE, ID=FALSE, fun="sum"))
      }
      
      # convert statistics to area
      qut <- data.frame("LAST"=(qut * (terra::xres(rV) * terra::yres(rV))))
      # assign names based on statistics
      if(length(band_names) == ll){
        names(qut) <- band_names[u]
      } else {
        names(qut) <- paste("C", "_", v, sep="")
      }
      # bind to vector data
      vecto <- cbind(vecto, qut)
      
      # remove temporary raster and clean workspace
      suppressWarnings(rm(list=c("rV","qut")))
      terra::tmpFiles(current=TRUE, orphan=TRUE, old=FALSE, remove=TRUE)
      invisible(gc())
      
    }
    
  } else {
    
    # cycle in requested statistics
    for(statis in extract_stats){
      
      s <- which(extract_stats == statis)
      message(paste("Processing statistic ", s, "/", length(extract_stats), sep=""))
      
      # query raster on spatial polygons
      if(extract_method == "exactextractr"){
        qut <- exactextractr::exact_extract(x = r, y = vecto, weights="area", default_weight=0, fun=statis, progress=FALSE)
        qut <- data.frame(qut)
      }
      if(extract_method == "terra"){
        qut <- tryCatch(suppressWarnings(terra::extract(r, vecto, weights=tex_fast, exact=tex_exact, cells=FALSE, ID=FALSE, fun=statis, na.rm=TRUE)), error=function(e){
          message("Error when using highest accuracy for raster zonal statistics. Going on with lower accuracy (force 'Fast' option).")
          suppressWarnings(terra::extract(r, vecto, weights=tex_exact, exact=tex_fast, cells=FALSE, ID=FALSE, fun=statis, na.rm=TRUE))
        })
      }
      if(extract_method %in% c("fast","rough")){
      statis_ext <- statis
        if(extract_method == "rough"){
          if(statis == "stdev"){
            statis_ext <- "sd"
          }
          if(statis == "variance"){
            statis_ext <- "var"
          }
        }
        qut <- suppressWarnings(terra::extract(r, vecto, weights=tex_fast, exact=tex_exact, cells=FALSE, ID=FALSE, fun=statis_ext, na.rm=TRUE))
      }
      
      statis_n <- statis
      if(statis == "stdev"){
        statis_n <- "std"
      }
      if(statis == "variance"){
        statis_n <- "var"
      }
      if(statis == "median"){
        statis_n <- "med"
      }
      
      # assign names based on statistics
      names(qut) <- paste(names(r), "_", as.character(statis_n), sep="")
      # bind to vector data
      vecto <- cbind(vecto, qut)
    }
  }
}

# reproject back data if necessary
if(vecto_epsg != raster_epsg){
  vecto <- sf::st_transform(vecto, vecto_epsg)
}

# ###################################
# Export result as vector file
message("-------------------------------------")
message("Writing output file ...")

if(out_vformat == "ESRI_Shapefile"){
  out_vfmt <- "ESRI Shapefile"
}
if(out_vformat == "GeoPackage"){
  out_vfmt <- "GPKG"
}

sf::st_write(vecto, dsn=output_file, driver=out_vfmt, quiet=TRUE)

# export results to QGIS
# Output <- sat
message("-------------------------------------")
message("Done")
message("-------------------------------------")
message(paste(c("Ended at: "), Sys.time(), sep=""))
message(paste("Elapsed time ", as.character(paste(as.integer(as.numeric(proc.time() - ptm)[3]/3600), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/3600 - as.integer(as.numeric(proc.time() - ptm)[3]/3600)) * 60)), ":", sprintf("%02i", as.integer((as.numeric(proc.time() - ptm)[3]/60 - as.integer(as.numeric(proc.time() - ptm)[3]/60)) * 60)), sep="")), " hours", sep=""))

# clean workspace and exit
rm(list=ls())
invisible(gc())
# q("no")
