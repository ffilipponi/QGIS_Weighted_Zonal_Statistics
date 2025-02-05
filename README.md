# QGIS Weighted Zonal Statistics
### QGIS processing tool to perform zonal statistics for raster data over vector polygons, weighted by fractional pixel area

Tool can be easily used in QGIS through the Processing Toolbox. It is a R script that makes use of 'terra' and 'exactextractr' to compute zonal statistics for raster data over vector polygons, weighted by fractional pixel area.

To install the processing tool in QGIS, first install [R cran](https://cran.r-project.org) and QGIS 'Processing R Provider' plugin. Then install the following required R cran libraries: 'sf', 'raster', 'terra', 'exactextractr' (>= 0.9.0). The following R command line should do the job:

`install.packages(pkgs = "exactextractr", dependencies = TRUE, repos = "http://cran.us.r-project.org")`

Configure R 'User library folder' in QGIS Processing Settings, in order to have the installed R libraries available in QGIS. Finally copy files available in [rscripts folder](https://github.com/ffilipponi/QGIS_Weighted_Zonal_Statistics/tree/main/rscripts) inside folder set in QGIS Processing Settings under 'R scripts folder'.
