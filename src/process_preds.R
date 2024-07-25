#' process_preds
#'
#' This function reads in a list of spatial predictor file paths and processes them for a niche
#' model. Specifically, it reads in raw files, sets their crs and extent to a define 'template'
#' file, and rasterizes the file if necessary.
#' 
#' @param input_path File pathway to the raw predictor spatial file
#' @param template_path File pathway to the 'template' file that will be used to standardize crs and ext
#' @param resolution Desired resolution for processed raster file output
#' @param save Whether to save (TRUE) the resulting dataframe (as .tif) or not (FALSE)
#' @param output_path If `save = TRUE`, the file path to save the dataframe.
#'
#' @return A spatraster with crs and ext that is consistent with template and res that is set to desired resolution

process_file <- function(input_path, template_path, resolution, save = FALSE, output_path = "data/input_processed") {
  # Read in the template .shp and save as an sf object
  template <- read_sf(template_path)
  assign("template_sf", template, envir = .GlobalEnv)
  
  # Initialize variables to hold the processed objects
  rast_shapefile <- NULL
  raster_file <- NULL
  
  if (grepl("\\.shp$", input_path, ignore.case = TRUE)) {
    # Read in shapefile
    shapefile <- st_read(input_path)
    
    # Set CRS to match the template
    shapefile <- st_transform(shapefile, crs = crs(template))
    
    # Set extent to match the template (crop if necessary)
    shapefile <- st_crop(shapefile, st_bbox(template))
    
    ## Rasterize the shapefile
    # Create an empty raster template with the specified resolution
    raster_template <- rast(ext(template), resolution = resolution)
    # Rasterize shapefile 
    rast_shapefile <- terra::rasterize(shapefile, raster_template, fun = mean)
    # Set CRS to match the shapefile
    crs(rast_shapefile) <- st_crs(shapefile)$proj4string
    
    # Extract filename without extension for naming
    file_name <- tools::file_path_sans_ext(basename(input_path))
    
    # Save as an independent object
    assign(file_name, rast_shapefile, envir = .GlobalEnv)
    
  } else if (grepl("\\.tif$", input_path, ignore.case = TRUE)) {
    # Read in raster and template files
    raster_file <- rast(input_path)
    
    # Set template temporarily to raster CRS
    shapefile_temp <- st_transform(template, crs = crs(raster_file))
    raster_file <- crop(raster_file, shapefile_temp)
    
    # Project cropped raster into template CRS
    raster_file <- project(raster_file, crs(template), res = resolution)
    
    # Extract filename without extension for naming
    file_name <- tools::file_path_sans_ext(basename(input_path))
    
    # Save as an independent object
    assign(file_name, raster_file, envir = .GlobalEnv)
    
  } else {
    stop("Unsupported file type.")
  }
  
  # Save the processed objects to the specified directory if 'save' is TRUE
  if (save) {
    if (!dir.exists(output_path)) {
      dir.create(output_path, recursive = TRUE)
    }
    
    if (!is.null(rast_shapefile)) {
      # Write raster file
      writeRaster(rast_shapefile, filename = file.path(output_path, paste0(file_name, ".tif")), overwrite = TRUE)
    }
    
    if (!is.null(raster_file)) {
      # Write raster file
      writeRaster(raster_file, filename = file.path(output_path, paste0(file_name, ".tif")), overwrite = TRUE)
      
    }
  }
}