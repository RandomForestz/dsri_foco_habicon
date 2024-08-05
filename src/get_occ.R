#' Retrieve Occurrence Data
#'
#' This function is a wrapper around `spocc::occ()` that retrieves and cleans occurrences for a set of species.
#' Note: depending on how many occurrences are returned this step could take a while.
#' 
#' @param species Vector of species scientific names
#' @param aoi An sf object of the area to search occurrences for 
#' @param dates A vector of min and max date to search for, in YYYY-MM-DD format
#' @param limit The maximum number of occurrences to return for a single species. Default is 100000.
#' @param type  Either 'sf' to return an sf object or 'df' to return a dataframe. Default is a dataframe.
#' @param save Whether to save (TRUE) the resulting dataframe or not (FALSE)
#' @param output_path If `save = TRUE`, the file path to save the dataframe or sf object.
#'
#' @return  Either an sf object or tabular dataframe of all cleaned species occurrences.
get_occ <- function(species,
                    aoi,
                    dates = c("2014-01-01", "2024-12-31"),
                    limit = 100000,
                    type = "df",
                    save = FALSE,
                    output_path = "data/input_occ") {
  
  ## get bounds for spocc search
  bounds <- aoi %>% 
    # make sure in WGS 84
    st_transform(crs = 4326) %>% 
    # get bounds 
    st_bbox() %>% 
    as.vector()
  
  # can't get occ() filters to work so will do that after the original call
  sp_occ <- spocc::occ(
    query = species,
    from = c(
      "gbif",
      "inat",
      "ecoengine",
      "vertnet",
      "bison",
      "ala",
      "idigbio",
      "ebird"
    ),
    geometry = bounds,
    date = c("2014-01-01", "2024-12-31"), # doesn't take date as a user argument, has to be hardcoded. Will have to filter later
    has_coords = TRUE,
    limit = limit
  ) %>%
    spocc::occ2df()
  
  
  # clean the dataset
  sp_occ_clean <- sp_occ %>%
    # clean up species names between providers
    mutate(species = stringr::word(name, 1, 2)) %>%
    # filter date range
    #filter(date >= dates[1] & date <= dates[2]) %>% 
    # remove duplicates
    distinct(species, longitude, latitude, .keep_all = TRUE) %>% 
    # remove missing coordinates
    filter(!is.na(latitude), !is.na(longitude)) %>%
    # make spatial and filter to aoi (vertnet skips that part)
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    # convert to crs of aoi
    st_transform(crs = crs(aoi)) %>%
    st_filter(aoi)
  
  
  # return object based on 'type'
  if(type == "df") {
    
    final_occ <- sp_occ_clean %>% 
      # convert back to coordinates
      st_transform(crs = 4326) %>% 
      mutate(lat = st_coordinates(.)[, "Y"],
             long = st_coordinates(.)[, "X"]) %>% 
      st_drop_geometry()
      
      
  } else {
    
    final_occ <- sp_occ_clean %>% 
      # convert back to decimal degrees
      st_transform(crs = 4326)
    
  }
  
 # save the file
  if (save) {
    
    if (!dir.exists(output_path)) {
      dir.create(output_path, recursive = TRUE)
    }
    
    if (inherits(final_occ, "sf")) {
      # Write shapefile
      st_write(final_occ,
               file.path(output_path, "retrieved_occurrences.shp"),
               append = FALSE)
      
    } else {
      # save dataframe
      write_csv(final_occ,
                file.path(output_path, "retrieved_occurrences.csv"))
      
    }

  }
  
  return(final_occ)
  
}