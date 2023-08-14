#' Extract Geotag Data from Photos
#'
#' This function extracts file name, latitude, longitude, and date and time from photos taken on a digital device with GPS capability. If the device used to take the photos did not have GPS capability, then only the file name and data and time take will be extracted and a message will print in the console.
#' The interactive map contains points for each photo containing GPS coordinates. When you hover over the point, the image file name is displayed, and if you click on a point, the image appears in a popup window.  
#'
#' @param folder The path to where the photos are located. The folder must contain files that have ".jpg" extensions.
#' @param data.file The name of the output data file including .csv extension. If set to NULL, no file is written.
#' @param plot.map If TRUE, then an interactive map will appear in the plot window (see details).
#' @param map.file The name of the output map widget file including .html extension. If set to NULL, then no file is written. 
#'                 If you do not supply the full file path, the file will be saved to the supplied photos folder (see details). 
#'
#' @details Most smart phones record GPS coordinates when taking photos and a lot of newer digital cameras also have GPS capability. However, most devices have to be set up to record this information. For instance, on an iPhone you need to have location services turned on or it won't record the location of a photo. 
#'          Another thing to note is that `plot.map` and `map.file` arguments work independently, so you can produce a map file without plotting in R or plot in R and don't produce a file.
#'
#' @return A tibble containing the following variables: FileName, latitude, longitude, date, time.
#'
#' @examples
#'  
#' \dontrun{
#' extract_photo_data(folder = "V:\\Library\\Media\\Photos\\2018 Trips\\Coho Baseline\\Glacier Cr",
#'                    data.file = "PhotoData.csv",
#'                    plot.map = TRUE,
#'                    map.file = "PhotoData_map.html")
#' 
#' }
##' @export
extract_photo_data <- function(folder, data.file = NULL, plot.map = TRUE, map.file = NULL){
  
  files <- list.files(folder, pattern = "*.jpg", full.names = TRUE, ignore.case = TRUE) #Get file names with .JPG extension
  
  data0 <- lapply(files, function(file){
    
    my.dat <- exifr::read_exif(path = file,  recursive = FALSE) #Recursive doesn't work so looping through files instead
    
    attach(my.dat)
    
    existsGPS <- exists("GPSPosition")
      
    detach()
    
    if(existsGPS){
    
      my.dat %>%  
        dplyr::select(FileName, GPSPosition, DateTimeOriginal)
      
    }else{
      
      cat(paste0("File '", file, "' does not contain GPS coordinates.\n"))
      
      my.dat %>%  
        dplyr::select(FileName, DateTimeOriginal) %>% 
        dplyr::mutate(GPSPosition = NA)

    }
    
  }) 
  
  data <- data0 %>% 
    dplyr::bind_rows()
  
  output <- data %>% 
    dplyr::select(FileName, GPSPosition, DateTimeOriginal) %>% 
    tidyr::separate(GPSPosition, into = c("latitude", "longitude"), sep = " ") %>% 
    dplyr::mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude)) %>% 
    tidyr::separate(DateTimeOriginal, into = c("date", "time"), sep = " ") %>%
    dplyr::mutate(date = gsub(":", replacement = "-", x = date) %>% as.Date())
   
  if(!is.null(data.file)){
    
    readr::write_csv(output, file = paste0(folder, "/", data.file))
    
  }
  
  if(sum(!is.na(output$latitude)) == 0){
    
    stop("None of the files contained GPS data.")
    
  }
  
  if(plot.map == TRUE | !is.null(map.file)){
    
    images <- paste0(folder, "/", output %>% dplyr::filter(!is.na(latitude)) %>% pull(FileName)) # These are the images to add as popups
    
    pts <- sf::st_as_sf(output %>% dplyr::filter(!is.na(latitude)), coords = c("longitude", "latitude"), crs = 4326) # Create sf object for plotting points
    
    map <- leaflet::leaflet() %>%
      leaflet::addProviderTiles(leaflet::providers$Esri.WorldTopoMap) %>% # Providers is a list from leaflet with all available base maps from different providers.
      leaflet::addCircleMarkers(data = pts, group = "points", color = "red", fillOpacity = 100, radius = 1.25, label = ~FileName, clusterOptions = leaflet::markerClusterOptions(removeOutsideVisibleBounds = F, spiderfyDistanceMultiplier=1.5)) %>%
      leafpop::addPopupImages(images, group = "points", width = 400) %>% # This adds the images as popups. Image rotation is not always correct in the plot window, but appears to be correct in .html file. 
      leaflet::addScaleBar(position = "topright", options = leaflet::scaleBarOptions(metric = FALSE, imperial = TRUE)) %>% 
      leaflet::addMiniMap(position = "bottomright", 
                          tiles = leaflet::providers$Esri.WorldTopoMap, 
                          mapOptions = list())  
     }
      
  if(!is.null(map.file)){
      
      htmlwidgets::saveWidget(map, file = paste0(folder, "/", map.file))
    
  }
  
  if(plot.map){print(map)} 
  
  invisible(output)
  
}