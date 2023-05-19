extract_photo_data <- function(folder, data.file = NULL, plot.map = TRUE, map.file = NULL){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # This function extracts file name, latitude, longitude, and date and time when the photo was taken from photos taken on a digital device. 
  #
  # If the device used to take the photos did not have GPS capability, then only the file name and data and time take will be extracted and a message will print in the console.
  #  Note: most smart phones record GPS coordinates when taking photos and a lot of newer digital cameras also have GPS capability.
  #
  # The function returns a tibble of the data, and, depending on how you set up the arguments, writes a .csv file of the data and an interactive leaflet map.
  # The interactive map contains points for each photo containing GPS coordinates. When you hover over the point, the image file name is displayed, and if you click on a point, 
  # the image appears in a popup window.  
  #
  #  Note: the rotation of the popup image is not always correct in the plot window; however, it appears to be correct in .html file. Not sure why or how to fix this issue.
  #
  # Inputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   folder - the path to where the photos are located. The folder must contain files that have ".jpg" extensions.
  #   data.file - the name of the output data file including .csv extension. If set to NULL, no file is written.
  #   plot.map - if TRUE, then an interactive map will appear in the plot window.
  #   map.file - the name of the output map widget file including .html extension. If set to NULL, then no file is written. 
  #              If you do not supply the full file path, the file will be saved to the supplied photos folder. 
  #              
  #
  #   **Note: `plot.map` and `map.file` work independently.**
  #
  # Outputs~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   This function returns a tibble containing the following variables: FileName, latitude, longitude, date, time.
  #   If plot.map = TRUE an interactive map appears in the plot window.
  #   If map.file is supplied, then an interactive map html is written.
  #
  # Example~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  # extract_photo_data(folder = "V:\\Library\\Media\\Photos\\2018 Trips\\Coho Baseline\\Glacier Cr", 
  #                       data.file = "PhotoData.csv", 
  #                       plot.map = TRUE, 
  #                       map.file = "PhotoData_map.html")
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
  setwd(folder)
  
  files <- list.files(folder, pattern = "*.jpg", full.names = FALSE, ignore.case = TRUE) #Get file names with .JPG extension
  
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
    
    readr::write_csv(output, file = data.file)
    
  }
  
  if(sum(!is.na(output$latitude)) == 0){
    
    stop("None of the files contained GPS data.")
    
  }
  
  if(plot.map == TRUE | !is.null(map.file)){
    
    images <- paste0(folder, "/", output$FileName) # These are the images to add as popups
    
    pts <- sf::st_as_sf(output %>% dplyr::filter(!is.na(latitude)), coords = c("longitude", "latitude"), crs = 4326) # Create sf object for plotting points
    
    map <- leaflet::leaflet() %>%
      leaflet::addProviderTiles(leaflet::providers$Esri.WorldTopoMap) %>% # Providers is a list from leaflet with all available base maps from different providers.
      leaflet::addCircleMarkers(data = pts, group = "points", color = "red", fillOpacity = 100, radius = 1.25, label = ~FileName) %>%
      leafpop::addPopupImages(images, group = "points", width = 400) %>% # This adds the images as popups. Image rotation is not always correct in the plot window, but appears to be correct in .html file. 
      leaflet::addScaleBar(position = "topright", options = leaflet::scaleBarOptions(metric = FALSE, imperial = TRUE)) %>% 
      leaflet::addMiniMap(position = "bottomright", 
                          tiles = leaflet::providers$Esri.WorldTopoMap, 
                          mapOptions = list())  
     }
      
  if(!is.null(map.file)){
      
      htmlwidgets::saveWidget(map, file = map.file)
    
  }
  
  if(plot.map){print(map)} 
  
  return(output)
  
}

