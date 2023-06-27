#' Plot Collections on a Map
#' 
#' This function plots collection locations on an interactive map.  
#' 
#' @param file the file path to output the html map with .html extension (default = NULL) 
#' @param input a `tibble` containing the following variables:
#'  \itemize{
#'      \item{number}{unique collection number; if NULL, this will default to 1:length(location); this will appear when you hover your mouse pointer over a waypoint}
#'      \item{location}{location name for each collection; this will pop up when you click on a waypoint}
#'      \item{group_name}{reporting group name for each collection}
#'      \item{group_color}{the group color for each collection (see details)}
#'      \item{latitude}{latitude in decimal degrees for each collection (see details)}
#'      \item{longitude}{longitude in decimal degrees for each collection (see details)}
#'  }
#' @param select.basemap if TRUE, select from a list of base maps in the console. If FALSE (default), basemap will be set to leaflet::providers$OpenStreetMap
#' @param png save map as png (default: FALSE)
#'  
#' @details
#'  If file is set to NULL, the map will not be save to an html file. Either R color names or hexidecimal colors can be used by this function. This function only works with coordinates in decimal degrees format with a negative sign to indicate western longitudes (e.g. 61.1596, -149.8876). 
#'  
#' @examples
#'  input <- tibble::tibble(number = 1:4,
#'                  location = c("Juneau Creek", "Russian River", "Killey River", "Funny River"),
#'                  group_name = c("Upper Kenai", "Upper Kenai", "Lower Kenai", "Lower Kenai"),
#'                  group_color = c("red", "red", "blue", "blue"),
#'                  latitude = c(60.49188, 60.4477, 60.478727, 60.42944),
#'                  longitude = c(-149.88454, -149.9865, -150.646294, -150.78141)
#'  )
#'  
#'  collections_map(input)
#'  
#' @export    
collections_map <- function(input, file = NULL, select.basemap = FALSE, png = FALSE){
  
  # Convert colors to hexadecimal - some R colors don't work in leaflet
  # Add 360 to longitudes < 0 so all points will show up on the map when centered on 180 degrees
  my.input <- input %>% 
    dplyr::mutate(group_color = GCLr::col2hex(group_color)$hex,
           longitude = dplyr::case_when(longitude < 0 ~ longitude + 360,
                                 TRUE~longitude))
  
  options <- names(leaflet::providers)
  
  mapname <- if(select.basemap == TRUE){select.list(options, "Select a basemap:", multiple = FALSE)}else{"OpenStreetMap"}
  
  minlong <- min(my.input$longitude)
  
  maxlong <- max(my.input$longitude)
  
  minlat <- min(my.input$latitude)
  
  maxlat <- max(my.input$latitude)
  
  map <- my.input %>% 
      leaflet::leaflet() %>%
      leaflet::addProviderTiles(leaflet::providers[[mapname]], options =  leaflet::leafletOptions(worldCopyJump = FALSE)) %>% 
      leaflet::addCircleMarkers(lng = my.input$longitude, lat = my.input$latitude, color = ~group_color, popup = ~location, fillOpacity = 100, radius = 2, label = ~number) %>% 
      leaflet::addScaleBar(position = "bottomleft",options = leaflet::scaleBarOptions(imperial = FALSE)) %>% 
      leaflet::addMiniMap( position = "bottomright") %>% 
      leaflet::addLegend("topright", colors = my.input$group_color %>% unique(), labels = my.input$group_name %>% unique,
                opacity = 1) %>% 
    leaflet::setMaxBounds(lng1 = minlong, lat1 = minlat, lng2 = maxlong, lat2 = maxlat)
  
  if(!is.null(file)){htmlwidgets::saveWidget(map, file = file)}
  
  if(!is.null(file) & png == TRUE){
    
    file <- gsub(".html", ".png", x = file)
    
    mapview::mapshot(map, file = file)
    
  }
    
  map
  
}