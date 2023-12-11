#' Get Sample Inventory Plots
#' 
#' This function pulls a sample inventory from Loki and produces 3 plots.
#' 
#' @param year the years you want to look at. Years must be greater than 2016 or the function will give an error (see details).
#' 
#' @param dir the directory where you want the plot .png files to be saved. **Note: the files are named automatically, so just the path to the folder here.**
#' 
#' @param username your state user name
#' 
#' @param password your password used to access LOKI; see Eric Lardizabal if you need to set up a password
#'
#' @details
#'  This function requires a build time stamp for for each collection to work. Build time stamps were not added to Loki until 2017, so this function will only work for years >= 2017.
#'  This function creates three plots:
#'  
#' \describe{
#' 
#'   \item{Plot 1}{A barplot of the number of samples in Loki by species and year}
#'   \item{Plot 2}{A barplot of the overall number of samples in Loki by species for the years supplied}
#'   \item{Plot 3}{A heatmap of the number of sample by month and year, with the overall average by month for the years supplied}
#'   
#' }
#' 
#' @examples
#' \dontrun{
#' 
#' get_fish_inventory(years = 2018:2023, username = "myusername", password = "mypassword")
#' 
#' }
#' 
#' @export
get_fish_inventory <- function(years, dir, username, password){
  
  if(min(years) < 2017){
    
    stop("Loki does not contain build time stamp information for collections prior to 2017. Change your years arguement and try again.")
    
  }
  
  options(java.parameters = "-Xmx4g")
  
  url <- GCLr:::loki_url() #This is a function that gets the correct URL to access the database on the oracle cloud
  
  drvpath <- system.file("java", "ojdbc8.jar", package = "GCLr")
  
  drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver", classPath = drvpath, " ")
  
  con <- RJDBC::dbConnect(drv, url = url, user = username, password = password)
  
  gnoqry <- paste("SELECT * FROM AKFINADM.FISH_INVENTORY_V") 
  
  dataAll0 <- RJDBC::dbGetQuery(con, gnoqry)  #Pulling data from LOKI using the connection and genotype query
  
  discon <- RJDBC::dbDisconnect(con)
  
  dataAll <- dataAll0 %>% 
    dplyr::filter(!is.na(BUILD_TIMESTAMP)) %>% 
    dplyr::mutate(Date = lubridate::as_date(BUILD_TIMESTAMP)) %>% 
    tidyr::separate(Date, into = c("year", "month", "day"), remove = FALSE) %>% 
    dplyr::mutate(dplyr::across( c("year", "month", "day"), as.numeric)) %>% 
    tibble::as_tibble() %>% 
    dplyr::filter(year %in% years)
    
  #For naming files and plot titles
  
  if(length(years)== 1){
    
    year_range <- years
    
  }else{
    
    year_range <- paste0(min(years), "-", max(years))
 
  }
  
  #For adjusting plot width
  width <- 4.5+length(years)*0.5
   
  #Number samples by species and year
  plot1_dat <- dataAll %>%
    dplyr::mutate(COMMON_NAME = dplyr::case_when(!COMMON_NAME %in% c("Salmon, Chinook", "Salmon, Sockeye", "Salmon, Coho", "Salmon, Pink", "Salmon, Chum") ~ "Other",
                                   TRUE~COMMON_NAME)) %>% 
    dplyr::group_by(year, COMMON_NAME) %>% 
    dplyr::summarize(N = sum(SAMPLE_COUNT)) %>% 
    dplyr::mutate(species = gsub(pattern = "Salmon, ", replacement = "", x = COMMON_NAME))
      
  cols <- c(Chinook = "blue", Chum = "green", Coho = "gray", Pink = "pink", Sockeye = "red", Other = "black")
  
 plot1 <-  plot1_dat %>% 
   dplyr::mutate(species = factor(species, levels = names(cols))) %>% 
   ggplot2::ggplot(ggplot2::aes(y = N, x = year, fill = species)) + 
   ggplot2::geom_bar(stat = "identity", position = "dodge")+
   ggplot2::scale_fill_manual(name = "Species", values = cols)+
   ggplot2::ylab("Number of Samples")+
   suppressWarnings(ggplot2::scale_x_discrete(name ="Year", limits = years)) +
   ggplot2::ggtitle(label = "Number of Samples Built in Loki by Species and Year") +
   ggplot2::theme(axis.title = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 12))
 
 ggplot2::ggsave(plot1, filename = paste0(dir, "/Number samples by species and year for ", year_range, ".png"), device = "png", height = 5, width = width)
 
 plot2 <-  plot1_dat %>%
   dplyr::group_by(species) %>% 
   dplyr::summarise(N = sum(N)) %>% 
   dplyr::mutate(species = factor(species, levels = names(cols))) %>% 
   ggplot2::ggplot(ggplot2::aes(y = N/100000, x = species, fill = species)) + 
   ggplot2::geom_bar(stat = "identity", position = "dodge")+
   ggplot2:: scale_fill_manual(name = "Species", values = cols)+
   ggplot2::ylab("Number of Samples (100,000's)")+
   ggplot2::ggtitle(label = paste0("Overall Number of Samples Built in Loki by Species: ", year_range))+
   ggplot2::theme(axis.title = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 12))
 
 ggplot2::ggsave(plot2, filename = paste0(dir, "/Overall number samples by species ", year_range, ".png"), device = "png", height = 5, width = 6)
 
 
 #Samples by year and month
 plot3_dat <- dataAll %>%
   dplyr::group_by(year, month) %>% 
   dplyr::summarize(N = sum(SAMPLE_COUNT), .groups = "drop") %>% 
   dplyr::mutate(month = factor(month.name[month], levels = month.name), year = as.character(year))
 
 if(length(years)>1){
   
   average <- plot3_dat %>% 
     dplyr::group_by(month) %>% 
     dplyr::summarize(N = mean(N)) %>% 
     dplyr::mutate(year = "Average")
   
   plot3_dat <- bind_rows(plot3_dat, average)
   
 }
 
 plot3 <-  plot3_dat %>% 
   dplyr::mutate(month = factor(month, levels = month.name[12:1]), year = factor(year, levels = unique(year))) %>% 
   ggplot2::ggplot(ggplot2::aes(y = month, x = year, fill = N)) +
   ggplot2::geom_tile()+
   ggplot2::scale_fill_gradient(low = "green", high = "red", name = "Number of Samples")+
   ggplot2::geom_text(ggplot2::aes(y = month, x = year, label = round(N, digits = 0)), size = 2.5)+
   ggplot2::ylab("Month")+
   ggplot2::xlab("Year")+
   ggtitle(label = "Number of Samples Built in Loki by Month and Year") +
   ggplot2::theme(axis.title = ggplot2::element_text(size = 14), axis.text = ggplot2::element_text(size = 12), axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) 
 
  ggplot2::ggsave(plot3, filename = paste0(dir, "/Number samples by month for ", year_range, ".png"), device = "png", height = 5, width = width)
 
  print(paste0("Your plots have been saved to: ", dir))
  
}