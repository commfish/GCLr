#' Convert GENEPOP File to ".gcl" Objects
#'
#' Create ".gcl" and `LocusControl` objects from a GENEPOP file.
#'
#' @param filename The full file path of the GENEPOP file.
#' 
#' @details 
#' An allele conversion argument may be added in the future if there is a need (mainly for SNPs). This would require a conversion input file from the user to define the corresponding alleles for each locus (e.g., 1 = A, 2 = C, 3 = G, 4 = T).
#' 
#' @return Assigns ".gcl" and LocusControl objects to the current workspace and returns a vector of the ".gcl" objects created without the ".gcl" extension (i.e., `sillyvec`).
#' 
#' @examples
#' GCLr::genepop2gcl(filename = system.file("genepop", "ex_genepop.txt", package = "GCLr"))
#' 
#' @export
genepop2gcl <- function(filename){

  if(strsplit(filename, ".gen")==filename){
    
    stop("Genepop file needs a '.gen' extention, hoser!!!")
    
    }

  if(exists("LocusControl")){
    
    stop("LocusControl already exists! Run this function in a workspace where no LocusControl exists.")
    
  }

  rawdat <- scan(filename, what = "", sep = "\n")

  len <- length(rawdat)

  popind <- sapply(rawdat, function(lin){
    
    strsplit(lin, " ")[[1]][1]=="Pop"
    
    }) %>% as.vector()

  npops <- sum(popind)
 
  ORD <- order(popind, decreasing = TRUE)[1:npops]

  loci <- sapply(rawdat[2:(ORD[1]-1)], function(str){
    
    strsplit(str, " ")[[1]][1]
    
    })

  nloc <- length(loci)
  
  # Attribute variables
  attr <- 
    c(
      "FK_FISH_ID",
      "COLLECTION_ID",
      "SILLY_CODE",
      "PLATE_ID",
      "PK_TISSUE_TYPE",
      "CAPTURE_LOCATION",
      "CAPTURE_DATE",
      "END_CAPTURE_DATE",
      "MESH_SIZE",
      "MESH_SIZE_COMMENT",
      "LATITUDE",
      "LONGITUDE",
      "AGENCY",
      "VIAL_BARCODE",
      "DNA_TRAY_CODE",
      "DNA_TRAY_WELL_CODE",
      "DNA_TRAY_WELL_POS",
      "CONTAINER_ARRAY_TYPE_ID",
      "SillySource"
    )
  
  # locus variables
  loc_vars <- c(loci, paste0(loci, ".1")) %>% 
    as.character() %>% 
    sort()
  
  # Put data in .gcl object tibble format
  # create initial data tibble
  dat0 <- rawdat[-(1:(nloc+1))][rawdat[-(1:(nloc+1))]!="Pop"] %>%
    tibble::as_tibble() %>% 
    tidyr::separate(col = value, sep = " ,  ", into = c("SillySource", "geno")) %>% 
    tidyr::separate(col = geno, sep = " ", into = loci) 
  
  if(length(grep(x = dat0$SillySource, pattern = "_[0-9]*$")) == length(dat0$SillySource)){
    
    dat0 <- dat0 %>% 
      tidyr::separate(col = SillySource, sep = "_", into = c("SILLY_CODE", "FK_FISH_ID"), remove = FALSE) 
    
  }else{
    
    dat0 <- dat0 %>%
      dplyr::rename(SILLY_CODE = SillySource) %>% 
      dplyr::group_by(SILLY_CODE) %>% 
      dplyr::mutate(FK_FISH_ID = 1:length(SILLY_CODE)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(SillySource = paste0(SILLY_CODE, "_", FK_FISH_ID))
    
    message("New FK_FISH_ID's were created because none were supplied in the GENEPOP file.")

  }
  
 dat0 <- dat0 %>% 
    dplyr::mutate(COLLECTION_ID = NA_integer_, 
           PLATE_ID = NA_character_, 
           PK_TISSUE_TYPE = NA_character_, 
           CAPTURE_LOCATION = NA_character_, 
           CAPTURE_DATE = as.Date(NA_character_), 
           END_CAPTURE_DATE = as.Date(NA_character_), 
           MESH_SIZE = NA_character_, 
           MESH_SIZE_COMMENT = NA_character_,
           LATITUDE = NA_real_,
           LONGITUDE = NA_real_,
           AGENCY = NA_character_,
           VIAL_BARCODE = NA_character_,
           DNA_TRAY_CODE = NA_character_,
           DNA_TRAY_WELL_CODE = NA_integer_,
           DNA_TRAY_WELL_POS = NA_character_,
           CONTAINER_ARRAY_TYPE_ID = NA_integer_
           )
  
  # Separate genotypes into two variables
  nchar <- dat0 %>% 
    dplyr::select(dplyr::all_of(loci)) %>% 
    dplyr::summarize(dplyr::across(dplyr::everything(), ~base::nchar(.) %>% max())) %>% 
    max()
  
  for(locus in loci){
    
    dat0 <- dat0 %>%
      tidyr::separate(!!dplyr::sym(locus), into = c(locus, paste0(locus, ".1")), sep = nchar/2, remove = FALSE) %>% 
      dplyr::mutate(dplyr::across(dplyr::all_of(c(locus, paste0(locus, ".1"))), ~as.numeric(.))) %>% 
      dplyr::mutate(dplyr::across(dplyr::all_of(c(locus, paste0(locus, ".1"))), ~as.character(.)))
    
  }
  
  # Final data tibble
  dat <- dat0 %>% 
    dplyr::mutate(dplyr::across(dplyr::all_of(loc_vars), ~dplyr::na_if(x = ., y = "0"))) %>% 
    dplyr::select(dplyr::all_of(c(attr, loc_vars)))
  
  # Assign silly objects to workspace
  sillyvec <- dat$SILLY_CODE %>% 
    unique()
  
  for(silly in sillyvec){
    
    my.dat <- dat %>% dplyr::filter(SILLY_CODE == silly)
    
    assign(x = paste0(silly, ".gcl"), value = my.dat, pos = 1)
    
  }
  
  message(paste0("A total of ", length(sillyvec), " '.gcl' objects were created from the GENEPOP file and assgined to the current workspace."))
  
  # Make LocusControl
  # alleles
  alleles <- lapply(loci, function(locus){
    
    tibble::tibble(call = dat %>% 
             dplyr::select(dplyr::starts_with(locus)) %>% 
             unlist() %>% 
             unique() %>% 
             na.omit() %>% 
             as.numeric() %>% 
             sort() %>% 
             as.character()) %>% 
      dplyr::mutate(allele = seq(length(call))) %>% 
      dplyr::select(allele, call)
  }) %>% purrr::set_names(loci) 
  
  # nalleles
  nalleles <- lapply(loci, function(locus){
    
    dim(alleles[[locus]])[[1]]
    
  }) %>% 
    unlist %>% 
    purrr::set_names(loci)
  
  # ploidy
  ploidy <-  dplyr::if_else (
    condition = dat0[ ,paste0(loci, ".1")] %>%is.na() %>% apply(., 2, sum)  == 0, 
    true = 2, 
    false = 1) %>% purrr::set_names(loci)
  
  assign("LocusControl", tibble::tibble(MarkerSuite = "From GENEPOP", 
                                        locusnames = loci,
                                        Publishedlocusnames = NA_character_, 
                                        nalleles = nalleles, 
                                        ploidy = ploidy,
                                        alleles = alleles), pos = 1)
  
  message("LocusControl created from GENEPOP file")
  
  return(sillyvec)
    
}