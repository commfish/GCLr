#' @title Plot Allele Frequency & Fis for SNPs
#'
#' @description 
#' This function plots ("*.pdf") and returns allele frequencies, Fis values, and HWE p-values for each locus.
#'
#' @param sillyvec A character vector of silly codes without the ".gcl" extension, objects must exist in your environment.
#' @param loci A character vector of locus names. Note, only loci with ploidy = 2 and nalleles = 2 will be plotted.
#' @param groupvec An optional numeric vector indicating the group affiliation of each silly in `sillyvec`, 
#' must be the same length as `sillyvec` (default = `NULL`). If `NULL`, `groupvec` will be `1:length(sillyvec)`
#' @param alpha A numeric vector of length 1 specifying the critical HWE p-value for plotting.
#' @param group_col An optional character vector of colors corresponding to each group in `groupvec` (default = `NULL`).
#' @param file A character vector of length 1 with the full file path including ".pdf" extension where plots will be saved.
#' @param group.pch A vector of point shape `pch` values corresponding to each group with `length = max(groupvec)`. 
#' If only 1 `pch` is supplied it will be recycled for all groups (default `pch` = 19).
#' @param point.size A numeric vector of length 1 specifying the size of points on the plot, used by [ggplot2::geom_point()] 
#' (default = 1).
#' @param pval.cex A numeric vector of length 1 specifying the character expansion (i.e. size of text) for HWE p-values (default = 3).
#' @param pval.digits A numeric vector of length 1 specifying the number of digits to show for HWE p-values (default = 2).
#' @param line.width A numeric vector of length 1 specifying the width of the line connecting the points, used by [ggplot2::geom_line()] 
#' default = 0.5).
#' @param ncores A numeric value for the number of cores to use in a \pkg{foreach} `%dopar%` loop (default = 4). 
#' If the number of cores exceeds the number on your device ([parallel::detectCores()]), then all cores will be used.
#' @param LocusCtl an object created by [GCLr::create_locuscontrol()] (default = LocusControl).
#'
#' @returns A tibble of HWE p-values, Fis values, and allele frequencies and a (".pdf") of locus-specific plots:
#'     \itemize{
#'       \item \code{output}: a data.frame with 7 columns containing the full HWE output:
#'         \itemize{
#'           \item \code{population}: silly code
#'           \item \code{locus}: locus name
#'           \item \code{pval}: HWE p-value as calculated by [HardyWeinberg::HWExactMat()]
#'           \item \code{fis}: Fis estimate as calculated by [HardyWeinberg::HWf()]
#'           \item \code{allele freq}: allele frequency (proportion) of the 1st allele in `LocusCtl` as calculated by [GCLr::calc_freq_pop()]
#'           }
#'       \item \code{file}: a (".pdf") of stacked allele frequency and Fis plots for all sillys in `sillyvec` 
#'                          with one page per locus in `loci`
#'       }
#' 
#' @details
#' HWE p-values are printed if less than alpha. Loci with a ploidy of 1 (e.g., mitochondrial SNPs or microhaplotypes) are excluded.
#' The allele frequency is that of the 1st allele in `LocusCtl` (i.e., alphabetical alleles).
#' 
#' @seealso 
#' [HardyWeinberg::HardyWeinberg-package()]
#' [HardyWeinberg::HWExactMat()]
#' [HardyWeinberg::HWf()]
#' [GCLr::calc_freq_pop()]
#' [GCLr::plot_allele_freq()]
#'
#' @examples
#' sillys <- GCLr::base2gcl(GCLr::ex_baseline)
#' loci <- GCLr::ex_LocusControl$locusnames[-c(97,98)]
#' GCLr::plot_freq_fis_4snps(sillyvec = sillys, loci = loci, LocusCtl = GCLr::ex_LocusControl)
#' 
#' @export
plot_freq_fis_4snps <-
  function(sillyvec,
           loci,
           groupvec = NULL,
           alpha = 0.05,
           groupcol = NULL,
           file = NULL,
           group.pch = 21,
           point.size = 3,
           pval.cex = 3,
           pval.digits = 2,
           line.width = 0.7,
           ncores = 4,
           LocusCtl = LocusControl) {

  start.time <- Sys.time()

  if(is.null(file)){

    file <- paste0(getwd(),"/FreqFisPlot.pdf")

  }
  
  if(ncores > parallel::detectCores()) {ncores = parallel::detectCores()}

  ploidy <- LocusCtl$ploidy[loci]

  nalleles <- LocusCtl$nalleles[loci]
  
  if(length(setdiff(loci, loci[ploidy==2 & nalleles==2]))>=1){
    
    message(paste0("The following loci were excluded from the plots becuase they have a ploidy of 1 or more than 2 alleles: ", paste0(setdiff(loci, loci[ploidy==2 & nalleles==2]), collapse = ", ")))
    
    }
  
  loci <- loci[ploidy==2 & nalleles==2]

  get_counts <- function(my.gcl, loci){
    
    lapply(loci, function(locus){
      
      my.alleles <- LocusCtl$alleles[[locus]] %>% 
        pull(call)
      
      my.gcl %>% 
        dplyr::select(contains(locus)) %>% 
        dplyr::mutate_all(.funs = factor, levels = my.alleles) %>% 
        dplyr::mutate_all(.funs = as.numeric) %>% 
        dplyr::mutate_all(.funs = as.character) %>%
        dplyr::mutate(counts = rowSums(.=="1")) %>% 
        dplyr::pull(counts)
    
    }) %>% 
      dplyr::bind_cols() %>% 
      purrr::set_names(loci)

  }
  
  gcls <- lapply(sillyvec, function(silly){get(paste0(silly,".gcl"), pos = 1)}) %>% 
    purrr::set_names(sillyvec)

  cl <- parallel::makePSOCKcluster(ncores)
  
  doParallel::registerDoParallel(cl, cores = ncores)  
  
  `%dopar%` <- foreach::`%dopar%`
  
  #Start parallel loop
  HWE_df <- foreach::foreach(silly = sillyvec, .packages = c("tidyverse", "HardyWeinberg")) %dopar% {  # , .export = "LocusCtl"

    my.gcl <- gcls[[silly]]

    genecounts0 <- 3 - get_counts(my.gcl, loci) 

    genecounts <- t(apply(genecounts0, 2, tabulate, nbins = 3))

    HWEpval <- suppressWarnings(HardyWeinberg::HWExactMat(genecounts, verbose = FALSE))$pvalvec

    FixBOOL <- apply(genecounts, 1, function(gcounts){sum(gcounts==0)>1})

    HWEfis <- suppressWarnings(base::apply(genecounts[!FixBOOL, ], 1, HardyWeinberg::HWf)) %>% 
      tibble::as_tibble() %>% 
      dplyr::rename(fis = value) %>% 
      dplyr::mutate(locus = dimnames(genecounts[!FixBOOL, ])[[1]])
   
    tibble::tibble(pop_no = match(silly, sillyvec) %>% na.omit() %>% as.numeric(), population = silly, locus = loci, pval = HWEpval) %>% 
      dplyr::left_join(HWEfis, by = "locus") %>% 
      dplyr::mutate(fis = tidyr::replace_na(fis, 0))
    
  } %>% 
    dplyr::bind_rows() %>% 
    dplyr::mutate(population = factor(population, levels = sillyvec))#silly
  
  parallel::stopCluster(cl) #End parallel loop
  
  if (is.null(groupvec)) {
    
    groupvec <- seq_along(sillyvec)
    
  }

  G <- dplyr::n_distinct(groupvec)

  if (is.null(groupcol)) {
    
    groupcol <- rainbow(G)
    
  }
  
  if (is.numeric(groupcol)) {
    
    groupcol <- colors()[groupcol]
    
  }  

  PopCol <- groupcol[groupvec] 
  
  Freq <- suppressMessages(GCLr::calc_freq_pop(sillyvec = sillyvec, loci = loci, ncores = ncores, LocusCtl = LocusCtl))
    
  q <- Freq %>% 
    dplyr::filter(allele_no == 1) %>% 
    dplyr::pull(proportion)
  
  HWE_df$q <- q
  
  if(length(group.pch)==1){

    group.pch <- rep(group.pch, G)

  }
  
  PopPch <- group.pch[groupvec]
  
  fisylim0 <- max(abs(HWE_df$fis), na.rm = TRUE)

  fisylim <- c(-fisylim0, fisylim0)

  pdf(file, width = 11, height = 8.5, family = "Helvetica", pointsize = 12)

  for(locus in loci){
  
    my.dat <- HWE_df %>% 
      dplyr::filter(locus==!!locus)
    
    freqplot <- my.dat %>%
      ggplot2::ggplot(ggplot2::aes(y = q, x = pop_no)) +
      ggplot2::geom_hline(yintercept = c(0, 1), size = 0.5) +
      ggplot2::geom_line(color = "black",
                         linetype = "dashed",
                         linewidth = line.width) +
      ggplot2::geom_point(
        fill = PopCol,
        color = "black",
        shape = PopPch,
        size = point.size
      ) +
      ggplot2::ylim(0, 1) +
      ggplot2::ylab("Frequency") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank()
      ) +
      ggplot2::ggtitle(label = unique(locus))
    
    labels <- as.character(round(my.dat$pval , pval.digits))
    
    labels[labels=="0"] <- paste(c("<0.", rep(0, pval.digits), "5"), collapse = "")
    
    labels[my.dat$pval>alpha] <- ""
    
    lab_adjust <- fisylim0*0.1
    
    fisplot <- my.dat %>%
      dplyr::mutate(labels = !!labels) %>%
      ggplot2::ggplot(ggplot2::aes(y = fis, x = pop_no, label = labels)) +
      ggplot2::geom_abline(slope = 0, size = 0.5) +
      ggplot2::geom_line(color = "black",
                         linetype = "dashed",
                         linewidth = line.width) +
      ggplot2::geom_point(
        fill = PopCol,
        color = "black",
        shape = PopPch,
        size = point.size
      ) +
      ggplot2::ylim(fisylim * 1.1) +
      ggplot2::ylab(expression(italic(F)[IS])) +
      ggplot2::xlab("Population") +
      ggplot2::geom_text(ggplot2::aes(y = (fis + lab_adjust * sign(fis)), label = labels),
                         color = PopCol,
                         size = pval.cex) +
      ggplot2::theme_bw()
    
    grid::grid.newpage()
    grid::grid.draw(rbind(
      ggplot2::ggplotGrob(freqplot),
      ggplot2::ggplotGrob(fisplot),
      size = "last"
    ))
    
  }

  dev.off()
  
  print(Sys.time()-start.time)

  return(HWE_df %>% dplyr::select(population, locus, pval, fis, allele_freq = q))
  
}