#' This function takes the genetic likelihood profile produced by rubias::self_assign() and summarizes it into tibbles at the group-by-group, group-by-pop, and pop-by-pop level. 
#' 
#' The function provides insights into how well populations/reporting groups are differentiated by the genetic markers and can help in marker selection.
#'
#' @param reference A rubias baseline object created by create_rubias_baseline().
#' @param gen_start_col The first column containing genetic data in the baseline object. The default is 5, unless the baseline object is modified with additional variables.
#' @param output A character vector specifying the summary tibbles to produce. Must be one or more of the following: "pop_pop", "pop_group", "group_group". The default is "group_group".
#'
#' @return A named list of tibbles in long format:
#'   \describe{
#'     \item{group_group}{A tibble containing the following variables: repunit (known repunit/ reporting group), inferred_repunit (inferred repunit), and mean_group_group_scaled_like (average probability of each individual from a group originating from a pop in that group).}
#'     \item{pop_group}{A tibble containing the following variables: collection (known population), inferred_repunit (inferred repunit), and mean_pop_group_scaled_like (average probability of each individual from a pop originating from a pop in that group).}
#'     \item{pop_pop}{A tibble containing the following variables: pop (known population), inferred_pop (inferred repunit), and mean_pop_pop_scaled_like (average probability of each individual from a pop originating from a pop).}
#'   }
#'
#' @examples
#' baseline <- read_csv("V:/Analysis/5_Coastwide/Chum/NPen2WA_Chum_baseline/rubias/baseline/NPen2Wa_Chum_227pops_91loci_base.csv")
#' ConfusionMatrices_out <- ConfusionMatrices.GCL(reference = baseline, gen_start_col = 5, output = c("group_group", "pop_group", "pop_pop"))

ConfusionMatrices.GCL <- function(reference, gen_start_col = 5, output = c("group_group", "pop_group", "pop_pop")[1]){
    
  start <- Sys.time()
  
  rubias_sa <- rubias::self_assign(reference = baseline, gen_start_col = gen_start_col) %>% 
    select(-missing_loci, -n_miss_loci, -n_non_miss_loci, -z_score, -log_likelihood)
  
  pop_pop <- function(rubias_sa){
    
    rubias_sa %>% 
      dplyr::group_by(indiv, collection, inferred_collection) %>% 
      dplyr::summarise(pop_scale_like = sum(scaled_likelihood), .groups = "drop_last") %>% 
      dplyr::group_by(collection, inferred_collection) %>% 
      dplyr::summarise(mean_pop_pop_scaled_like = mean(pop_scale_like), .groups = "drop_last") %>% 
      dplyr::rename(pop = collection, inferred_pop = inferred_collection)
    
  }
  
  pop_group <-  function(rubias_sa){
    
    rubias_sa %>% 
      dplyr::group_by(indiv, collection, inferred_repunit) %>% 
      dplyr::summarise(pr_scale_like = sum(scaled_likelihood), .groups = "drop_last") %>% 
      dplyr::group_by(collection, inferred_repunit) %>% 
      dplyr::summarise(mean_pop_group_scaled_like = mean(pr_scale_like), .groups = "drop_last") %>% 
      dplyr::rename(pop = collection, inferred_group = inferred_repunit)
    
  }
  
  group_group <- function(rubias_sa){ 
    
  rubias_sa %>% 
      dplyr::group_by(indiv, collection, repunit, inferred_repunit) %>%
      dplyr::summarise(repu_scaled_like = sum(scaled_likelihood), .groups = "drop_last") %>% 
      dplyr::group_by(collection, repunit, inferred_repunit) %>% 
      dplyr::summarise(mean_repu_scaled_like_col = mean(repu_scaled_like), .groups = "drop_last") %>% 
      dplyr::group_by(repunit, inferred_repunit) %>% 
      dplyr::summarise(mean_group_group_scaled_like = mean(mean_repu_scaled_like_col), .groups = "drop_last") %>% 
      dplyr::rename(group = repunit, inferred_group = inferred_repunit)
    
  }

  out_sum <- lapply(output, function(out){
    
    func <- get(out)
    
    func(rubias_sa)
    
  }) %>% purrr::set_names(output)
  
  print(Sys.time()-start)
  
  return(out_sum)

}