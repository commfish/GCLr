#' Create confusion matrices 
#'
#' This function is a wrapper for rubias::self_assign(). It takes the genetic likelihood profile produced by rubias::self_assign() and summarizes it into tibbles at the group by group, group by pop, and pop by pop level. This genetic likelihood profile is intended to show how well populations/reporting groups are differentiated by the markers at hand. These probabilities can give an idea of potential misallocation and where. This function can aid in marker selection.
#'
#' @param reference A rubias baseline object created by create_rubias_base().
#' 
#' @param gen_start_col The first column containing genetic data in the baseline object. This will usually be the default unless the baseline object is modified with additional variables.
#'
#' @param output A character vector of summary tibbles you want to produce. Must be one or more of the following: "pop_pop", "pop_group", "group_group". The default is "group_group".
#'
#' @return Named list of tibble(s) in long format:
#'   \itemize{
#'     \item \code{group_group}: A tibble containing the following variables:
#'       \itemize{
#'         \item \code{repunit}: The known repunit (aka reporting group).
#'         \item \code{inferred_repunit}: The inferred repunit.
#'         \item \code{mean_group_group_scaled_like}: The average probability of each individual from a group originating from a pop in that group.
#'       }
#'     \item \code{pop_group}: A tibble containing the following variables:
#'       \itemize{
#'         \item \code{collection}: The known population.
#'         \item \code{inferred_repunit}: The inferred repunit.
#'         \item \code{mean_pop_group_scaled_like}: The average probability of each individual from a pop originating from a pop in that group.
#'       }
#'     \item \code{pop_pop}: A tibble containing the following variables:
#'       \itemize{
#'         \item \code{pop}: The known population.
#'         \item \code{inferred_pop}: The inferred repunit.
#'         \item \code{mean_pop_pop_scaled_like}: The average probability of each individual from a pop originating from a pop.
#'       }
#'   }
#'
#' @examples
#' baseline <- read_csv("V:/Analysis/5_Coastwide/Chum/NPen2WA_Chum_baseline/rubias/baseline/NPen2Wa_Chum_227pops_91loci_base.csv")
#' ConfusionMatrices_out <- confusion_matrix(reference = baseline , gen_start_col = 5, output = c("group_group", "pop_group", "pop_pop"))
#'
#' @import rubias
#' @import dplyr
#' @import purrr
#'
#' @export

confusion_matrix <- function(reference, gen_start_col = 5, output = c("group_group", "pop_group", "pop_pop")[1]){
  
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