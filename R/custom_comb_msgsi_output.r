#' Custom Combine \pkg{Ms.GSI} Output
#'
#' This function computes summary statistics from \pkg{Ms.GSI} output, similar to `custom_comb_rubias_output()`.
#' It is primarily used to "roll up" reporting groups from fine-scale to broad-scale.
#' Output is a single tibble with `mixture` as a column. It can take either the `mdl_out` list 
#' object from [Ms.GSI::msgsi_mdl()], OR it can read in the .csv output files created by [Ms.GSI::msgsi_mdl()].
#'
#' @param mdl_out Output list object from [Ms.GSI::msgsi_mdl()].
#' @param path Path to \pkg{Ms.GSI} directory containing `mix` folder with output .csv files specified by `file_path` argument in [Ms.GSI::msgsi_mdl()].
#' @param mix Character vector of length 1 with mixture silly, used to read in output .csv files if `mdl_out = NULL`.
#' @param new_pop_info Tibble with 2 columns: `repunit` with existing fine-scale groups and `new_repunit` for rolled up broad-scale groups. 
#'                     Same object as [Ms.GSI::stratified_estimator_msgsi()] input.
#' @param group_names_new Optional character vector of new group_names for ordering as factor.
#' @param nreps Total number of iterations (includes burn-ins) used in [Ms.GSI::msgsi_mdl()]. Only used if not providing `path` to folder containing output .csv files.
#' @param nburn Number of warm-up iterations used in [Ms.GSI::msgsi_mdl()]. Only used if not providing `path` to folder containing output .csv files.
#' @param thin Frequency to thin iterations in the output used in [Ms.GSI::msgsi_mdl()]. Only used if not providing `path` to folder containing output .csv files.
#' @param nchains Number of independent MCMC chains used in [Ms.GSI::msgsi_mdl()]. Only used if not providing `path` to folder containing output .csv files.
#' @param keep_burn Boolean to save the burn-in iterations or not used in [Ms.GSI::msgsi_mdl()], default is `FALSE`. Only used if not providing `path` to folder containing output .csv files.
#' 
#' @return A tibble with 9 fields.
#'   - \code{group}: Character vector of reporting groups specified in `new_pop_info$new_repunit` (factor for ordering, plotting purposes if `group_names_new` provided).
#'   - \code{mean}: Mean stock composition.
#'   - \code{median}: Median stock composition.
#'   - \code{sd}: Standard deviation.
#'   - \code{ci.05}: Lower bound of 90% credibility interval.
#'   - \code{ci.95}: Upper bound of 90% credibility interval.
#'   - \code{P=0}: The proportion of the stock comp distribution that was below "threshold", either 5e-7 or ~1/2 fish if `harvest` is provided (i.e., posterior probability that stock comp = 0).
#'   - \code{GR}: Gelman-Rubin diagnostic used to assess convergence among MCMC chains, GCL standard is < 1.2.
#'   - \code{n_eff}: Effective sample size is an estimate of independent sample size of the posterior sample and is used to assess MCMC convergence. No official "threshold", but generally larger is better.
#'
#' @examples
#' \dontrun{
#' AY_2024_GAPS_comparison_34RG_msgsi_estimates <- sapply(AY_2024_GAPS_comparison_mixnames, function(mix) {
#' custom_comb_msgsi_output(
#'  path = "msgsi/output/",
#'  mix = mix,
#'  new_pop_info = GTseq_GAPS_RG_conversion %>% 
#'    dplyr::select(repunit = `58_SEAK_CW_AK`, new_repunit = hybrid_GAPS_33_pub) %>% 
#'    dplyr::distinct(),
#'  group_names_new = groups_33_gaps_33
#' ) %>%
#'  dplyr::mutate(mixture = mix) %>%
#'  dplyr::relocate(mixture)
#' }, simplify = FALSE) %>%
#'  dplyr::bind_rows() %>% 
#'  dplyr::mutate(mixture = factor(x = mixture, levels = AY_2024_GAPS_comparison_mixnames))
#' }
#'
#' @export
custom_comb_msgsi_output <- function(mdl_out = NULL,
                                     path = NULL,
                                     mix = NULL,
                                     new_pop_info,
                                     group_names_new = NULL,
                                     nreps = NULL,
                                     nburn = NULL,
                                     thin = NULL,
                                     nchains = NULL,
                                     keep_burn = NULL
) {
  
  
  if (is.null(mdl_out) & is.null(path)) {
    stop("`Either provide a Ms.GSI output or a valid `path` to folders containing output .csv files for each mixture.")
  }
  
  # get msgsi_specs and trace_comb
  if (is.null(mdl_out)) {
    grp_info <- readr::read_csv(file = file.path(path, mix, "comb_groups.csv"),
                                show_col_types = FALSE)
    
    msgsi_specs <- readr::read_csv(file = file.path(path, mix, "msgsi_specs.csv"),
                                   show_col_types = FALSE) %>%
      tibble::column_to_rownames("name")
    nreps <- msgsi_specs["nreps",]
    nburn <- msgsi_specs["nburn",]
    thin <- msgsi_specs["thin",]
    nchains <- msgsi_specs["nchains",]
    keep_burn <- as.logical(msgsi_specs["keep_burn",])
    
    trace_comb <- readr::read_csv(file = file.path(path, mix, "trace_comb.csv"),
                                  show_col_types = FALSE)
    
  } else {
    grp_info <- mdl_out$comb_groups
    trace_comb <- mdl_out$trace_comb
    
    if(any(sapply(list(nreps, nburn, thin, nchains, keep_burn), is.null))) {
      stop("`If not providing `path` to folders containing output .csv files for each mixture, need to specify `nreps`, `nburn`, `thin`, `nchains`, and `keep_burn`.")
    }
    
  }
  
  if (is.null(group_names_new)) {
    group_names_new <- unique(new_pop_info$new_repunit)
  }  # used to order groups as a factor
  
  trace_comb_new_grp <- trace_comb %>% 
    tidyr::pivot_longer(cols = -c(itr, chain), names_to = "collection", values_to = "rho") %>% 
    dplyr::left_join(y = grp_info, by = "collection") %>%  # get combined collections to repunit
    dplyr::left_join(y = new_pop_info, by = "repunit") %>%  # get new repunit from new_pop_info, same as in Ms.GSI::stratified_estimator_msgsi
    dplyr::mutate(repunit = new_repunit) %>% 
    dplyr::select(-new_repunit, -grpvec) %>% 
    dplyr::group_by(chain, itr, repunit) %>%  # roll up from collection (populations) to new repunits
    dplyr::summarise(rho = sum(rho), .groups = "drop") %>% 
    tidyr::pivot_wider(names_from = repunit, values_from = rho)
  
  # need to create `p_combo` from trace_comb so I can make `mc_pop`
  p_combo <- lapply(seq(nchains), function(chn) {
    trace_comb_new_grp %>% 
      dplyr::filter(chain == chn) %>% 
      dplyr::select(-itr, -chain)
  })
  
  keep_list <- ((nburn*keep_burn + 1):(nreps - nburn * isFALSE(keep_burn)))[!((nburn*keep_burn + 1):(nreps - nburn * isFALSE(keep_burn))) %% thin] / thin  # these are from Ms.GSI::msgsi_mdl input parameters, also saved as msgsi_specs.csv
  
  mc_pop <- coda::as.mcmc.list(
    lapply(p_combo, function(rlist) coda::mcmc(rlist[keep_list,]))
  )
  
  summ_func <- function(combo_file, keeplist, mc_file, groupnames, n_ch) {
    
    lapply(combo_file, function(rlist) rlist[keeplist, ]) %>%
      dplyr::bind_rows() %>%
      tidyr::pivot_longer(cols = tidyr::everything()) %>%
      dplyr::summarise(
        mean = mean(value),
        median = stats::median(value),
        sd = stats::sd(value),
        ci.05 = stats::quantile(value, 0.05),
        ci.95 = stats::quantile(value, 0.95),
        p0 = mean(value < 5e-7),  # need to update to incorporate harvest, like Bobby said; Numeric constant specifying how low stock comp is before assuming 0, used for `P=0` calculation, default is from BAYES/rubias.
        .by = name
      ) %>%
      dplyr::mutate(
        GR = { if (n_ch > 1) {
          coda::gelman.diag(mc_file,
                            transform = FALSE,
                            autoburnin = FALSE,
                            multivariate = FALSE)$psrf[,"Point est."]
        } else NA },
        n_eff = coda::effectiveSize(mc_file) # in alphabetical order like tidyverse summarise()
      ) %>%
      dplyr::mutate(name_fac = factor(name, levels = groupnames)) %>%
      dplyr::arrange(name_fac) %>%
      dplyr::select(-name_fac) %>%
      dplyr::rename(group = name)
    
  }  # summary function
  
  summ_pop_new_grp <- summ_func(
    combo_file = p_combo,
    keeplist = keep_list,
    mc_file = mc_pop,
    groupnames = group_names_new,
    n_ch = nchains
  ) %>% 
    dplyr::mutate(group = factor(x = group, levels = group_names_new))
  
  return(summ_pop_new_grp)
  
}  # end