#' Bootstrap confidence intervals of the estimates obtained from Emulate_Trial.
#'
#' @param data A wide format data frame
#' @param grace The length of the grace period.
#' @param time Fixed time point to assess survival and survival difference.
#' @param formula_RHS Right-hand side of formula for the relationship between
#'                    if and when a patient is treated and his covariates.
#'                    Is used to model the artifical censoring to obtain weights.
#' @param var_list A named list naming the variables in data that should be used
#'                 time_OS, status_OS, time_treatment, status_treatment.
#' @param use.weights True false indicating whether the trial emulation should use
#' weights to adjust for the selection bias/informative censoring introduced by the
#' artificial censoring from the cloning strategy.
#' @param WC_limit The Webster-Clark limit to specify the which patients should be weigthed in
#' the experimental arm based on time_treatment. Only those treated after the WC_limit are
#' weighed. If NULL, all patients are weighted.
#' @param ncores Number of cores used for parrelized computations.
#' @param mboot Number of bootstrap repitions used.
#' @param CI_type The choice of confidence interval used, either "none", "normal" or "quantile"
#' @param a corresponding the an a\%-confidence interval.
#' @return Confidence intervals of the survival estimates
#' @include Internal_functions.R
#' @examples
#' library(survival)
#' library(dplyr)
#' library(CCWTrialEmulation)
#'
#'
#'Emulate_Trial(
#'  data = jasa |> mutate(
#'    futime = ifelse(futime == 0, 1, futime),
#'    #does not accept events at time 0
#'    wait.time = ifelse(wait.time == 0, 1, wait.time)
#'  ),
#'  grace = 60,
#'  time = 1 * 365,
#'  formula_RHS = "~ age + surgery",
#'  var_list = list(
#'    "time_OS" = "futime",
#'    "status_OS" = "fustat",
#'    "time_treatment" = "wait.time",
#'    "status_treatment" = "transplant"
#'  )
#')
#'
#'Emulate_Trial_Bootstrap(
#'  jasa |> mutate(
#'    futime = ifelse(futime == 0, 1, futime),
#'    #does not accept events at time 0
#'    wait.time = ifelse(wait.time == 0, 1, wait.time)
#'  ),
#'  grace = 60,
#'  time = 1*365,
#'  formula_RHS = "~ age + surgery",
#'  var_list = list(
#'    "time_OS" = "futime",
#'    "status_OS" = "fustat",
#'    "time_treatment" = "wait.time",
#'    "status_treatment" = "transplant"
#'  ),
#'  use.weights = T,
#'  ncores = 5,
#'  mboot = 100,
#'  CI_type = "gaussian",
#'  a = 95
#')
#'
#'
#'
#'
#' @export
#' @import parallel utils survival dplyr

Emulate_Trial_Bootstrap <- function(data,
                                    grace,
                                    time,
                                    formula_RHS = NULL,
                                    var_list,
                                    use.weights = T,
                                    WC_limit = NULL,
                                    ncores = 1,
                                    mboot,
                                    CI_type = "gaussian",
                                    a = 95) {

  data <- check_and_format(data = data,
                           grace = grace,
                           time = time,
                           formula_RHS = formula_RHS,
                           var_list = var_list,
                           use.weights = use.weights,
                           WC_limit = WC_limit,
                           plot = F)

  #OBS
  #Make checks for ncores, mboot, CI_type, a.
  #eg. warning if ncores > mboot

  #Check a
  if(!(is.numeric(a) & length(a) == 1)){
    stop("a should be a numeric of length 1")
  }

  if(!(a>0 & a < 100)){
    stop("a should be strictly between 0 and 100")
  }


  #Check mboot
  if(!(is.numeric(mboot) & length(mboot) == 1)){
    stop("mboot should be an integer")
  }


  #Check ncores
  if(!(is.numeric(ncores) & length(ncores) == 1)){
    stop("ncores should be an integer")
  }

  max_cores <- parallel::detectCores()

  if(ncores > mboot){
    warning("ncores is larger than the number of bootstrap repitions (mboot). ncores is reduced to mboot")
    ncores <- mboot
  }

  if(ncores > max_cores){
    warning("ncores is larger than the number of cores available. All cores are used")
    ncores <- max_cores
  }

  cloned_data <- make_cloned_data(data, grace)
  n <- nrow(data)


  if (ncores == 1) { #No parrelization

run <- function(){
  indicies <- sample(1:n, size = n, replace = T)

  control_arm <- cloned_data[1:n, ]
  experimental_arm <- cloned_data[(n + 1):(2*n), ]

  cloned_data <- rbind(control_arm[indicies, ], experimental_arm[indicies, ])

  if (use.weights) {
    data_to_use <- weight_cloned_data(cloned_data, formula_RHS = formula_RHS, grace = grace, WC_limit = WC_limit)
  } else{
    data_to_use <- cloned_data |> dplyr::mutate(Tstart = 0,
                                          Tstop = time_OS,
                                          weight = 1)  |>  dplyr::select(Tstart, Tstop, status_OS, arm, weight)
  }
  perform_analysis(data = data_to_use,
                   time = time,
                   plot = F)
}
    boot_res <- replicate(n = mboot, expr = run())

    if(CI_type == "none"){
      return(boot_res)
    } else{
      CI <- get_CI(boot_res, CI_type, a)
      return(CI)
    }
  } else{
   # stop("ncores > 1 not yet implemented")

        cluster <- parallel::makeCluster(ncores)
        parallel::clusterExport(cl = cluster,
                                unclass(utils::lsf.str(
                                  envir = asNamespace("CCWTrialEmulation"),
                                  all = T
                                )),
                                envir = as.environment(asNamespace("CCWTrialEmulation")))
        parallel::clusterExport(envir = environment(),
                                cl = cluster,
                                ls(environment()))

        boot_res <- parallel::parSapply(
          cl = cluster,
          1:1000,
          FUN = function(i) { #i isnt used, this functions are a paralelized replicate.

              indicies <- sample(1:n, size = n, replace = T)

              control_arm <- cloned_data[1:n, ]
              experimental_arm <- cloned_data[(n + 1):(2*n), ]

              cloned_data <- rbind(control_arm[indicies, ], experimental_arm[indicies, ])

              if (use.weights) {
                data_to_use <- weight_cloned_data(cloned_data, formula_RHS = formula_RHS, grace = grace, WC_limit = WC_limit)
              } else{
                data_to_use <- cloned_data %>% mutate(Tstart = 0,
                                                      Tstop = time_OS,
                                                      weight = 1) %>% select(Tstart, Tstop, status_OS, arm, weight)
              }
              perform_analysis(data = data_to_use,
                               time = time,
                               plot = F)
          }
        )
        parallel::stopCluster(cl = cluster)
        rm(cluster)
        if(CI_type == "none"){
          return(boot_res)
        } else{
          CI <- get_CI(boot_res, CI_type, a)
          return(CI)
        }
  }
}
