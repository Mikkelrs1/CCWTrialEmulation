

#' Emulate trial with grace period using CCW approach.
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
#' @param plot True/False indicating whether Kaplan-Meier plots of the arms
#'                    should be returned aswell.
#' @return Results from the survival analysis evaluating the effect of the experimental treatment
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
#'  ),
#'  plot = T
#')
#' @export
#' @import survival dplyr survminer

Emulate_Trial <- function(data,
                          grace,
                          time,
                          formula_RHS,
                          var_list,
                          use.weights = T,
                          WC_limit = NULL,
                          plot = F){

  data <- check_and_format(data = data,
                           grace = grace,
                           time = time,
                           formula_RHS = formula_RHS,
                           var_list = var_list,
                           use.weights = use.weights,
                           WC_limit = WC_limit,
                           plot = plot)



  cloned_data <- make_cloned_data(data,
                                  grace
  )

  if(use.weights){
    data_to_use <- weight_cloned_data(cloned_data = cloned_data,
                                        formula_RHS = formula_RHS,
                                        grace = grace,
                                        WC_limit = WC_limit
    )
  } else{

    data_to_use <- cloned_data  |>  dplyr::mutate(Tstart = 0,
                                          Tstop = time_OS,
                                          weight = 1) |> dplyr::select(Tstart, Tstop, status_OS, arm, weight)
  }


  perform_analysis(data = data_to_use,
                   time = time,
                   plot = plot
  )
}
