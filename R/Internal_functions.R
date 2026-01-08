################################################################################################################
#### R script containing all internal functions used in the TrialEmulation R package ####
################################################################################################################

#######################################################################################
#### Internal function used for performing checks to ensure input is as expected   ####
########################################################################################
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
#' @param plot TRUE/FALSE specifying whether a plot should be returned.
#' @return A dataframe with all individuals from the original cloned (twice the number of rows) and where an addition "arm" (control/experimental) has been added.
#' @noRd
check_and_format <- function(data,
                             grace,
                             time,
                             formula_RHS,
                             var_list,
                             use.weights,
                             WC_limit,
                             plot) {

  #Check data
  if (!is.data.frame(data)) {
    stop("data should be a dataframe")
  }


  #Check grace
  if (!is.numeric(grace)) {
    stop("grace should be a numeric")
  }
  if (! grace > 0) {
    stop("grace should be strictly positive")
  }


    #Check use.weights
  if(!((is.logical(use.weights)) & length(use.weights) == 1)){
  stop("use.weights must be TRUE/FALSE")
  }

  #Check plot
  if(!((is.logical(plot)) & length(plot) == 1)){
    stop("plot must be TRUE/FALSE")
  }

    #Check formula_RHS
  if(use.weights){
    if(is.null(formula_RHS)){
      stop("formula_RHS must be specified when use.weights is TRUE")
    } else{
      if (!is.character(formula_RHS)) {
        stop("Formula_RHS should be a character string")
      }
      if (length(formula_RHS) != 1) {
        stop("Formula_RHS should have length 1")
      }

      split1 <- unlist(strsplit(formula_RHS, split = c("[+]")))
      split2 <- unlist(strsplit(split1, split = c("[:]")))
      split3 <- unlist(strsplit(split2, split = c("[*]")))

      covariates <- gsub(pattern = " |~", "", split3)

      if (!all(covariates %in% colnames(data))) {
        stop("Not all covariates indicated by formula_RHS could be found in the dataset")
      }
    }
  } else{
    if(!is.null(formula_RHS)){
      warning("formula_RHS is not used when use.weights is FALSE")
      covariates <- NULL
    }
  }

  #Check var_list

  if (!is.list(var_list)) {
    stop("var_list should be a list of the form list(time_OS = 'variable1', status_OS = 'variable2', time_treatment =  'variable3', status_treatment =  'variable4')")
  }

  if (length(var_list) != 4) {
    stop("var_list should have length 4")
  }

  if (!(all(names(var_list) %in% c("time_OS", "status_OS", "time_treatment", "status_treatment")) &
      all(c("time_OS", "status_OS", "time_treatment", "status_treatment") %in% names(var_list)))){
    stop("var_list should have names time_OS, status_OS, time_treatment, status_treatment")
  }

  if(!all(unlist(var_list) %in% colnames(data))){
    stop("Variable names specified in var_list should exists in data")
  }

  data$time_OS <- data[, var_list$time_OS]
  data$status_OS <- data[, var_list$status_OS]
  data$time_treatment <- data[, var_list$time_treatment]
  data$status_treatment <- data[, var_list$status_treatment]

  if(!is.numeric(data$time_OS)){
    stop("time_OS should be numeric")
  }

  if(!all(data$time_OS > 0)){
    stop("time_OS should be strictly positive")
  }


  if(!is.numeric(data$status_OS)){
    stop("status_OS should be numeric")
  }
  if(!all(data$status_OS %in% c(0,1))){
    stop("status_OS should only take values 0/1 (censored/death)")
  }

  if(!is.numeric(data$time_treatment)){
    stop("time_treatment should be numeric")
  }

  if(!all(data$time_treatment[!is.na(data$time_treatment)] > 0)){
    stop("time_treatment should be strictly positive (or NA)")
  }

  if(!all(is.na(data$time_treatment[data$status_treatment == 0]))){
    warning("time_treatment is set to NA when status_treatment is 0")
    data$time_treatment[data$status_treatment == 0] <- NA
  }


  if(!is.numeric(data$status_treatment)){
    stop("status_treatment should be numeric")
  }
  if(!all(data$status_treatment %in% c(0,1))){
    stop("status_treatment should only take values 0/1 (not treated/treated)")
  }



  #CHECKs OF THE WEBSTER-CLARK limit WC_limit

  event_times <- data$time_OS[data$status_OS == 1]

  if(!is.null(WC_limit)) {
    if(is.numeric(WC_limit) & length(WC_limit) == 1 & WC_limit > 0 & WC_limit < grace){
      stop("The Webster-Clark limit (WC_limit) should be a single positive numeric value smaller than the length of the grace window (grace)")
    }


    if(sum((grace > event_times & event_times > WC_limit)) == 0){
      stop("No event times occour within the grace window (grace) but after the Webster-Clark limit (WC_limit)")
    }
    if(sum((grace > event_times & event_times > Wc_limit)) < 10){
      warning("Less than 10 patients have within eventtimes in the grace window (grace) but after the Webster-Clark limit (WC_limit)")
    }
    if( !use.weights){
      warning("A Webster-Clark limit is specified (WC_limit) although weights are not used (use.weights)")
    }


  }


  #Last consistency checks


  if(grace < min(data$time_treatment, na.rm = T)){
    stop("All treatments are given after the grace period")
  }




  if(sum((grace > event_times))/length(event_times) > 0.5){
    warning("Over half of deaths ocours within the grace period. Consider a shorter grace period if possible")
  }


  if(time > max(data$time_OS)){
    stop("The provided timepoint to evaluate survival probabilities is beyond the available follow-up.")
  }

  #Prepare formatted data

  data <- data[, c("time_OS", "status_OS", "time_treatment", "status_treatment", covariates)]
  return(data)
}


#######################################################################################
#### Internal function used for cloning the data         ####
########################################################################################

#' Clone each individuals of the original data into a version in the experimental arm
#' and a version in the control arm.
#'
#' @param data Wide format dataframe containing the variables death (status indicator with values 0/1) and follow-up time
#' @param grace The length of the grace period.
#' @return A dataframe with all individuals from the original cloned (twice the number of rows) and where an addition "arm" (control/experimental) has been added.
#' @noRd
make_cloned_data <- function(data, grace) {
  #Create artificial arms
  control_arm <- data |>  dplyr::mutate(arm = "Control")
  experimental_arm <- data |>  dplyr::mutate(arm = "Experimental")


  #Add variables for OS outcome handling the new artifical censoring (AC)

  #Case 1: Patients who received  treatment within the grace window
  #Case 2: Otherwise (didnt receive treatment within the grace window)
  control_arm <- control_arm |>  dplyr::mutate(
    status_OS2 = dplyr::case_when(
      time_treatment <= grace ~ 0,
      T ~ status_OS
    ),
    time_OS2 = dplyr::case_when(
      time_treatment <= grace ~ time_treatment,
      T ~ time_OS
    )
  )


  #Case 1: Patients who received treatment within the grace window or died during the grace window
  #Case 2: Otherwise (Patient who did not receive treatment in but are still alive after grace window)
  experimental_arm <- experimental_arm |>  dplyr::mutate(
    status_OS2 = dplyr::case_when(
      (status_treatment == 1 & time_treatment <= grace) |
        time_OS < grace ~ status_OS,
      T ~ 0
    ),
    time_OS2 = dplyr::case_when(
      (status_treatment == 1 & time_treatment <= grace) |
        time_OS < grace ~ time_OS,
      T ~ grace
    )
  )


  #Add variables for time to artificial censoring (AC) due to cloning

  #Case 1: Receive treatment in grace window
  #Case 2: Dies/lost to follow-up in grace window before receiving treatment
  #Case 3: Not treated in the grace period with follow-up beyond the grace period.
  control_arm <- control_arm |>  dplyr::mutate(
    status_AC = dplyr::case_when(time_treatment <= grace ~ 1, time_OS <= grace ~ 0, T ~ 0),
    time_AC = dplyr::case_when(
      time_treatment <= grace ~ time_treatment,
      time_OS <= grace ~ time_OS,
      T ~ grace
    )
  )


  #Case 1: Receive treatment in grace window
  #Case 2: Dies/lost to follow-up in grace window before receiving treatment
  #Case 3: Not treated in the grace period with follow-up beyond the grace period
  experimental_arm <- experimental_arm |>  dplyr::mutate(
    status_AC = dplyr::case_when(time_treatment <= grace ~ 0, time_OS <= grace ~ 0, T ~ 1),
    time_AC = dplyr::case_when(
      time_treatment <= grace ~ time_treatment,
      time_OS <= grace ~ time_OS,
      T ~ grace
    )
  )

  cloned_data <- rbind(control_arm, experimental_arm)
  return(cloned_data |>  dplyr::select(-c(time_OS, status_OS, status_treatment)) |>   dplyr::rename(time_OS = time_OS2, status_OS = status_OS2))
}


#####################################################################################
#### Internal function used to split and weight the cloned data ####
#####################################################################################

#' Obtain patients that are dropped to the left in the tree based on the split
#'
#' @param cloned_data Data frame outputtet from the make_cloned_data function
#' which must contain the predictors indicated in the foruma_AC as variables.
#' @param formula_AC A formula object indicating the relationship between patient
#' characteristics and treatment decision (if a patient is treated and when).
#' @param grace The length of the grace period.
#' @param WC_limit The Webster-Clark limit to specify the which patients should be weigthed in
#' the experimental arm based on time_treatment. Only those treated after the WC_limit are
#' weighed. If NULL, all patients are weighted.
#' @return Long format data frame containing all clones splitting in accordance
#' with their changes IPCWs. Newly added variables include Tstart, Tstop and weight.
#' @noRd
weight_cloned_data <- function(cloned_data, formula_RHS , grace, WC_limit){

  cloned_data$Tstart <- 0
  cloned_data$Tstop <- cloned_data$time_OS

  ######################
  # External arm       #
  ######################
  experimental_arm <- cloned_data |>
    dplyr::filter(arm == "Experimental") |>
    survival::survSplit(cut = grace, end = "Tstop",
              start = "Tstart", event = "status_OS",id = "ID")


  #Add weights to the remnant population after artifical censoring at the end of the grace period
  if(is.null(WC_limit)){
    formula_AC_logistic <- as.formula(paste0("status_AC", formula_RHS))
    fit_logistic <- glm(formula_AC_logistic,
                        data = experimental_arm |>  dplyr::filter(time_OS >= grace),
                        family = "binomial")


    wh <- experimental_arm$ID |>  duplicated() #Utilizing the ordering of the data from survSplit, otherwise order to Tstart first

    weights = 1/(1 - predict(fit_logistic, newdata = experimental_arm[wh,], type = "response"))

    experimental_arm$weight = ifelse(wh, weights, 1)

    #EXACTLY THE SAME WEIGHTS AS IN THE ORIGINIAL IMPLEMENTATION
    # fit_cox <- coxph(Surv(Tstart, Tstop, status_AC) ~ age+sex+emergency+
    #                       stage+deprivation+charlson+perf, data = experimental_arm |>  dplyr::filter(time_OS >= grace))
    #
    #
    # experimental_arm$weight <- 1/predict(fit_cox, newdata = experimental_arm |>  dplyr::mutate(Tstart = 0, Tstop = Tstop - 0.001), type = "survival")
  } else{
    formula_AC_logistic <- as.formula(paste0("status_AC", formula_RHS))
    fit_logistic <- glm(
      formula_AC_logistic,
      data = experimental_arm |>  dplyr::filter(
        time_OS >= grace,
        is.na(time_treatment) | time_treatment >= WC_limit
      ),
      family = "binomial"
    )

    wh <- experimental_arm$ID |>  duplicated()

    weights = 1/(1 - predict(fit_logistic, newdata = experimental_arm[wh,], type = "response"))
    weights = ifelse(experimental_arm$time_treatment >= WC_limit, weights, 1)
    experimental_arm$weight = ifelse(wh, weights, 1)
  }


  ######################
  # Control arm        #
  ######################
  #browser()
  control_arm <- cloned_data |>
    dplyr::filter(arm == "Control")

  times_of_interest <- sort(unique(control_arm |>  dplyr::filter(status_AC == 1) |>  dplyr::pull(time_AC)))

  control_arm <- survival::survSplit(control_arm, cut = times_of_interest, end="Tstop",
                           start="Tstart", event="status_OS",id="ID")


  tmp <- control_arm |>  dplyr::group_by(ID) |>  dplyr::summarise(index = which.max(Tstart),
                                                    status_AC = dplyr::first(status_AC))

  control_arm$status_AC <- 0
  control_arm$status_AC[tmp$index |>  cumsum()] <- tmp$status_AC

  formula_AC_Cox <- as.formula(paste0("Surv(Tstart, Tstop, status_AC)", formula_RHS))
    fit_cox <- survival::coxph(formula_AC_Cox,
                   ties = "efron",
                   data = control_arm)

  control_arm$weight <- 1/predict(fit_cox, newdata = control_arm |>  dplyr::mutate(Tstart = 0, Tstop = Tstop - 0.001), type = "survival")


  data_weighted <- rbind(experimental_arm, control_arm) |>  dplyr::select(Tstart, Tstop, status_OS, arm, weight)

  return(data_weighted)
}

#####################################################################################
#### Internal functions used to perform analyses of interest ####
#####################################################################################

#' A function to perform the survival analysis on the data.
#'
#' @param weighted_data A dataframe as outputted by weight_cloned_data
#' @param time Fixed timepoint to assess survival and survival difference.
#' @return A list of results from the performed analysis.
#' @noRd
perform_analysis <- function(data, time, plot){
  KM <- survival::survfit(Surv(Tstart, Tstop, status_OS) ~ arm,
                data = data, weights = weight)
  OS <- summary(KM, times = time, extend = T)$surv
  OS[3] <- OS[2] - OS[1]



  RMS <- summary(KM, rmean = time)$table[,"rmean"]
  RMS[3] <- RMS[2] - RMS[1]


  HR <- survival::coxph(Surv(Tstart, Tstop, status_OS) ~ arm, data = data, weights = weight) |>  coef() |>  exp() |>  unname()


  res <- c(OS, RMS, HR)
  names(res) <- c("OS_Control", "OS_Experimental", "OS_Difference",
                  "RMS_Control", "RMS_Experimental", "RMS_Difference",
                  "Hazard ratio")


  if(plot){
    splot <- survminer::ggsurvplot(fit = KM,
                                   data = data,
                                   pval = F,
                                   censor = F,
                                   fontsize = 5,
                                   size = 0.2,
                                   conf.int = FALSE,
                                   legend.title = "",
                                   legend.labs = c("Control arm", "Experimental arm"),
                                   main = "",
                                   xlab = "Time from diagnosis",
                                   ylab = "Survival probability")

    out <- list(res, splot$plot)
  } else{
    out <- res
  }

  return(out)
}




#####################################################################################
#### Internal function used to obtain confidence intervals of interest based on bootstrap samples
#####################################################################################

#' @return a%-confidence intervals for each metric
#' @noRd
get_CI <- function(boot_res, CI_type, a){
  lower_quantile <- (1-a/100)/2
  upper_quantile <- 1 - lower_quantile
  if(CI_type == "gaussian"){
    apply(boot_res, MARGIN = 1, FUN = function(sample) {
      mean <- mean(sample)
      se <- sd(sample)

      lower <- mean + qnorm(lower_quantile) * se
      upper <- mean + qnorm(upper_quantile) * se
      return(c("lower" = lower, "upper" = upper))
    })
  } else{ #CI_type = "quantile"
    apply(boot_res, MARGIN = 1, FUN = quantile, probs = c(lower_quantile, upper_quantile))
  }
}
