#' @title Distribute liver drug
#'
#' @param timestep current timestep
#' @param variables the available variables
#' @param events the available events
#' @param parameters model parameters
#' @param correlations model correlations
#' @param renderer object for model outputs
#' @noRd
create_liverdrug_listener <- function(variables, events, parameters, correlations, renderer) {
  function(timestep) {
    time_index = which(parameters$liverdrug_timesteps == timestep)
    ages_d<-get_age(variables$birth$get_values(), timestep)
    target <- which(ages_d > parameters$liverdrug_min_age & ages_d < parameters$liverdrug_max_age)
    to_vaccinate <- target[sample_intervention(
      target,
      'liverdrug',
      parameters$liverdrug_coverages[[time_index]],
      correlations
    )]
    renderer$render('n_vaccinated_liverdrug', length(to_vaccinate), timestep)
    if (length(to_vaccinate) > 0) {
      variables$liverdrug_vaccinated$queue_update(
        timestep,
        to_vaccinate
      )
    }
    if (time_index < length(parameters$liverdrug_timesteps)) {
      events$liverdrug_vaccination$schedule(
        parameters$liverdrug_timesteps[[time_index + 1]] - timestep
      )
    }
  }
}



calculate_liverdrug_efficacy <- function(t,parameters) {
  parameters$liverdrug_prophylaxis[t]
}

