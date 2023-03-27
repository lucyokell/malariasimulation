#' @title Preset parameters for the DHA-PQP drug
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @details Use a list of preset parameters for the DHA-PQP drug (dihydroartemisinin-piperaquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.95, drug_rel_c: 0.09434, drug_prophylaxis_shape: 4.4, drug_prophylaxis_scale: 28.1, user_prophylaxis (list)=NA
#' @export
DHA_PQP_params <- list(.95, 0.09434, 4.4, 28.1,NA)

#' @title Preset parameters for the AL drug
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @details Use a list of preset parameters for the AL drug (artemether-lumefantrine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.95, drug_rel_c: 0.05094, drug_prophylaxis_shape: 11.3, drug_prophylaxis_scale: 10.6, user_prophylaxis=NA
#' @export
AL_params <- list(.95, 0.05094, 11.3, 10.6,NA)

#' @title Preset parameters for the SP-AQ drug
#' @description From SI of Thompson et al (2022). doi: 10.1016/S2214-109X(22)00416-8
#' @details Use a list of preset parameters for the SP-AQ drug (sulphadoxine-pyrimethamine and amodiaquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.9, drug_rel_c: 0.32, drug_prophylaxis_shape: 4.3, drug_prophylaxis_scale: 38.1, user_prophylaxis=NA
#' @export
SP_AQ_params <- list(0.9, 0.32, 3.4, 39.34,NA)

#' @title Parameterise drugs to use in the model
#'
#' @param parameters the model parameters
#' @param drugs a list of drug parameters, can be set using presets
#' @export
set_drugs <- function(parameters, drugs) {
  keys <- c(
    'drug_efficacy',
    'drug_rel_c',
    'drug_prophylaxis_shape',
    'drug_prophylaxis_scale',
    'user_prophylaxis'
  )
  for (drug in seq_along(drugs)) {  ## add new loop which resets drugs rather than adding to existing ones.
    for (i in seq_along(drugs[[drug]])) {
      parameters[[keys[[i]]]]<-numeric(0)
    }
  }
  for (drug in seq_along(drugs)) {
    for (i in seq_along(drugs[[drug]])) {
      if(i<5) { parameters[[keys[[i]]]] <- c(parameters[[keys[[i]]]], drugs[[drug]][[i]])
      } else {
        parameters[[keys[[i]]]] <- append(parameters[[keys[[i]]]], list(drugs[[drug]][[i]]))
      }
    }
  }
  parameters
}

#' @title Prophylaxis: obtain probability of protection against new infection due to prior treatment
#'
#' @param time_since_drug vector of times since individuals were treated
#' @param drug_prophylaxis_shape vector of Weibull shape parameters describing prophylactic profile of drug
#' @param drug_prophylaxis_scale vector of Weibull scale parameters describing prophylactic profile of drug
#' @param user_prophylaxis list of user-defined dataframes describing prophylactic protection. If not missing, this is used instead of drug_prophylaxis_shape and drug_prophylaxis_scale
#' @export
#' 
get_prophylaxis <- function(time_since_drug,drug_prophylaxis_shape,drug_prophylaxis_scale,user_prophylaxis) {
  prophylaxis<-rep(0,length(time_since_drug))
  is.weibull<-is.na(user_prophylaxis)
  if(any(is.weibull)) {
    prophylaxis[is.weibull] <- weibull_survival(
      time_since_drug[is.weibull],
      drug_prophylaxis_shape[is.weibull],
      drug_prophylaxis_scale[is.weibull]
    )
  }
  if(any(!is.weibull)) {
    indexes<-which(!is.weibull)
    for(i in indexes) {
      prophylaxis[i]<-ifelse(time_since_drug[i]>max(user_prophylaxis[[i]]$day),0,
                             user_prophylaxis[[i]]$prophylaxis[which(user_prophylaxis[[i]]$day==time_since_drug[i])])
    }
  }
  return(prophylaxis)
}


#' @title Parameterise clinical treatment
#'
#' @param parameters the model parameters
#' @param drug the index of the drug to set as a treatment
#' @param timesteps vector of timesteps for each coverage change
#' @param coverages vector of coverages for this drug
#' @export
set_clinical_treatment <- function(parameters, drug, timesteps, coverages) {
  n_drugs <- length(parameters$drug_efficacy)
  if (drug < 1 | drug > n_drugs) {
    stop('Drug index is invalid, please set drugs using set_drugs')
  }
  drug_index <- which(parameters$clinical_treatment_drugs == drug)
  if (length(drug_index) == 0) {
    drug_index <- length(parameters$clinical_treatment_drugs) + 1
  }
  parameters$clinical_treatment_drugs[[drug_index]] <- drug
  parameters$clinical_treatment_timesteps[[drug_index]] <- timesteps
  parameters$clinical_treatment_coverages[[drug_index]] <- coverages
  last_timestep <- max(unlist(parameters$clinical_treatment_timesteps))

  for (t in seq(last_timestep)) {
    if (sum(get_treatment_coverages(parameters, t)) > 1) {
      stop('The sum of drug coverages cannot be greater than 1 at any timestep')
    }
  }
  parameters
}

get_treatment_coverages <- function(parameters, timestep) {
  vnapply(
    seq_along(parameters$clinical_treatment_drugs),
    function(drug_index) {
      previous <- which(
        parameters$clinical_treatment_timesteps[[drug_index]] <= timestep
      )
      if (length(previous) == 0) {
        return(0)
      }
      last_set <- max(previous)
      parameters$clinical_treatment_coverages[[drug_index]][[last_set]]
    }
  )
}

