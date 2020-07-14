/*
 * mosquito_emergence.cpp
 *
 *  Created on: 15 May 2020
 *      Author: gc1610
 */

#include <Rcpp.h>
#include <individual.h>
#include "mosquito_biology.h"

//' @title Mosquito births
//' @description
//' This is the process for mosquito birth, it defines how many new early stage
//' larvae are created on each timestep.
//' @param mosquito the mosquito individual
//' @param susceptable the susceptable mosquito state
//' @param infected the infected mosquito state
//' @param unborn the unborn mosquito state
//' @param early_larval_stage the early stage larval state
//[[Rcpp::export]]
Rcpp::XPtr<process_t> create_egg_laying_process_cpp(
    std::string mosquito,
    std::string susceptable,
    std::string incubating,
    std::string infected,
    std::string unborn,
    std::string early_larval_stage
    ) {
    auto process = [=](ProcessAPI& api) {
        auto n_M = api.get_state(mosquito, susceptable).size() +
            api.get_state(mosquito, incubating).size() +
            api.get_state(mosquito, infected).size();
        const auto& u = api.get_state(mosquito, unborn);
        const auto& parameters = api.get_parameters();
        if (n_M > 0) {
            auto n_eggs = parameters.at("beta")[0] * n_M;
            if (n_eggs > u.size()) {
                Rcpp::stop("Run out of mosquitos");
            }
            if (n_eggs >= 1) {
                auto target = std::vector<size_t>(n_eggs);
                auto it = u.cbegin();
                for (auto i = 0; i < target.size(); ++i) {
                    target[i] = *it;
                    ++it;
                }
                api.queue_state_update(mosquito, early_larval_stage, target);
            }
        }
    };
    return Rcpp::XPtr<process_t>(new process_t(process));
}


void deaths(const individual_index_t& source, std::vector<size_t>& target, double rate) {
    auto uniform = Rcpp::runif(source.size());
    auto i = 0u;
    for (auto s : source) {
        if (uniform[i] < rate) {
            target.push_back(s);
        }
        ++i;
    }
}

//[[Rcpp::export]]
Rcpp::XPtr<process_t> create_larval_death_process_cpp(
    std::string mosquito,
    std::string early_larval_stage,
    std::string late_larval_stage,
    std::string unborn,
    double K0,
    double R_bar
    ) {
    auto process = [=](ProcessAPI& api) {
        const auto timestep = api.get_timestep();
        const auto& parameters = api.get_parameters();
        const auto& early_larval = api.get_state(mosquito, early_larval_stage);
        const auto& late_larval = api.get_state(mosquito, late_larval_stage);
        auto n = early_larval.size() + late_larval.size();
        auto k = carrying_capacity(timestep, parameters, K0);
        auto early_regulation = 1 + n / k;
        auto late_regulation = 1 + parameters.at("gamma")[0] * n / k;
        auto larval_deaths = std::vector<size_t>();
        deaths(early_larval, larval_deaths, parameters.at("me")[0] * early_regulation);
        deaths(late_larval, larval_deaths, parameters.at("ml")[0] * late_regulation);
        api.queue_state_update(mosquito, unborn, larval_deaths);
    };
    return Rcpp::XPtr<process_t>(new process_t(process));
}

