#ifndef _SOLVER_H
#define _SOLVER_H

#include <vector>
#include <memory>

#include "Objects.h"

void source_iteration( int& N_iter,
                       const double epsilon,
                       const std::vector<std::shared_ptr<Region>>& mesh,
                       const std::vector<std::shared_ptr<Region>>& region,
                       const std::vector<double>& mu, 
                       const std::vector<double>& w,
                       const std::shared_ptr<BC>& BC_left,
                       const std::shared_ptr<BC>& BC_right,
                       std::vector<double>& phi,
                       std::vector<std::vector<double>>& psi,
                       std::vector<double>& rho,
                       const std::string space_method,
                       const std::string accelerator_type );

void source_iteration_TD_implicit( 
        const double epsilon,
        const std::vector<std::shared_ptr<Region>>& mesh,
        const std::vector<std::shared_ptr<Material>>& material,
        const std::vector<std::shared_ptr<Region>>& region,
        const std::vector<double>& mu, 
        const std::vector<double>& w,
        const std::shared_ptr<BC>& BC_left,
        const std::shared_ptr<BC>& BC_right,
        const double speed, const double dt,
        const std::vector<std::vector<double>>& psi_initial,
        std::vector<std::vector<double>>& phi,
        const std::string accelerator_type,
        const std::vector<double>& time );

void source_iteration_TD_MB( 
        const double epsilon,
        const std::vector<std::shared_ptr<Region>>& mesh,
        const std::vector<std::shared_ptr<Material>>& material,
        const std::vector<std::shared_ptr<Region>>& region,
        const std::vector<double>& mu, 
        const std::vector<double>& w,
        const std::shared_ptr<BC>& BC_left,
        const std::shared_ptr<BC>& BC_right,
        const double speed, const double dt,
        const std::vector<std::vector<double>>& psi_initial,
        std::vector<std::vector<double>>& phi,
        const std::string accelerator_type,
        const std::vector<double>& time );

#endif // _SOLVER_H

