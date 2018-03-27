#ifndef _SOLVER_H
#define _SOLVER_H

#include <vector>
#include <memory>

#include "Objects.h"
#include "Accelerator.h"

void source_iteration( const double epsilon,
                       const std::vector<std::shared_ptr<Region>>& mesh,
                       const std::vector<double>& mu, 
                       const std::vector<double>& w,
                       const std::shared_ptr<BC>& BC_left,
                       const std::shared_ptr<BC>& BC_right,
                       AcceleratorDSA& DSA,
                       std::vector<double>& phi,
                       std::vector<double>& rho );

void source_iteration_TD( const double epsilon,
                          const std::vector<std::shared_ptr<Region>>& mesh,
                          const std::vector<double>& mu, 
                          const std::vector<double>& w,
                          const std::shared_ptr<BC>& BC_left,
                          const std::shared_ptr<BC>& BC_right,
                          AcceleratorDSA& DSA, 
                          const double speed, const double dt, const int K,
                          const std::vector<std::vector<double>>& psi_initial,
                          std::vector<std::vector<double>>& phi );


#endif // _SOLVER_H

