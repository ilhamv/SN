#ifndef _ACCELERATOR_H_
#define _ACCELERATOR_H_

#include <memory>
#include <vector>

#include "Objects.h"


class AcceleratorDSA
{
    private:
        int J;
        std::vector<double> A;
        std::vector<double> B;
        std::vector<double> C;
        std::vector<double> f;
    public:
        AcceleratorDSA( const std::vector<std::shared_ptr<Region>>& mesh, 
                        const std::shared_ptr<BC>& BC_left,
                        const std::shared_ptr<BC>& BC_right );
        ~AcceleratorDSA() {};

        void accelerate( const std::vector<std::shared_ptr<Region>>& mesh,
                         const std::vector<double>& phi_old, 
                               std::vector<double>& phi );
};

#endif // _ACCELERATOR_H
