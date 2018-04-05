#ifndef _ACCELERATOR_H_
#define _ACCELERATOR_H_

#include <memory>
#include <vector>

#include "Objects.h"


class Accelerator
{
    private:
        const std::string a_type;
    public:
        Accelerator( const std::string t ): a_type(t) {};
        ~Accelerator() {};


        virtual std::string type() { return a_type; }
        virtual void accelerate( const std::vector<std::shared_ptr<Region>>& mesh,
                                 const std::vector<double>& phi_old, 
                                 std::vector<double>& phi ) = 0;
};


//=============================================================================
// NONE
//=============================================================================

class AcceleratorNONE : public Accelerator
{
    private:
        int J;
        std::vector<double> A;
        std::vector<double> B;
        std::vector<double> C;
        std::vector<double> f;
    public:
        AcceleratorNONE(): Accelerator("NONE") {};
        ~AcceleratorNONE() {};

        void accelerate( const std::vector<std::shared_ptr<Region>>& mesh,
                         const std::vector<double>& phi_old, 
                               std::vector<double>& phi );
};
//=============================================================================
// DSA
//=============================================================================

class AcceleratorDSA : public Accelerator
{
    private:
        int J;
        std::vector<double> A;
        std::vector<double> B;
        std::vector<double> C;
        std::vector<double> f;
    public:
        AcceleratorDSA( const std::vector<std::shared_ptr<Region>>& mesh, 
                        std::shared_ptr<BC> BC_left,
                        std::shared_ptr<BC> BC_right );
        ~AcceleratorDSA() {};

        void accelerate( const std::vector<std::shared_ptr<Region>>& mesh,
                         const std::vector<double>& phi_old, 
                               std::vector<double>& phi );
};
//=============================================================================
// IDSA
//=============================================================================

class AcceleratorIDSA : public Accelerator
{
    private:
        int J;
        std::vector<double> A;
        std::vector<double> B;
        std::vector<double> C;
        std::vector<double> f;
    public:
        AcceleratorIDSA( const std::vector<std::shared_ptr<Region>>& mesh, 
                         std::shared_ptr<BC> BC_left,
                         std::shared_ptr<BC> BC_right );
        ~AcceleratorIDSA() {};

        void accelerate( const std::vector<std::shared_ptr<Region>>& mesh,
                         const std::vector<double>& phi_old, 
                               std::vector<double>& phi );
};

#endif // _ACCELERATOR_H
