#ifndef _OBJECTS_H
#define _OBJECTS_H

#include <vector>
#include <string>
#include <memory>


//==============================================================================
// Material
//==============================================================================

class Material
{
    private:
        const int         m_id;
        const std::string m_name;
        const double      m_SigmaS;
        const double      m_SigmaT;
    public:
        Material( const int i, const std::string n, const double t, 
                  const double s ): 
            m_id(i), m_name(n), m_SigmaS(s), m_SigmaT(t) {};
        ~Material() {};
        int id() { return m_id; }
        std::string name() { return m_name; }
        double SigmaT() { return m_SigmaT; }
        double SigmaS() { return m_SigmaS; }
};


//==============================================================================
// Region
//==============================================================================

class Region
{
    private:
        const std::shared_ptr<Material> M;
        const double r_dz;
        const double r_Q;

    public:
        Region( const std::shared_ptr<Material>& m, const double dz,
                const double Q ): M(m), r_dz(dz), r_Q(Q) {};
        ~Region() {};

        double SigmaT() { return M->SigmaT(); }
        double SigmaS() { return M->SigmaS(); }
        double dz() { return r_dz; }
        double Q() { return r_Q; }
        std::shared_ptr<Material> material() { return M; }
};


//==============================================================================
// Boundary Condition
//==============================================================================

class BC
{
    private:
        const std::string bc_type;
    public:
        BC( const std::string t ): bc_type(t) {};
        ~BC() {};
        virtual void set_boundary( std::vector<double>& psi ) = 0;
        virtual std::string type() final { return bc_type; }
};

class BCVacuum : public BC
{
    public:
        BCVacuum(): BC("Vacuum") {};
        ~BCVacuum() {};
        void set_boundary( std::vector<double>& psi );
};

class BCReflective : public BC
{
    public:
        BCReflective(): BC("Reflective") {};
        ~BCReflective() {};
        void set_boundary( std::vector<double>& psi );
};

class BCIsotropic : public BC
{
    private:
        const double magnitude;

    public:
        BCIsotropic( const double m ): BC("Isotropic"), magnitude(m) {};
        ~BCIsotropic() {};
        void set_boundary( std::vector<double>& psi );
};

class BCMonoDirectional : public BC
{
    private:
        const std::vector<double> psi_b;

    public:
        BCMonoDirectional( const std::vector<double>& p ): 
            BC("Isotropic"), psi_b(p) {};
        ~BCMonoDirectional() {};
        void set_boundary( std::vector<double>& psi );
};

#endif // OBJECTS_H
