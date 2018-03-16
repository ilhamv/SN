#include <vector>
#include <string>

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
