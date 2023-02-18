#include "entropy_variables.hpp"

#include <cmath>

using namespace std;




vec_qq entropy_var_Custom(vec_qq const& u)
{
    return u;
}


vec_qq entropy_var_Transport(vec_qq const& u)
{
    return u;
}


vec_qq entropy_var_Burger(vec_qq const& u)
{
    return u;
}


vec_qq entropy_var_Euler(vec_qq const& u)
{
    double rho;
    double p;
    double s;

    vec_qq v;

    for(int l = 0; l < order; l++)
    {
        for(int m = 0; m < order; m++)
        {
            rho = u(l,m,0);

            p = (gamma_-1)*(u(l,m,3)-0.5*((u(l,m,1)*u(l,m,1)+u(l,m,2)*u(l,m,2))/u(l,m,0)));

            s = log(p) - gamma_*log(rho);

            v(l,m,0) = (gamma_-s)/(gamma_-1) - (u(l,m,1)*u(l,m,1)+u(l,m,2)*u(l,m,2))/(2*u(l,m,0)*p);
            v(l,m,1) = u(l,m,1)/p;
            v(l,m,2) = u(l,m,2)/p;
            v(l,m,3) = -u(l,m,0)/p;
        }
    }
    return v;
}
