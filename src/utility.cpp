#include "utility.hpp"
#include "omp.h"

#include <iostream>

using namespace std;


Array<Ny,Nx, System_dim> quadrature_average_2D(grid_vec_qq const& u)
{
    Array<Ny,Nx, System_dim> u_bar;

    #pragma omp parallel for collapse(2) if(omp_get_num_threads() == 1)
    for(size_t i=0; i<Ny; i++)
    {
        for(size_t j=0; j<Nx; j++)
        {
            for(int k=0; k<System_dim; k++)
            {
                for(int l=0; l<order; l++)
                {
                    for(int m=0; m<order; m++)
                    {
                        u_bar(i,j,k) += w_ref(l)*w_ref(m)*u(i+1,j+1,l,m,k);
                    }
                }
                u_bar(i,j,k) *= 0.25;
            }
        }
    }

    return u_bar;
}


Array<Ny,Nx> quadrature_average_2D(grid_point_qq const& nu)
{
    Array<Ny,Nx> nu_bar;

    #pragma omp parallel for collapse(2) if(omp_get_num_threads() == 1)
    for(size_t i=0; i<Ny; i++)
    {
        for(size_t j=0; j<Nx; j++)
        {
            for(int l=0; l<order; l++)
            {
                for(int m=0; m<order; m++)
                {
                    nu_bar(i,j) += w_ref(l)*w_ref(m)*nu(i+1,j+1,l,m);
                }
            }
            nu_bar(i,j) *= 0.25;
        }
    }

    return nu_bar;
}