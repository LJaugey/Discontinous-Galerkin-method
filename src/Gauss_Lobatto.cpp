#define _USE_MATH_DEFINES

#include "Gauss_Lobatto.hpp"
#include <string>
#include <cmath>
#include <limits>
#include <iostream>

#include "write_data.hpp"
#include "utility.hpp"



using namespace std;

grid_vec_qq Gauss_Lobatto_2D(Array<Nx+1> x, Array<Ny+1> y, Array<System_dim> (*f)(double,std::string, double,std::string))
{
    std::string boundary_1 = "None";
    std::string boundary_2 = "None";

    Array<order> x_;
    Array<order> y_;
    
    grid_vec_qq u;

    Array<System_dim> f_temp;


    for(int i = 0; i<Ny; i++)
    {
        y_ = lgl_nodes();

        y_ = y(i) + 0.5*(y_-y_(0))*dy;
        

        for(int j = 0; j<Nx; j++)
        {
            x_ = lgl_nodes();
        
            x_ = x(j) + 0.5*(x_-x_(0))*dx;

            for(int m = 0; m<order; m++)
            {
                if(m == 0)              boundary_2 = "start";
                else if(m == order-1)   boundary_2 = "end";
                else                    boundary_2 = "None";

                for(int l = 0; l<order; l++)
                {
                    if(l == 0)              boundary_1 = "start";
                    else if(l == order-1)   boundary_1 = "end";
                    else                    boundary_1 = "None";

                    f_temp = f(x_(l),boundary_1 ,y_(m),boundary_2);
                    
                    for(int k = 0; k<System_dim; k++)
                    {
                        u(i+1,j+1,m,l,k) = f_temp(k);
                    }
                }
            }
            if((Save_data || Save_data_last) && i == 0)
            {
                write_data(x_, "x"+exp_name);
            }
        }
        if(Save_data || Save_data_last)
        {
            write_data(y_, "y"+exp_name);
        }
    }
    
    return u;
}








// The next two functions are inspired by the matlab function found here:
// https://www.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights
Array<order> lgl_nodes()
{
    // Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    Array<order> x;

    for(int i = 0; i<order;i++)
    {
        x(i) = cos(M_PI*i/(order-1));
    }

    // The Legendre Vandermonde Matrix
    
    Array<order,order> P;

    // Compute P_((order-1)) using the recursion relation
    // Compute its first and second derivatives and 
    // update x using the Newton-Raphson method.

    Array<order> xold;
    xold = 2.0;
    
    while((abs(x-xold).max())>numeric_limits<double>::epsilon())
    {
        xold=x;

        P[0] = 1;
        P[1] = x;

        for(int k=2;k<order;k++)
        {
            P[k]=( (2*k-1)*x*P[k-1]-(k-1)*P[k-2] )/k;
        }

        x = xold - ( x*P[order-1]-P[order-2] )/( order*P[order-1] );
    }

    return -x;
}



Array<order> lgl_weights(Array<order> x)
{
    Array<order,order> P;

    P[0].fill(1);
    P[1].fill(x);
        
    for(int k=2;k<order;k++)
    {
        P[k]=( (2*k-1)*x*P[k-1]-(k-1)*P[k-2] )/k;
    }

    return 2.0/((order-1)*order*P[order-1]*P[order-1]);
}


