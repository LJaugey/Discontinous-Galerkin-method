#include <iostream>
#include <cmath>
#include <string>
#include <cstring>
#include <chrono>
#include <omp.h>

//#include "parameters.hpp"
//#include "Artificial_viscosity.hpp"
#include "utility.hpp"
//#include "init.hpp"
#include "initial_condition.hpp"
//#include "Equation.hpp"
//#include "entropy_variables.hpp"
#include "Gauss_Lobatto.hpp"
#include "RHS.hpp"
#include "write_data.hpp"



using namespace std;


void apply_BC_2D(grid_vec_qq& u, string BC)
{
    if(BC == "outflow")
    {
        // Horizontal
        for(size_t i=1; i<u.size()-1; i++)
        {
            for(int l=0; l<order; l++)
            {
                for(int m=0; m<order; m++)
                {
                    for(int k=0; k<System_dim; k++)
                    {
                        // left
                        u(i,0,l,m, k) = u(i,1,l,0, k);
                        // right
                        u(i,u.size(1)-1,l,m, k) = u(i,u.size(1)-2,l,order-1, k);
                    }
                }
            }
        }
        // Vertical
        for(size_t j=1; j<u.size(1)-1; j++)
        {
            for(int l=0; l<order; l++)
            {
                for(int m=0; m<order; m++)
                {
                    for(int k=0; k<System_dim; k++)
                    {
                        // lower (before: upper)
                        u(0,j,l,m,k) = u(1,j,0,m,k);
                        // upper (before: lower)
                        u(u.size()-1,j,l,m,k) = u(u.size()-2,j,order-1,m,k);
                    }
                }
            }
        }
    }
    else if(BC == "periodic")
    {
        // Horizontal
        for(size_t i=1; i<u.size()-1; i++)
        {
            for(int l=0; l<order; l++)
            {
                for(int m=0; m<order; m++)
                {
                    for(int k=0; k<System_dim; k++)
                    {
                        // left
                        u(i,0,l,m,k) = u(i,u.size(1)-2,l,order-1,k);
                        // right
                        u(i,u.size(1)-1,l,m,k) = u(i,1,l,0,k);
                    }
                }
            }
        }
        // Vertical
        for(size_t j=1; j<u.size(1)-1; j++)
        {
            for(int l=0; l<order; l++)
            {
                for(int m=0; m<order; m++)
                {
                    for(int k=0; k<System_dim; k++)
                    {
                        // lower (before: upper)
                        u(0,j,l,m,k) = u(u.size()-2,j,order-1,m,k);
                        // upper (before: lower)
                        u(u.size()-1,j,l,m,k) = u(1,j,0,m,k);
                    }
                }
            }
        }
    }
}


double compute_dt(grid_vec_qq const& u, grid_point_qq const& nu)
{
    if(Equation_Name=="Custom")
    {
        return t_fin;
    }
    else if(Equation_Name=="Euler")
    {
        double dt_inv = 0.0;
        double v1;
        double v2;
        double p;
        double a;

        double Max_eig1 = 0.0;
        double Max_eig2 = 0.0;
        
        //#pragma omp parallel for private(v1,v2, p,a, Max_eig1,Max_eig2) reduction(max:dt_inv)
        #pragma omp parallel for collapse(4) private(v1,v2, p,a, Max_eig1,Max_eig2) reduction(max:dt_inv)
        for(size_t i=1; i<u.size()-1; i++)
        {
            for(size_t j=1; j<u.size(1)-1; j++)
            {
                for(int l=0; l<order; l++)
                {
                    for(int m=0; m<order; m++)
                    {
                        v1 = u(i,j,l,m,1)/u(i,j,l,m,0);
                        v2 = u(i,j,l,m,2)/u(i,j,l,m,0);

                        p = (gamma_-1)*(u(i,j,l,m,3)-0.5*((v1*v1 + v2*v2)*u(i,j,l,m,0)));
                        a = sqrt(gamma_*p/u(i,j,l,m,0));

                        Max_eig1 = abs(v1) + a;
                        Max_eig2 = abs(v2) + a;

                        dt_inv = max(dt_inv, pow(order-1, 2)*0.5*(Max_eig1/dx + Max_eig2/dy) + nu(i,j,l,m)*pow(order-1,4)/(dx*dy));
                    }
                }
            }
        }
        return CFL/dt_inv;
    }
    else
    {
        Array<Ny,Nx, System_dim> u_bar = quadrature_average_2D(u);

        double max_val = u_bar(0,0,0);

        for(size_t i=0; i<u_bar.size(); i++)
        {
            for(size_t j=0; j<u_bar.size(1); j++)
            {
                max_val = max(max_val, u_bar(i,j,0));
            }
        }
        return t_fin;
    }
    
    cerr<<"Warning: Equation name is not recognised. Stopping simulation..."<<endl;
    return t_fin;
}









int main(int argc, char **argv)
{
    double dt = t_fin;
    bool fixed_dt = false;
    /*
    if(argc>1)
    {
        if(strcmp(argv(1), "Nx") == 0)
        {
            Nx = stoi(argv(2));
            dx = (bx-ax)/Nx;
            exp_name = argv[argc-1];
            if(Artificial_Viscosity_Name!="Entropy")    dt = CFL*pow(dx,order/3.0);
            else                                        dt = CFL*pow(dx,order/2.0);
            fixed_dt = true;
        }
        else if(strcmp(argv(1), "Ny") == 0)
        {
            Ny = stoi(argv(2));
            dy = (by-ay)/Ny;
            exp_name = argv[argc-1];
            if(Artificial_Viscosity_Name!="Entropy")    dt = CFL*pow(dy,order/3.0);
            else                                        dt = CFL*pow(dy,order/2.0);
            fixed_dt = true;
        }
        if(strcmp(argv(3), "Nx") == 0)
        {
            Nx = stoi(argv(4));
            dx = (bx-ax)/Nx;
            exp_name = argv[argc-1];
            if(Artificial_Viscosity_Name!="Entropy")    dt = CFL*pow(min(dx,dy),order/3.0);
            else                                        dt = CFL*pow(min(dx,dy),order/2.0);
            fixed_dt = true;
        }
        else if(strcmp(argv(3), "Ny") == 0)
        {
            Ny = stoi(argv(4));
            dy = (by-ay)/Ny;
            exp_name = argv[argc-1];
            if(Artificial_Viscosity_Name!="Entropy")    dt = CFL*pow(min(dx,dy),order/3.0);
            else                                        dt = CFL*pow(min(dx,dy),order/2.0);
            fixed_dt = true;
        }
        else
        {
            fixed_dt = false;
        }
    }
    */
    
    cout<<"dx: "<<dx<<endl;
    cout<<"dy: "<<dy<<endl;



    if(Save_data || Save_data_last)
    {
        clear_data("x"+exp_name);
        clear_data("y"+exp_name);
        clear_data("t"+exp_name);
        clear_data("u"+exp_name);
        clear_data("nu"+exp_name);
    
        if(omp_get_max_threads()>1)
        {
            cerr<<"Warning: Artificial viscosity will not be saved (number of threads != 1)"<<endl;
        }
    }

    //init_Equation();

    Array<Nx+1> x_;
    Array<Ny+1> y_;

    for(int j = 0; j<Nx+1; j++)   x_(j) = ax + j*dx;
    for(int i = 0; i<Ny+1; i++)   y_(i) = ay + i*dy;

    // u[Ny,Nx, order,order, System_dim]


    t = t0;

    
    auto u0 = init_u0();

    grid_vec_qq u_0;
    grid_vec_qq u_1;
    grid_vec_qq u_2;
    grid_vec_qq q1;
    grid_vec_qq q2;
    grid_vec_qq v;
    grid_vec_qq res_temp;

    
    grid_vec_q v_hat_im;
    grid_vec_q v_hat_ip;
    grid_vec_q v_hat_mj;
    grid_vec_q v_hat_pj;
    
    grid_point_qq nu;


    grid_vec_ EC_temp;
    grid_vec_ dvm;
    grid_vec_ dvp;
    grid_vec_ qm_p;
    grid_vec_ q_hat_m;
    grid_vec_ f1_;
    grid_vec_ fs_m;
    grid_vec_ qp_m;
    grid_vec_ q_hat_p;
    grid_vec_ fs_p;
    grid_vec_ f2_;


    u_0 = Gauss_Lobatto_2D(x_, y_, u0);
    apply_BC_2D(u_0, BC);

    if(Save_data_first && !(Save_data && sampling==1))
    {
        if(fixed_dt)    write_data(u_0, "u"+exp_name);
        else            write_data(quadrature_average_2D(u_0), "u");
    
        write_data(t,"t"+exp_name);
    }

    double Save_it = false;

    int it = 0;
    
    int percent = 1;
    
    auto start = std::chrono::high_resolution_clock::now();


    u_old = u_0;

    while(t<t_fin)
    {
        dt_old = dt;
        
        if(viscosity_bool)
        {
            artificial_viscosity(u_0, nu);
            u_old = u_0;
        }

        if(!fixed_dt)   dt = compute_dt(u_0, nu);

        // save data if it is a sampling iteration or if it is the last iteration
        Save_it = (((it%sampling == sampling-1) && Save_data)) || (t+dt>=t_fin && Save_data_last);

        if(t+dt>=t_fin)
        {
            dt = t_fin-t;
            if(Save_data_last)  write_data(quadrature_average_2D(nu), "nu"+exp_name);
        }

        // u1
        RHS_2D(u_0,v,res_temp, nu, q1,q2, v_hat_im,v_hat_ip,v_hat_mj,v_hat_pj, EC_temp,dvm,dvp,qm_p,q_hat_m,f1_,fs_m,qp_m,q_hat_p,fs_p,f2_);
        u_1 = u_0 + dt*res_temp;
        apply_BC_2D(u_1, BC);



        // u2
        RHS_2D(u_1,v,res_temp, nu, q1,q2, v_hat_im,v_hat_ip,v_hat_mj,v_hat_pj, EC_temp,dvm,dvp,qm_p,q_hat_m,f1_,fs_m,qp_m,q_hat_p,fs_p,f2_);
        u_2 = 0.25*u_1 + 0.75*u_0 + (0.25*dt)*res_temp;
        apply_BC_2D(u_2, BC);

        // u
        RHS_2D(u_2,v,res_temp, nu, q1,q2, v_hat_im,v_hat_ip,v_hat_mj,v_hat_pj, EC_temp,dvm,dvp,qm_p,q_hat_m,f1_,fs_m,qp_m,q_hat_p,fs_p,f2_);
        u_0 = (1.0/3.0)*u_0 + (2.0/3.0)*u_2 + (2.0*dt/3.0)*res_temp;
        apply_BC_2D(u_0, BC);

        t += dt;

        if(floor(100*t/t_fin)>percent || t>= t_fin)
        {
            cout<<"iteration: "<<it<<"\t\t"<<"dt: "<<dt<<"\t\t"<<percent<<" %"<<endl;
            percent+=1;
        }
        if(Save_it)
        {
            if(fixed_dt)    write_data(u_0, "u"+exp_name);
            else            write_data(quadrature_average_2D(u_0), "u");
        
            write_data(t,"t"+exp_name);
        }
        
        it++;
    }

    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    

    cout<<endl<<"Code ran in "<<(double)duration.count()/1000.0<<" seconds"<<endl;

    return 0;
}
