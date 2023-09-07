
#include "Artificial_viscosity.hpp"
#include "utility.hpp"
#include <iostream>

using namespace std;



void dilation_based_viscosity_Euler(grid_vec_qq const& u, grid_point_qq & nu)
{
    #pragma omp parallel
    {

        ND::Array<order,order> speed1;
        ND::Array<order,order> speed2;

        double ds = 0.0;

        double temp1 = 0.0;
        double temp2 = 0.0;


        double v1 = 0;
        double v2 = 0;

        double P = 0;
        double a = 0;

        double df_max = 0.0;

        #pragma omp for collapse(2)
        for(size_t i=0; i<u.size(); i++)
        {
            for(size_t j=0; j<u.size(1); j++)
            {
                df_max = 0.0;

                for(int l = 0; l<order; l++)
                {
                    for(int m = 0; m<order; m++)
                    {
                        speed1(l,m) = u(i,j,l,m,1)/u(i,j,l,m,0);
                        speed2(l,m) = u(i,j,l,m,2)/u(i,j,l,m,0);


                        //=================================================//

                        v1 = u(i,j,l,m,1)/u(i,j,l,m,0);
                        v2 = u(i,j,l,m,2)/u(i,j,l,m,0);

                        P = (gamma_-1)*(u(i,j,l,m,3)-0.5*((v1*v1 + v2*v2)*u(i,j,l,m,0)));

                        a = sqrt(gamma_*P/u(i,j,l,m,0));

                        df_max = max(df_max,sqrt((abs(v1)+a)*(abs(v1)+a) + (abs(v2)+a)*(abs(v2)+a)));
                    }
                }
                

                for(int l = 0; l<order; l++)
                {
                    for(int m = 0; m<order; m++)
                    {
                        temp1 = 0.0;
                        temp2 = 0.0;
                        
                        for(int s = 0; s<order; s++)
                        {
                            temp1 += D(m,s)*speed1(l,s);
                        }

                        for(int r = 0; r<order; r++)
                        {
                            temp2 += D(l,r)*speed2(r,m);
                        }

                        // take the absolute value of the divergence
                        ds = abs(2*temp1/dx + 2*temp2/dy);

                        nu(i,j,l,m) = min(c_vis*ds*dx*dy/((order-1)*(order-1)), c_max*sqrt(dx*dy)*df_max/(order-1));
                    }
                }
            }
        }
    }
}





/*

void entropy_based_viscosity_Euler(grid_vec_qq const& u, grid_point_qq & nu)
{
    ND::Array<Ny,Nx> nu_av;

    ND::Array<Ny,Nx,order> nu_av_y;
    
    
//    #pragma omp parallel
//    {
//        Mat_2 vi(Mat_1(System_dim), order);
//        Mat_1 EC_temp(System_dim);
//        Mat_1 f_diff_0;
//        Mat_1 f_diff_k;
//        
//        double Q = 0.0;
//
//        double P_m_p;
//        double F_m_p;
//        double P_i_m;
//        double F_i_m;
//        double P_i_p;
//        double F_i_p;
//        double P_p_m;
//        double F_p_m;
//
//        double H_max;
//        double mu_E;
//        double v;
//        double P;
//        double a;
//        double df_max;
//
//        double mu_max;
//
//        #pragma omp for
//        for(size_t i=1; i<u.size()-1; i++)
//        {
//            vi = entropy_var(u(i));
//
//            
//            EC_flux(u[i-1][order-1],u(i,0), EC_temp);
//            f_diff_0 = int_flux(u[i-1][order-1],u(i,0))-EC_temp;
//
//            EC_flux(u(i)[order-1],u[i+1](0), EC_temp);
//            f_diff_k = int_flux(u(i)[order-1],u[i+1](0))-EC_temp;
//
//            Q = 0.0;
//
//            for(int i=0; i<order; i++)
//            {
//                Q += vi[order-1](i)*f_diff_k(i)-vi(0,i)*f_diff_0(i);
//            }
//
//            P_m_p = (gamma_-1)*(u[i-1][order-1](2)-0.5*(u[i-1][order-1](1)*u[i-1][order-1](1)/u[i-1][order-1](0)));
//            F_m_p = u[i-1][order-1](1)*log(P_m_p/pow(u[i-1][order-1](0), gamma_))/(gamma_-1);
//
//            P_i_m = (gamma_-1)*(u(i,0,2)-0.5*(u(i,0,1)*u(i,0,1)/u(i,0,0)));
//            F_i_m = u(i,0,1)*log(P_i_m/pow(u(i,0,0), gamma_))/(gamma_-1);
//
//            P_i_p = (gamma_-1)*(u(i)[order-1](2)-0.5*(u(i)[order-1](1)*u(i)[order-1](1)/u(i)[order-1](0)));
//            F_i_p = u(i)[order-1](1)*log(P_i_p/pow(u(i)[order-1](0), gamma_))/(gamma_-1);
//
//            P_p_m = (gamma_-1)*(u[i+1](0,2)-0.5*(u[i+1](0,1)*u[i+1](0,1)/u[i+1](0,0)));
//            F_p_m = u[i+1](0,1)*log(P_p_m/pow(u[i+1](0,0), gamma_))/(gamma_-1);
//
//
//            H_max = (order-1)*max(abs(F_m_p-F_i_m),abs(F_i_p-F_p_m))/dx;
//
//
//            mu_E = c_vis*dx*dx*max(Q,H_max)/((order-1)*(order-1));
//
//
//            v = 0.0;
//            P = 0.0;
//            a = 0.0;
//
//            df_max = 0.0;
//
//            for(int l=0; l<order; l++)
//            {
//                v = u(i,j,l,m,1)/u(i,j,l,m,0);
//
//                P = (gamma_-1)*(u(i,j,l,m,2)-0.5*((v*v)*u(i,j,l,m,0)));
//
//                a = sqrt(gamma_*P/u(i,j,l,m,0));
//
//                df_max = max(df_max, abs(v)+a);
//            }
//
//            mu_max = c_max*dx*df_max/(order-1);
//
//            nu_av(i) = min(mu_E, mu_max);
//        }
//
//        double slope;
//
//        #pragma omp for
//        for(size_t i=1; i<u.size()-2; i++)
//        {
//            slope = (nu_av[i+1]-nu_av(i))/dx;
//            for(int j=floor(order/2.0); j<order; j++)
//            {
//                nu(i,j) = nu_av(i) + (dx/2.0)*(x_ref(j)-0.0)*slope;
//            }
//            for(int j=0; j<floor(order/2.0); j++)
//            {
//                nu[i+1](j) = nu_av[i+1] - (dx/2.0)*(0.0-x_ref(j))*slope;
//            }
//        }
//    }
    
    
    
    #pragma omp parallel
    {
        ND::Array<order,order> E;
        ND::Array<order,order> E_old;
        ND::Array<order,order> F_1;
        ND::Array<order,order> F_2;
        ND::Array<order,order> F_old_1;
        ND::Array<order,order> F_old_2;

        double p, p_old;
        
        double df_num;
        double df_num_old;

        double df_num_1;
        double df_num_2;

        double df_num_old_1;
        double df_num_old_2;

        double Q_max = 0.0;

        double pm_R, pp_L;
        double Em_R, Ep_L;
        double F_1_m_R, F_1_p_L;
        double F_2_m_R, F_2_p_L;
        double H_L = 0.0;
        double H_R = 0.0;
        double H_max = 0.0;

        double mu_E = 0.0;

        double v1;
        double v2;

        double P;
        double a;
        double df_max = 0.0;
        
        double h = sqrt(dx*dy);

        double mu_max;

        #pragma omp for collapse(2)
        for(size_t i=1; i<u.size()-1; i++)
        {
            for(size_t j=1; j<u.size(1)-1; j++)
            {
                for(int l=0; l<order; l++)
                {
                    for(int m=0; m<order; m++)
                    {
                        p = (gamma_-1)*(u(i,j,l,m,3)-0.5*(u(i,j,l,m,1)*u(i,j,l,m,1)+u(i,j,l,m,2)*u(i,j,l,m,2))/u(i,j,l,m,0));
                        p_old = (gamma_-1)*(u_old(i,j,l,m,3)-0.5*(u_old(i,j,l,m,1)*u_old(i,j,l,m,1)+u_old(i,j,l,m,2)*u_old(i,j,l,m,2))/u_old(i,j,l,m,0));

                        E(l,m) = -u(i,j,l,m,0)*(log(p)-gamma_*log(u(i,j,l,m,0)))/(gamma_-1);
                        E_old(l,m) = -u_old(i,j,l,m,0)*(log(p_old)-gamma_*log(u_old(i,j,l,m,0)))/(gamma_-1);

                        F_1(l,m) = u(i,j,l,m,1)/u(i,j,l,m,0)*E(l,m);
                        F_old_1(l,m) = u_old(i,j,l,m,1)/u_old(i,j,l,m,0)*E_old(l,m);

                        F_2(l,m) = u(i,j,l,m,2)/u(i,j,l,m,0)*E(l,m);
                        F_old_2(l,m) = u_old(i,j,l,m,2)/u_old(i,j,l,m,0)*E_old(l,m);
                    }
                }
                
                df_max = 0.0;
                Q_max = 0.0;

                for(int l=0; l<order; l++)
                {
                    for(int m=0; m<order; m++)
                    {
                        df_num = 0.0;
                        df_num_old = 0.0;

                        df_num_1 = 0.0;
                        df_num_2 = 0.0;
                        df_num_old_1 = 0.0;
                        df_num_old_2 = 0.0;

                        for(int r=0; r<order; r++)
                        {
                            df_num_1 += D(l,r)*F_1(l,r);
                            df_num_old_1 += D(l,r)*F_old_1(l,r);
                        }
                        for(int s=0; s<order; s++)
                        {
                            df_num_2 += D(l,s)*F_2(s,m);
                            df_num_old_2 += D(l,s)*F_old_2(s,m);
                        }

                        df_num = 2.0*(df_num_1/dx + df_num_2/dy);
                        df_num_old = 2.0*(df_num_old_1/dx + df_num_old_2/dy);

                        Q_max = max(Q_max,abs((E(l,m)-E_old(l,m))/dt_old + 0.5*(df_num+df_num_old)));
                        //Q_max = max(Q_max,abs((E(l)-E_old(l))/dt_old + df_num));

                        // df_max

                        v1 = u(i,j,l,m,1)/u(i,j,l,m,0);
                        v2 = u(i,j,l,m,2)/u(i,j,l,m,0);

                        P = (gamma_-1)*(u(i,j,l,m,3)-0.5*((v1*v1+v2*v2)*u(i,j,l,m,0)));

                        a = sqrt(gamma_*P/u(i,j,l,m,0));

                        df_max = max(df_max, max(abs(v1)+a,abs(v2)+a));
                    }
                }
                
                H_max = 0.0;
                
                for(int l=0; l<order; l++)
                {
                    // vertical
                    pm_R = (gamma_-1)*(u(i)[j-1](l)[order-1](3)-0.5*(u(i)[j-1](l)[order-1](1)*u(i)[j-1](l)[order-1](1)+u(i)[j-1](l)[order-1](2)*u(i)[j-1](l)[order-1](2))/u(i)[j-1](l)[order-1](0));
                    pp_L = (gamma_-1)*(u(i,j+1,l,0,3)-0.5*(u(i,j+1,l,0,1)*u(i,j+1,l,0,1)+u(i,j+1,l,0,2)*u(i,j+1,l,0,2))/u(i,j+1,l,0,0));

                    Em_R = -u(i)[j-1](l)[order-1](0)*(log(pm_R)-gamma_*log(u(i)[j-1](l)[order-1](0)))/(gamma_-1);
                    Ep_L = -u(i,j+1,l,0,0)*(log(pp_L)-gamma_*log(u(i,j+1,l,0,0)))/(gamma_-1);


                    F_1_m_R = u(i)[j-1](l)[order-1](1)/u(i)[j-1](l)[order-1](0)*Em_R;
                    F_1_p_L = u(i,j+1,l,0,1)/u(i,j+1,l,0,0)*Ep_L;

                    H_L = -(order-1)*(F_1(l,0)-F_1_m_R)/h;
                    H_R = (order-1)*(F_1(l)[order-1]-F_1_p_L)/h;

                    H_max = max(H_max,max(abs(H_L),abs(H_R)));


                    // Horizontal
                    pm_R = (gamma_-1)*(u[i-1](j)[order-1](l,3)-0.5*((u[i-1](j)[order-1](l,1)*u[i-1](j)[order-1](l,1)+u[i-1](j)[order-1](l,2)*u[i-1](j)[order-1](l,2))/u[i-1](j)[order-1](l,0)));
                    pp_L = (gamma_-1)*(u(i+1,j,0,l,3)-0.5*((u(i+1,j,0,l,1)*u(i+1,j,0,l,1)+u(i+1,j,0,l,2)*u(i+1,j,0,l,2))/u(i+1,j,0,l,0)));

                    Em_R = -u[i-1](j)[order-1](l,0)*(log(pm_R)-gamma_*log(u[i-1](j)[order-1](l,0)))/(gamma_-1);
                    Ep_L = -u(i+1,j,0,l,0)*(log(pp_L)-gamma_*log(u(i+1,j,0,l,0)))/(gamma_-1);


                    F_2_m_R = u[i-1](j)[order-1](l,2)/u[i-1](j)[order-1](l,0)*Em_R;
                    F_2_p_L = u(i+1,j,0,l,2)/u(i+1,j,0,l,0)*Ep_L;

                    H_L = -(order-1)*(F_2(0,l)-F_2_m_R)/h;
                    H_R = (order-1)*(F_2[order-1](l)-F_2_p_L)/h;

                    H_max = max(H_max,max(abs(H_L),abs(H_R)));
                }
                
                mu_E = c_vis*max(Q_max, H_max)*(dx*dy/((order-1)*(order-1)));
                
                mu_max = c_max*dx*df_max/(order-1);

                nu_av(i,j) = min(mu_E, mu_max);
            }
        }

        double slope;
        
        //#pragma omp for collapse(2)
        //for(size_t i=1; i<u.size()-1; i++)
        //{
        //    for(size_t j=1; j<u.size(1)-1; j++)
        //    {
        //        for(int l=0; l<order; l++)
        //        {
        //            for(int m=0; m<order; m++)
        //            {
        //                nu(i,j,l,m) = nu_av(i,j);
        //            }
        //        }
        //    }
        //}
        #pragma omp for collapse(2)
        for(size_t i=1; i<u.size()-1; i++)
        {
            for(size_t j=1; j<u.size(1)-1; j++)
            {
                slope = 0.5*(nu_av[i+1](j)-nu_av[i-1](j))/dy;

                for(int l=0; l<order; l++)
                {
                    nu_av_y(i,j,l) = 0.5*(nu_av[i-1](j)+nu_av(i,j)) + (dy/2.0)*(x_ref(l)+1)*slope;
                }
            }
        }
        #pragma omp for collapse(2)
        for(size_t i=1; i<u.size()-1; i++)
        {
            for(size_t j=1; j<u.size(1)-1; j++)
            {
                for(int l=0; l<order; l++)
                {
                    slope = 0.5*(nu_av_y(i)[j+1](l)-nu_av_y(i)[j-1](l))/dx;

                    for(int m=0; m<order; m++)
                    {
                        nu(i,j,l,m) = 0.5*(nu_av_y[i-1](j,l)+nu_av_y(i,j,l)) + (dx/2.0)*(x_ref(m)+1)*slope;
                    }
                }
            }
        }
        #pragma omp for
        for(size_t i=0; i<u.size(); i++)
        {
            for(int l=0; l<order; l++)
            {
                for(int m=0; m<order; m++)
                {
                    nu(i,0,l,m) = nu_av(i,0);
                    nu(i)[u.size(1)-1](l,m) = nu_av(i)[u.size(1)-1];
                }
            }
        }
        #pragma omp for
        for(size_t j=0; j<u.size(1); j++)
        {
            for(int l=0; l<order; l++)
            {
                for(int m=0; m<order; m++)
                {
                    nu(0,j,l,m) = nu_av(0,j);
                    nu[u.size()-1](j,l,m) = nu_av[u.size()-1](j);
                }
            }
        }
    }
}

*/