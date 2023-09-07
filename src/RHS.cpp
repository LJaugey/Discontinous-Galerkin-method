#include "RHS.hpp"
#include "utility.hpp"

#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

void RHS_2D(const grid_vec_qq u,grid_vec_qq v,grid_vec_qq result, const grid_point_qq nu, grid_vec_qq q1,grid_vec_qq q2, grid_vec_q v_hat_im, grid_vec_q v_hat_ip, grid_vec_q v_hat_mj, grid_vec_q v_hat_pj, grid_vec_ EC_temp, grid_vec_ dvm, grid_vec_ dvp, grid_vec_ qm_p, grid_vec_ q_hat_m, grid_vec_ f1_, grid_vec_ fs_m, grid_vec_ qp_m, grid_vec_ q_hat_p, grid_vec_ fs_p, grid_vec_ f2_)
{
    #pragma omp parallel
    {
        #pragma omp for collapse(2)
        for(size_t i = 0; i<result.size(); i++)
        {
            for(size_t j = 0; j<result.size(1); j++)
            {
                v[i][j] = entropy_var(u[i][j]);
            }
        }

        #pragma omp for collapse(2)
        for(size_t i = 1; i<result.size()-1; i++)  //do not take first/last cell because of BCs
        {
            for(size_t j = 1; j<result.size(1)-1; j++)  //do not take first/last cell because of BCs
            {
                RHS_2D_cell(result[i][j], u[i-1][j],u[i][j-1],u[i][j],u[i][j+1],u[i+1][j],v[i-1][j],v[i][j-1],v[i][j],v[i][j+1],v[i+1][j], nu[i][j], q1[i][j],q2[i][j], v_hat_im[i][j], v_hat_ip[i][j], v_hat_mj[i][j], v_hat_pj[i][j], EC_temp[i][j], dvm[i][j],dvp[i][j],qm_p[i][j],q_hat_m[i][j],f1_[i][j],fs_m[i][j],qp_m[i][j],q_hat_p[i][j],fs_p[i][j],f2_[i][j]);
            }
        }
    }
}

void RHS_2D_cell(vec_qq result, const vec_qq umj, const vec_qq uim, const vec_qq uij, const vec_qq uip, const vec_qq upj, const vec_qq vmj, const vec_qq vim, const vec_qq vij, const vec_qq vip, const vec_qq vpj, const point_qq nu, vec_qq q1,vec_qq q2, vec_q v_hat_im, vec_q v_hat_ip, vec_q v_hat_mj, vec_q v_hat_pj, vec_ EC_temp, vec_ dvm, vec_ dvp, vec_ qm_p, vec_ q_hat_m, vec_ f1_, vec_ fs_m, vec_ qp_m, vec_ q_hat_p, vec_ fs_p, vec_ f2_)
{
    result.fill(0.0);
    q1.fill(0.0);
    q2.fill(0.0);
    dvm.fill(0.0);
    dvp.fill(0.0);


    // Horizontal

    for(int l=0; l<order; l++)
    {
        for(int k=0; k<System_dim; k++)
        {
            // left
            v_hat_im(l,k) = 0.5*(vim(l,order-1,k)+vij(l,0,k));
            // right
            v_hat_ip(l,k) = 0.5*(vij(l,order-1,k)+vip(l,0,k));
        }
    }

    // Vertical

    for(int m=0; m<order; m++)
    {
        for(int k=0; k<System_dim; k++)
        {
            // lower
            v_hat_mj(m,k) = 0.5*(vmj(order-1,m,k)+vij(0,m,k));
            // upper
            v_hat_pj(m,k) = 0.5*(vij(order-1,m,k)+vpj(0,m,k));
        }
    }
    

    
    
    if(viscosity_bool)
    {
        for(int l=0; l<order; l++)
        {
            for(int m=0; m<order; m++)
            {
                for(int s=0; s<order; s++)
                {
                    q1[l][m] += D(m,s)*vij[l][s];
                }

                for(int r=0; r<order; r++)
                {
                    q2[l][m] += D(l,r)*vij[r][m];
                }
                
                q1[l][m] *= nu(l,m)*2.0/dx;
                q2[l][m] *= nu(l,m)*2.0/dy;
            }
        }
    

        // Horizontal
        for(int l=0; l<order; l++)
        {
            //left boundary
            q1[l][0] += nu(l,0)*2*(vij[l][0]-v_hat_im[l])/(dx*w_ref[0]);
            //right boundary
            q1[l][order-1] += nu[l][order-1]*2*(v_hat_ip[l]-vij[l][order-1])/(dx*w_ref[order-1]);
        }

        

        // Vertical
        for(int m=0; m<order; m++)
        {
            //lower boundary
            q2[0][m] += nu(0,m)*2*(vij[0][m]-v_hat_mj[m])/(dy*w_ref[0]);
            //upper boundary
            q2[order-1][m] += nu[order-1][m]*2*(v_hat_pj[m]-vij[order-1][m])/(dy*w_ref[order-1]);
        }
        // q done

        for(int l=0; l<order; l++)
        {
            for(int m=0; m<order; m++)
            {
                for(int s=0; s<order; s++)
                {
                    result[l][m] += (2.0/dx)*D(m,s)*q1[l][s];
                }
                for(int r=0; r<order; r++)
                {
                    result[l][m] += (2.0/dy)*D(l,r)*q2[r][m];
                }
            }
        }
    }


    
    for(int l=0; l<order; l++)
    {
        for(int m=0; m<order; m++)
        {
            EC_flux1(uij[l][m],uij[l][m], EC_temp);
            
            result[l][m] -= (4.0/dx)*D(m,m)*EC_temp;

            for(int s=m+1; s<order; s++)
            {
                EC_flux1(uij[l][m],uij[l][s], EC_temp);

                result[l][m] -= (4.0/dx)*D(m,s)*EC_temp;
                result[l][s] -= (4.0/dx)*D(s,m)*EC_temp;
            }
            


            EC_flux2(uij[l][m],uij[l][m], EC_temp);
            
            result[l][m] -= (4.0/dy)*D(l,l)*EC_temp;

            for(int r=l+1; r<order; r++)
            {
                EC_flux2(uij[l][m],uij[r][m], EC_temp);

                result[l][m] -= (4.0/dy)*D(l,r)*EC_temp;
                result[r][m] -= (4.0/dy)*D(r,l)*EC_temp;
            }
        }
    }
    
    
    for(int l=0; l<order; l++)
    {
        for(int k = 0; k<System_dim; k++)
        {
            dvm[k] = 0.0;
            dvp[k] = 0.0;
        }

        // left boundary
        for(int s=0; s<order; s++)
        {
            dvm += D[order-1][s]*vim[l][s];
        }
        qm_p = nu(l,0)*(2.0/dx)*(dvm + (v_hat_im[l]-vim[l][order-1])/w_ref[order-1]);
        q_hat_m = 0.5*(qm_p+q1[l][0]);
        f1_ = f1(uij[l][0]);
        fs_m = int_flux1(uim[l][order-1],uij[l][0]);
        result[l][0] += (2.0/dx)*(fs_m + q1[l][0] - f1_ - q_hat_m)/w_ref[0];

        
        // right boundary
        for(int r=0; r<order; r++)
        {
            dvp += D(0,r)*vip[l][r];
        }
        qp_m = nu[l][order-1]*(2.0/dx)*(dvp + (vip[l][0]-v_hat_ip[l])/w_ref[0]);
        q_hat_p = 0.5*(q1[l][order-1]+qp_m);
        f1_ = f1(uij[l][order-1]);
        fs_p = int_flux1(uij[l][order-1],uip[l][0]);
        result[l][order-1] += (2.0/dx)*(f1_ + q_hat_p - fs_p - q1[l][order-1])/w_ref[order-1];
    }


    for(int m=0; m<order; m++)
    {
        for(int k = 0; k<System_dim; k++)
        {
            dvm[k] = 0.0;
            dvp[k] = 0.0;
        }
        
        // upper boundary
        for(int r=0; r<order; r++)
        {
            dvm += D[order-1][r]*vmj[r][m];
        }
        qm_p = nu(0,m)*(2.0/dy)*(dvm + (v_hat_mj[m]-vmj[order-1][m])/w_ref[order-1]);
        q_hat_m = 0.5*(qm_p+q2[0][m]);
        f2_ = f2(uij[0][m]);
        fs_m = int_flux2(umj[order-1][m],uij[0][m]);
        result[0][m] += (2.0/dy)*(fs_m + q2[0][m] - f2_ - q_hat_m)/w_ref[0];
    
    
        // lower boundary
        for(int r=0; r<order; r++)
        {
            dvp += D(0,r)*vpj[r][m];
        }
        qp_m = nu[order-1][m]*(2.0/dy)*(dvp + (vpj[0][m]-v_hat_pj[m])/w_ref[0]);
        q_hat_p = 0.5*(q2[order-1][m]+qp_m);
        f2_ = f2(uij[order-1][m]);
        fs_p = int_flux2(uij[order-1][m],upj[0][m]);
        result[order-1][m] += (2.0/dy)*(f2_ + q_hat_p - fs_p - q2[order-1][m])/w_ref[order-1];
    }
}
