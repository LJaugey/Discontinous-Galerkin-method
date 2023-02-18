#ifndef RHS_HPP
#define RHS_HPP

#include "parameters.hpp"

void RHS_2D(const grid_vec_qq u,grid_vec_qq v,grid_vec_qq result, const Array<Ny+2,Nx+2, order,order> nu, grid_vec_qq q1,grid_vec_qq q2, grid_vec_q v_hat_im, grid_vec_q v_hat_ip, grid_vec_q v_hat_mj, grid_vec_q v_hat_pj, grid_vec_ EC_temp, grid_vec_ dvm, grid_vec_ dvp, grid_vec_ qm_p, grid_vec_ q_hat_m, grid_vec_ f1_, grid_vec_ fs_m, grid_vec_ qp_m, grid_vec_ q_hat_p, grid_vec_ fs_p, grid_vec_ f2_);
void RHS_2D_cell(vec_qq result, const vec_qq umj, const vec_qq uim, const vec_qq uij, const vec_qq uip, const vec_qq upj, const vec_qq vmj, const vec_qq vim, const vec_qq vij, const vec_qq vip, const vec_qq vpj, const point_qq nu, vec_qq q1,vec_qq q2, vec_q v_hat_im, vec_q v_hat_ip, vec_q v_hat_mj, vec_q v_hat_pj, vec_ EC_temp, vec_ dvm, vec_ dvp, vec_ qm_p, vec_ q_hat_m, vec_ f1_, vec_ fs_m, vec_ qp_m, vec_ q_hat_p, vec_ fs_p, vec_ f2_);

#endif