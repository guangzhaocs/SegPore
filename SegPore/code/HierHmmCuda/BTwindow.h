#ifndef BTWINDOW_H_INCLUDED
#define BTWINDOW_H_INCLUDED
#include <vector>
#include "Parameters.h"
using namespace std;


class BTwindow
{
   public:
      BTwindow(int _win_type, int _st_b_ind, int _en_b_ind, Parameters& _para):
          type(_win_type),
          st_b_ind(_st_b_ind),
          en_b_ind(_en_b_ind),
          global_para(_para)
          {}

      int get_sig_st(void);
      int get_sig_en(void);
      int get_sig_len(void);
      int get_b_len(void);
      float get_B_mu(void);
      float get_B_sigma(void);
      float get_B_hmm_loglik(void);
      float get_T_lr_loglik(void);
      void set_state_arr(void);

    //private:
      int type;
      int st_b_ind;
      int en_b_ind;
      Parameters& global_para;

};

int BTwindow::get_sig_st(void){
    if(type == 0){
        return global_para.b_st_sig_ind_arr[st_b_ind];
    }else{
        return global_para.t_st_sig_ind_arr[st_b_ind];
    }
}


int BTwindow::get_sig_en(void){
    if(type == 0){
        return global_para.b_en_sig_ind_arr[en_b_ind - 1];
    }else{
        return global_para.t_en_sig_ind_arr[en_b_ind - 1];
    }
}


int BTwindow::get_sig_len(void){
    if(type == 0){
        return global_para.b_en_sig_ind_arr[en_b_ind - 1] - global_para.b_st_sig_ind_arr[st_b_ind];
    }else{
        return global_para.t_en_sig_ind_arr[en_b_ind - 1] - global_para.t_st_sig_ind_arr[st_b_ind];
    }
}


int BTwindow::get_b_len(void){
    return en_b_ind - st_b_ind;
}


float BTwindow::get_B_mu(void){
    assert(type == 0);
    return global_para.B_mu_mat[st_b_ind][en_b_ind - 1 - st_b_ind];
}


float BTwindow::get_B_sigma(void){
    assert(type == 0);
    return global_para.B_sigma_mat[st_b_ind][en_b_ind - 1 - st_b_ind];
}


float BTwindow::get_B_hmm_loglik(void){
    assert(type == 0);
    return global_para.B_hmm_loglik_mat[st_b_ind][en_b_ind - 1 - st_b_ind];
}


float BTwindow::get_T_lr_loglik(void){
    assert(type == 1);
    return global_para.T_lr_loglik_arr[st_b_ind];
}


void BTwindow::set_state_arr(void){
    int state;
    if(type == 0){
        state = 0;
    }else{
        state = 1;
    }
    int sig_st = get_sig_st();
    int sig_ed = get_sig_en();
    for(int i=sig_st;i<sig_ed;i++){
        global_para.state_arr[i] = state;
    }
}


#endif // BTWINDOW_H_INCLUDED
