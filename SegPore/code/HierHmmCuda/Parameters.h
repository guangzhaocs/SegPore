#ifndef PARAMETERS_H_INCLUDED
#define PARAMETERS_H_INCLUDED
#include <vector>
using namespace std;
/////////////////////////////////////////////////////////////////
class Parameters
{
   public:
       Parameters(int _M, vector<float>& _signal_arr, vector<int>& _state_arr, vector<vector<float>>& _trans_mat,
                  int _N_B, vector<int>& _b_st_sig_ind_arr, vector<int>& _b_en_sig_ind_arr,
                  vector<vector<float>>& _B_mu_mat, vector<vector<float>>& _B_sigma_mat, vector<vector<float>>& _B_hmm_loglik_mat,
                  int _N_T, vector<int>& _t_st_sig_ind_arr, vector<int>& _t_en_sig_ind_arr, vector<float>& _T_lr_loglik_arr):
                  M(_M),
                  signal_arr(_signal_arr),
                  state_arr(_state_arr),
                  trans_mat(_trans_mat),
                  N_B(_N_B),
                  b_st_sig_ind_arr(_b_st_sig_ind_arr),
                  b_en_sig_ind_arr(_b_en_sig_ind_arr),
                  B_mu_mat(_B_mu_mat),
                  B_sigma_mat(_B_sigma_mat),
                  B_hmm_loglik_mat(_B_hmm_loglik_mat),
                  N_T(_N_T),
                  t_st_sig_ind_arr(_t_st_sig_ind_arr),
                  t_en_sig_ind_arr(_t_en_sig_ind_arr),
                  T_lr_loglik_arr(_T_lr_loglik_arr)
                  {}

       int M;
       vector<float>& signal_arr;         // 1 x M array
       vector<int>& state_arr;            // 1 x M array
       vector<vector<float>>& trans_mat;          // 2 x 2, 0th state is B, 1st state is T

       int N_B;
       vector<int>& b_st_sig_ind_arr;     // 1 x N_B array
       vector<int>& b_en_sig_ind_arr;     // 1 x N_B array
       vector<vector<float>>& B_mu_mat;          // N_B x N_B_MAX_MERGE
       vector<vector<float>>& B_sigma_mat;       // N_B x N_B_MAX_MERGE
       vector<vector<float>>& B_hmm_loglik_mat;  // N_B x N_B_MAX_MERGE

       int N_T;
       vector<int>& t_st_sig_ind_arr;     // 1 x N_T array
       vector<int>& t_en_sig_ind_arr;     // 1 x N_T array
       vector<float>& T_lr_loglik_arr;    // 1 x N_T array

};

/////////////////////////////////////////////////////////////////
#endif // PARAMETERS_H_INCLUDED
