#ifndef REGION_H_INCLUDED
#define REGION_H_INCLUDED
#include <vector>
#include <stack>
#include "float.h"
#include "Parameters.h"
#include "BTwindow.h"
using namespace std;

#define N_B_MAX_MERGE 5



void list_all_segmentations(int nrow, int irow, stack<int>& stk, vector<vector<int>>& res){


    if(irow >= nrow){
        vector<int> tmp_path(stk.size(), 0);
        stack<int> stk_copy(stk);
        int i = 0;
        while(!stk_copy.empty()){
            tmp_path[i] = stk_copy.top();
            stk_copy.pop();
            i += 1;
        }
        sort(tmp_path.begin(), tmp_path.end());
        res.push_back(tmp_path);
        return;

    }

    for(int i=0;i<N_B_MAX_MERGE;i++){
        if(irow + i < nrow){
            stk.push(irow * N_B_MAX_MERGE + i);
            list_all_segmentations(nrow, irow + i + 1, stk, res);
            stk.pop();
        }
    }
}



class Region
{
   public:
      Region(int _st_b_ind, int _en_b_ind, Parameters& _para):
          st_b_ind(_st_b_ind),
          en_b_ind(_en_b_ind),
          global_para(_para)
          {
              n_b = en_b_ind - st_b_ind;

              // init all_B_win_mat
              for(int i=0;i<n_b;i++){
                 vector<BTwindow> tmp_B_win_arr;
                 for(int j=0;j<N_B_MAX_MERGE;j++){
                     if(st_b_ind + i + j < en_b_ind){
                        BTwindow win = BTwindow(0, st_b_ind + i, st_b_ind + i + j + 1, global_para);
                        tmp_B_win_arr.push_back(win);
                     }else{
                        BTwindow win = BTwindow(0, -1, -1, global_para);
                        tmp_B_win_arr.push_back(win);
                     }
                }
                all_B_win_mat.push_back(tmp_B_win_arr);
              }

              // init all_B_win_mat
              for(int i=0;i<n_b-1;i++){
                BTwindow win = BTwindow(1, st_b_ind + i, st_b_ind + i + 1, global_para);
                all_T_win_arr.push_back(win);
              }

          }

      int st_b_ind;
      int en_b_ind;
      int n_b;
      Parameters& global_para;
      vector<vector<BTwindow>> all_B_win_mat;
      vector<BTwindow> all_T_win_arr;
      vector<int> best_path;
      float max_loglik = -9999;

      int get_sig_st(void);
      int get_sig_en(void);
      int get_sig_len(void);

      float cal_hmm_joint_loglik(vector<int>& path, vector<vector<float>>& trans_mat);
      void find_best_path(void);
      void set_test(float);
      void set_best_path_state(int);


};


int Region::get_sig_st(void){
    return global_para.b_st_sig_ind_arr[st_b_ind];
}


int Region::get_sig_en(void){
    return global_para.b_en_sig_ind_arr[en_b_ind - 1];
}


int Region::get_sig_len(void){
    return global_para.b_en_sig_ind_arr[en_b_ind - 1] - global_para.b_st_sig_ind_arr[st_b_ind];
}

float Region::cal_hmm_joint_loglik(vector<int>& path, vector<vector<float>>& trans_mat){

    int path_len = path.size();
    float res = (path_len - 1) * (std::log(trans_mat[0][1]) + std::log(trans_mat[1][0]));
    int p, i, j, tmp_b_sig_len, tmp_t_sig_len;

    for(p=0; p<path_len; p++){
        i = path[p] / N_B_MAX_MERGE;
        j = path[p] % N_B_MAX_MERGE;
        BTwindow tmpwin = all_B_win_mat[i][j];
        //cout << endl;
        //cout << tmpwin.st_b_ind << " -- "<<tmpwin.en_b_ind << endl;
        tmp_b_sig_len = tmpwin.get_sig_len();
        //cout << "tmp_b_sig_len : " << tmp_b_sig_len << endl;
        res = res + tmpwin.get_B_hmm_loglik() * tmp_b_sig_len + std::log(trans_mat[0][0]) * (tmp_b_sig_len - 1);
        if(p != path_len -1){  // if not the last element
            BTwindow tmp_t_win = all_T_win_arr[i + j];
            tmp_t_sig_len = tmp_t_win.get_sig_len();
            //cout << "get_T_lr_loglik : " << tmp_t_win.get_T_lr_loglik() << endl;
            res += tmp_t_win.get_T_lr_loglik() * tmp_t_sig_len + std::log(trans_mat[1][1]) * (tmp_t_sig_len - 1);
        }

    }
    return res;
}


void Region::find_best_path(void){


    vector<vector<int>> all_paths;
    stack<int> stk;
    list_all_segmentations(en_b_ind - st_b_ind, 0, stk, all_paths);

    int max_idn = 0;
    float tmp_loglik;

    for(int i = 0;i<all_paths.size();i++){

        tmp_loglik = cal_hmm_joint_loglik(all_paths[i], global_para.trans_mat);
//        cout << tmp_loglik <<  " === ";
//        for(auto p: all_paths[i]){
//            cout << p << ", ";
//        }
//        cout << endl;
        if(tmp_loglik > max_loglik){
            max_loglik = tmp_loglik;
            max_idn = i;
        }
//        cout << max_loglik <<  endl;
//        cout << endl;

    }

    for(auto p: all_paths[max_idn]){
        best_path.push_back(p);
    }
}


void Region::set_best_path_state(int center_mode){

    int n_path = best_path.size();
    assert(n_path != 0);

    vector<int> update_B_win_inds;
    vector<int> update_T_win_inds(n_path-1, -1);

    if (center_mode == 1){
       for(int i=1;i<n_path-1;i++){
         update_B_win_inds.push_back(best_path[i]);
       }
    }else if(center_mode == 0){
       for(auto p: best_path){
         update_B_win_inds.push_back(p);
       }
    }else{
        // add exception
    }

    for(int i=0;i<n_path-1;i++){
         int row = best_path[i] / N_B_MAX_MERGE;
         int coloum = best_path[i] % N_B_MAX_MERGE;
         update_T_win_inds[i] = row + coloum;
    }


    for(auto i:update_B_win_inds){
        int row = i / N_B_MAX_MERGE;
        int coloum = i % N_B_MAX_MERGE;
        all_B_win_mat[row][coloum].set_state_arr();

    }
    for(auto i:update_T_win_inds){
        all_T_win_arr[i].set_state_arr();

    }
}

#endif // REGION_H_INCLUDED
