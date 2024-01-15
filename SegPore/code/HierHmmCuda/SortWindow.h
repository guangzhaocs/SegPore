#ifndef SORTWINDOW_H_INCLUDED
#define SORTWINDOW_H_INCLUDED
#include <vector>
using namespace std;

#define N_B_MAX_MERGE 5

struct SortWindow{   // rename: WinIndex


    int b_st_sig_ind;
    int b_en_sig_ind;
    int sig_len;

    float B_mu;
    float B_sigma;

    float B_prev_mu;
    float B_prev_sigma;

    float B_next_mu;
    float B_next_sigma;

	int old_window_idx;
};


bool LessSort(SortWindow a,SortWindow b){
    return (a.sig_len<b.sig_len);
}

void generateWinRegion(int N_B, vector<int>& b_st_sig_ind_arr, vector<int>& b_en_sig_ind_arr,
                       vector<vector<float>>& B_mu_mat, vector<vector<float>>& B_sigma_mat,
                       vector<vector<float>>& B_prev_mu_mat, vector<vector<float>>& B_prev_sigma_mat,
                       vector<vector<float>>& B_next_mu_mat, vector<vector<float>>& B_next_sigma_mat,
                       vector<SortWindow>& sort_win_arr){

    for(int i=0;i<N_B;i++){
        for(int j=0;j<N_B_MAX_MERGE;j++){
            if(i + j >= N_B){
                continue;
            }
            SortWindow temp;
            temp.b_st_sig_ind = b_st_sig_ind_arr[i];
            temp.b_en_sig_ind = b_en_sig_ind_arr[i + j];
            temp.sig_len = b_en_sig_ind_arr[i + j] - b_st_sig_ind_arr[i];

            temp.B_mu = B_mu_mat[i][j];
            temp.B_sigma = B_sigma_mat[i][j];

            temp.B_prev_mu = B_prev_mu_mat[i][j];
            temp.B_prev_sigma = B_prev_sigma_mat[i][j];

            temp.B_next_mu = B_next_mu_mat[i][j];
            temp.B_next_sigma = B_next_sigma_mat[i][j];

            temp.old_window_idx = i * N_B_MAX_MERGE + j;

            sort_win_arr.push_back(temp);
        }
    }


    // sort
    sort(sort_win_arr.begin(), sort_win_arr.end(), LessSort);

}

#endif // SORTWINDOW_H_INCLUDED
