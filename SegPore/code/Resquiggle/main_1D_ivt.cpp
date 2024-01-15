#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <cassert>
#include <algorithm>
#include <math.h>
#include <float.h>
using namespace std;


void read2DFloatFile(string fileName, vector<vector<float>>& _mat){
    ifstream inFile(fileName);
    string line;
    while (getline(inFile, line))
    {
        stringstream ss(line);
        string str;
        vector<float> _row;
        while (getline(ss, str, ',')){
            float str_float;
            istringstream istr(str);
            istr >> str_float;
            _row.push_back(str_float);
        }
        _mat.push_back(_row);
    }
}


void read2DIntFile(string fileName, vector<vector<int>>& _mat){
    ifstream inFile(fileName);
    string line;
    while (getline(inFile, line))
    {
        stringstream ss(line);
        string str;
        vector<int> _row;
        while (getline(ss, str, ',')){
            int str_float;
            istringstream istr(str);
            istr >> str_float;
            _row.push_back(str_float);
        }
        _mat.push_back(_row);
    }
}


void read2DStrFile(string fileName, vector<vector<string>>& _mat){
    ifstream inFile(fileName);
    string line;
    while (getline(inFile, line))
    {
        stringstream ss(line);
        string str;
        vector<string> _row;
        while (getline(ss, str, ',')){
            _row.push_back(str);
        }
        _mat.push_back(_row);
    }
}


void readKmerFile(string fileName, map<string, vector<float> >& ref_kmer_dict){
    ifstream inFile(fileName);
    string line;
    while (getline(inFile, line))
    {
        stringstream ss(line);
        string kmer;
        while (getline(ss, kmer, ',')){
            break;
        }

        string str;
        vector<float> _row;
        while (getline(ss, str, ',')){
            float str_float;
            istringstream istr(str);
            istr >> str_float;
            _row.push_back(str_float);
        }
        ref_kmer_dict.insert({kmer, _row});
    }
}


void label2seq(vector<int>& label_arr, int contig_pos_start, vector<string>& kmer_str_arr, vector<int>& kmer_pos_arr){

    kmer_str_arr.push_back("-");
    kmer_pos_arr.push_back(-1);

    map<int, string> Id2Base;
    Id2Base[1] = 'A';
    Id2Base[2] = 'C';
    Id2Base[3] = 'G';
    Id2Base[4] = 'T';
    Id2Base[5] = 'N';

    int _pos = contig_pos_start;
    string _kmer;
    string null_c = "N";
    for(int i = 2;i<static_cast<int>(label_arr.size())-2; i++){
        _kmer = Id2Base[label_arr[i-2]] + Id2Base[label_arr[i-1]] + Id2Base[label_arr[i]] +
                Id2Base[label_arr[i+1]] + Id2Base[label_arr[i+2]];

         string::size_type idx;
         idx=_kmer.find(null_c);
         if(idx == string::npos){
             kmer_str_arr.push_back(_kmer);
             kmer_pos_arr.push_back(_pos);
         }
         _pos += 1;
    }
}


float log_score(map<string, vector<float> >& ref_kmer_dict, float mu_val, string kmer){

    static const float inv_sqrt_2pi = 0.3989422804014327;
    static const float unif_log_prob = -6.0;

    if(std::abs(mu_val) < 1e-10 || kmer == "-"){
        return unif_log_prob;
    }else{
        float kmer_mu1 = ref_kmer_dict[kmer][0];
        float kmer_sigma1 = ref_kmer_dict[kmer][1];
        float a = (mu_val - kmer_mu1) / kmer_sigma1;
        return std::log(inv_sqrt_2pi / kmer_sigma1 * std::exp(-0.5f * a * a));
    }
}


float getMax(float v1, float v2, float v3, float v4){
    float v_max = v1;
    if(v2 > v_max){
       v_max = v2;
    }
    if(v3 > v_max){
       v_max = v3;
    }
    if(v4 > v_max){
       v_max = v4;
    }
    return v_max;
}


int getFlag(float v_max, float v1, float v2, float v3, float v4){

    if (v_max == v1 && v_max != v2 && v_max != v3 && v_max != v4){
        return 1;
    }else if(v_max != v1 && v_max == v2 && v_max != v3 && v_max != v4){
        return 2;
    }else if(v_max != v1 && v_max != v2 && v_max == v3 && v_max != v4){
        return 3;
    }else if(v_max != v1 && v_max != v2 && v_max != v3 && v_max == v4){
        return 4;
    }else if(v_max == v1 && v_max == v2 && v_max != v3 && v_max != v4){
        return 5;
    }else if(v_max == v1 && v_max != v2 && v_max == v3 && v_max != v4){
        return 6;
    }else if(v_max == v1 && v_max != v2 && v_max != v3 && v_max == v4){
        return 7;
    }else if(v_max != v1 && v_max == v2 && v_max == v3 && v_max != v4){
        return 8;
    }else if(v_max != v1 && v_max == v2 && v_max != v3 && v_max == v4){
        return 9;
    }else if(v_max != v1 && v_max != v2 && v_max == v3 && v_max == v4){
        return 10;
    }else if(v_max == v1 && v_max == v2 && v_max == v3 && v_max != v4){
        return 11;
    }else if(v_max == v1 && v_max == v2 && v_max != v3 && v_max == v4){
        return 12;
    }else if(v_max == v1 && v_max != v2 && v_max == v3 && v_max == v4){
        return 13;
    }else if(v_max != v1 && v_max == v2 && v_max == v3 && v_max == v4){
        return 14;
    }else{
        return 15;
    }
}


void run_one_read(vector<float>& mu_arr, vector<float>& sigma_arr, vector<int>& label_arr,  vector<int>& start_pos_arr,
                  vector<int>& end_pos_arr, vector<int>& len_arr, vector<string>& read_info_arr,
                  map<string, vector<float> >& ref_kmer_dict, string save_file_dir){

    // generate full mu arr
    int n_mu = static_cast<int>(mu_arr.size());
    assert (n_mu == static_cast<int>(sigma_arr.size()));
    assert (n_mu == static_cast<int>(start_pos_arr.size()));
    assert (n_mu == static_cast<int>(end_pos_arr.size()));
    assert (n_mu == static_cast<int>(len_arr.size()));

    vector<float> mu_val_arr(n_mu + 1, 0.0);
    vector<float> sigma_val_arr(n_mu + 1, 0.0);
    vector<int> start_pos_val_arr(n_mu + 1, -1);
    vector<int> end_pos_val_arr(n_mu + 1, -1);
    vector<int> len_pos_val_arr(n_mu + 1, -1);

    int k = 1;
    for(int i=n_mu-1;i>=0;i--){  // reverse the current signal
        mu_val_arr[k] = mu_arr[i];
        sigma_val_arr[k] = sigma_arr[i];
        start_pos_val_arr[k] = start_pos_arr[i];
        end_pos_val_arr[k] = end_pos_arr[i];
        len_pos_val_arr[k] = len_arr[i];
        k ++;
    }

    // split read info
    string read_name = read_info_arr[0];
    string read_index = read_info_arr[1];
    string contig = read_info_arr[2];
    int contig_pos_start = 0;

    // generate full kmer arr
    vector<string> kmer_str_arr;
    vector<int> kmer_pos_arr;
    label2seq(label_arr, contig_pos_start, kmer_str_arr, kmer_pos_arr);
    assert(kmer_str_arr.size() == kmer_pos_arr.size());

    // dp matrix and trace back matrix
    int m = mu_val_arr.size();
    int n = kmer_str_arr.size();
    vector<vector<float> > dp_mat(m, vector<float> (n, 0.0));
    vector<vector<int> > trace_mat(m, vector<int> (n, 0));
    cout << " * dp mat size : " << m << " x " << n << endl;

    // first row
    for(int i=1; i<n; i++){
        dp_mat[0][i] = dp_mat[0][i - 1] + log_score(ref_kmer_dict, mu_val_arr[0], kmer_str_arr[i]);
        trace_mat[0][i] = 2;
    }

    // first column
    for(int i=1; i<m; i++){
        dp_mat[i][0] = dp_mat[i - 1][0] + log_score(ref_kmer_dict, mu_val_arr[i], kmer_str_arr[0]);
        trace_mat[i][0] = 3;
    }

    // v1:  M[i, j] = M[i-1, j-1] + s(mu_i, kmer_j)
    // v2:  M[i, j] = M[i, j-1] + s(0.0, kmer_j)
    // v3:  M[i, j] = M[i-1, j] + s(mu_i, "-")
    // v4:  M[i, j] = M[i-1, j] + s(mu_i, kmer_j)  repeat the last kmer
    // float v_max=0, v1, v2, v3, v4;
    float v1, v4;
    for(int i=1; i<m; i++){
        for(int j=1; j<n; j++){
            v1 = dp_mat[i - 1][j - 1] + log_score(ref_kmer_dict, mu_val_arr[i], kmer_str_arr[j]);
            //v2 = dp_mat[i][j - 1] + log_score(ref_kmer_dict, mu_val_arr[0], kmer_str_arr[j]);
            //v3 = dp_mat[i - 1][j] + log_score(ref_kmer_dict, mu_val_arr[i], kmer_str_arr[0]);
            v4 = dp_mat[i - 1][j] + log_score(ref_kmer_dict, mu_val_arr[i], kmer_str_arr[j]);

            //v_max = getMax(v1, v2, v3, v4);
            //dp_mat[i][j] = v_max;
            //trace_mat[i][j] = getFlag(v_max, v1, v2, v3, v4);

            if(v1 > v4){
                dp_mat[i][j] = v1;
                trace_mat[i][j] = 1;
            }else if (v1 < v4){
                dp_mat[i][j] = v4;
                trace_mat[i][j] = 4;
            }else{
                dp_mat[i][j] = v1;
                trace_mat[i][j] = 7;
            }
        }
    }

    cout << " * dp end ! " << endl;

//    for(int i=0; i<10; i++){
//        for(int j=0; j<10; j++){
//            cout << dp_mat[i][j] << " ";
//        }
//        cout << endl;
//    }
//
//    for(int i=0; i<10; i++){
//        for(int j=0; j<10; j++){
//            cout << trace_mat[i][j] << " ";
//        }
//        cout << endl;
//    }

    // trace back
    vector<float> pred_mu_val_arr;
    vector<float> pred_sigma_val_arr;
    vector<int> pred_start_pos_val_arr;
    vector<int> pred_end_pos_val_arr;
    vector<int> pred_len_pos_val_arr;
    vector<string> pred_kmer_str_arr;
    vector<int> pred_kmer_pos_arr;


    //float max_dp_start_v = - DBL_MAX;
    //int max_start_idx = 0;
    //for(int d = 0; d < n; d++){
    //    if(dp_mat[1][d] > max_dp_start_v){
    //        max_dp_start_v = dp_mat[1][d];
    //        max_start_idx = d;
    //    }
    //}
    //cout << " * max start value is in " << max_start_idx << endl;


    float max_dp_v = - DBL_MAX;
    int max_idx = 400;
    for(int d = 400; d < n; d++){
        if(dp_mat[m-1][d] > max_dp_v){
            max_dp_v = dp_mat[m-1][d];
            max_idx = d;
        }
    }
    cout << " * max value is in " << max_idx << endl;


    int i = m - 1;
    //int j = n - 1;
    int j = max_idx;
    while(i > 0 && j > 0){
        if(trace_mat[i][j] == 2){
            pred_mu_val_arr.push_back(mu_val_arr[0]);
            pred_sigma_val_arr.push_back(sigma_val_arr[0]);
            pred_start_pos_val_arr.push_back(start_pos_val_arr[0]);
            pred_end_pos_val_arr.push_back(end_pos_val_arr[0]);
            pred_len_pos_val_arr.push_back(len_pos_val_arr[0]);
            pred_kmer_str_arr.push_back(kmer_str_arr[j]);
            pred_kmer_pos_arr.push_back(kmer_pos_arr[j]);
            j = j - 1;
        }else if(trace_mat[i][j] == 3 || trace_mat[i][j] == 8){
            pred_mu_val_arr.push_back(mu_val_arr[i]);
            pred_sigma_val_arr.push_back(sigma_val_arr[i]);
            pred_start_pos_val_arr.push_back(start_pos_val_arr[i]);
            pred_end_pos_val_arr.push_back(end_pos_val_arr[i]);
            pred_len_pos_val_arr.push_back(len_pos_val_arr[i]);
            pred_kmer_str_arr.push_back(kmer_str_arr[0]);
            pred_kmer_pos_arr.push_back(kmer_pos_arr[0]);
            i = i - 1;
        }else if(trace_mat[i][j] == 4 || trace_mat[i][j] == 9 || trace_mat[i][j] == 10 || trace_mat[i][j] == 14){
            pred_mu_val_arr.push_back(mu_val_arr[i]);
            pred_sigma_val_arr.push_back(sigma_val_arr[i]);
            pred_start_pos_val_arr.push_back(start_pos_val_arr[i]);
            pred_end_pos_val_arr.push_back(end_pos_val_arr[i]);
            pred_len_pos_val_arr.push_back(len_pos_val_arr[i]);
            pred_kmer_str_arr.push_back(kmer_str_arr[j]);
            pred_kmer_pos_arr.push_back(kmer_pos_arr[j]);
            i = i - 1;
        }else{ // trace_mat[i][j] in [1, 5, 6, 7, 11, 12, 13, 15]
            pred_mu_val_arr.push_back(mu_val_arr[i]);
            pred_sigma_val_arr.push_back(sigma_val_arr[i]);
            pred_start_pos_val_arr.push_back(start_pos_val_arr[i]);
            pred_end_pos_val_arr.push_back(end_pos_val_arr[i]);
            pred_len_pos_val_arr.push_back(len_pos_val_arr[i]);
            pred_kmer_str_arr.push_back(kmer_str_arr[j]);
            pred_kmer_pos_arr.push_back(kmer_pos_arr[j]);
            i = i - 1;
            j = j - 1;
        }
    }

    if(i == 0 && j != 0){
        for(int idx=j;idx>0;idx--){
            pred_mu_val_arr.push_back(mu_val_arr[0]);
            pred_sigma_val_arr.push_back(sigma_val_arr[0]);
            pred_start_pos_val_arr.push_back(start_pos_val_arr[0]);
            pred_end_pos_val_arr.push_back(end_pos_val_arr[0]);
            pred_len_pos_val_arr.push_back(len_pos_val_arr[0]);
            pred_kmer_str_arr.push_back(kmer_str_arr[idx]);
            pred_kmer_pos_arr.push_back(kmer_pos_arr[idx]);
        }
    }else if(i != 0 && j == 0){
        for(int idx=i;idx>0;idx--){
            pred_mu_val_arr.push_back(mu_val_arr[idx]);
            pred_sigma_val_arr.push_back(sigma_val_arr[idx]);
            pred_start_pos_val_arr.push_back(start_pos_val_arr[idx]);
            pred_end_pos_val_arr.push_back(end_pos_val_arr[i]);
            pred_len_pos_val_arr.push_back(len_pos_val_arr[i]);
            pred_kmer_str_arr.push_back(kmer_str_arr[0]);
            pred_kmer_pos_arr.push_back(kmer_pos_arr[0]);
        }
    }else{}

    reverse(pred_mu_val_arr.begin(), pred_mu_val_arr.end());
    reverse(pred_sigma_val_arr.begin(), pred_sigma_val_arr.end());
    reverse(pred_start_pos_val_arr.begin(), pred_start_pos_val_arr.end());
    reverse(pred_end_pos_val_arr.begin(), pred_end_pos_val_arr.end());
    reverse(pred_len_pos_val_arr.begin(), pred_len_pos_val_arr.end());
    reverse(pred_kmer_str_arr.begin(), pred_kmer_str_arr.end());
    reverse(pred_kmer_pos_arr.begin(), pred_kmer_pos_arr.end());

    int n_pred  = static_cast<int>(pred_mu_val_arr.size());
    assert(n_pred == static_cast<int>(pred_sigma_val_arr.size()));
    assert(n_pred == static_cast<int>(pred_start_pos_val_arr.size()));
    assert(n_pred == static_cast<int>(pred_end_pos_val_arr.size()));
    assert(n_pred == static_cast<int>(pred_kmer_str_arr.size()));
    assert(n_pred == static_cast<int>(pred_kmer_pos_arr.size()));
    assert(n_pred == static_cast<int>(pred_len_pos_val_arr.size()));

    ofstream myFile;
    myFile.open(save_file_dir, ios::app);
    for(int i=0;i<static_cast<int>(pred_mu_val_arr.size());i++){
        myFile << read_index << "\t";
        myFile << read_name << "\t";
        myFile << contig << "\t";
        myFile << pred_kmer_str_arr[i] << "\t";
        myFile << pred_kmer_pos_arr[i] << "\t";
        myFile << pred_mu_val_arr[i] << "\t";
        myFile << pred_sigma_val_arr[i] << "\t";
        myFile << pred_start_pos_val_arr[i] << "\t";
        myFile << pred_end_pos_val_arr[i] << "\t";
        myFile << pred_len_pos_val_arr[i] << "\t";
        myFile << "\n";
    }
    myFile.close();
}


int main(int argc, char ** argv)
{
    // file name
    string root_dir = argv[1];
    cout << root_dir << endl;

	string eventalign_fileName = argv[2];
    cout << eventalign_fileName << endl;

    string ref_kmer_file = argv[3];
    cout << ref_kmer_file << endl;


    string mu_dir = root_dir + "/mu.csv";
    string sigma_dir = root_dir + "/sigma.csv";
    string len_dir = root_dir + "/len.csv";
    string label_dir = root_dir + "/label.csv";
    string start_pos_dir = root_dir + "/start_pos.csv";
    string end_pos_dir = root_dir + "/end_pos.csv";
    string read_info_dir = root_dir + "/readname.csv";
    string save_file_dir = root_dir + "/" + eventalign_fileName;
    string ref_kmer_dir = "/scratch/cs/nanopore/chengg1/base_calling/dataset/xpore/0_reference/model_kmer_motif/" + ref_kmer_file;

    // read signal files and border files
    vector<vector<float> > all_mu_arr;
    vector<vector<float> > all_sigma_arr;
    vector<vector<int> > all_label_arr;
    vector<vector<int> > all_start_pos_arr;
    vector<vector<int> > all_end_pos_arr;
    vector<vector<int> > all_len_arr;
    vector<vector<string> > all_read_info_arr;
    map <string, vector<float>> ref_kmer_dict;

    read2DFloatFile(mu_dir, all_mu_arr);
    read2DFloatFile(sigma_dir, all_sigma_arr);
    read2DIntFile(label_dir, all_label_arr);
    read2DIntFile(start_pos_dir, all_start_pos_arr);
    read2DIntFile(end_pos_dir, all_end_pos_arr);
    read2DIntFile(len_dir, all_len_arr);
    read2DStrFile(read_info_dir, all_read_info_arr);
    readKmerFile(ref_kmer_dir, ref_kmer_dict);

    int n_read  = static_cast<int>(all_mu_arr.size());
    assert(n_read == static_cast<int>(all_sigma_arr.size()));
    assert(n_read == static_cast<int>(all_label_arr.size()));
    assert(n_read == static_cast<int>(all_end_pos_arr.size()));
    assert(n_read == static_cast<int>(all_start_pos_arr.size()));
    assert(n_read == static_cast<int>(all_len_arr.size()));
    assert(n_read == static_cast<int>(all_read_info_arr.size()));

    for(int idx=0; idx<n_read; idx++){
        cout << " --- run idx : " << idx << " --- " << endl;
        run_one_read(all_mu_arr[idx], all_sigma_arr[idx], all_label_arr[idx],  all_start_pos_arr[idx],
                     all_end_pos_arr[idx], all_len_arr[idx], all_read_info_arr[idx], ref_kmer_dict, save_file_dir);

    }
    return 0;
}
