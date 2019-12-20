#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<Rcpp.h>


//[[Rcpp::plugins(cpp14)]]
//[[Rcpp::export]]

bool sequenceToFastaConts(std::vector<int>& starts,std::vector<int>& widths,std::string& sequence,std::string& newFasta,std::string& nameTag){
    bool x = true;
    if(std::ifstream(newFasta)){
        x = false;
    }
    std::ofstream outfile (x ? std::ofstream(newFasta):std::ofstream(newFasta,std::ios::app));
    if(outfile.is_open()){
        for(int i = 0;i < (int) starts.size();i++){
            outfile << ">"+nameTag << std::endl;
            outfile << sequence.substr(starts[i],widths[i]) << std::endl;
        }
        outfile.close();
        return true;
    }
    else{
        return false;
    }
}
