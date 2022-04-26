#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "KargvaModel.h"
#include "../ModelDatabase.h"

using namespace std;
class KargvaDatabase : public ModelDatabase {
    public:
        KargvaDatabase();
        void SNPInfo();
};

KargvaDatabase::KargvaDatabase(){
    SNPInfo();
}
void KargvaDatabase::SNPInfo(){
    ifstream snpsearch;
    snpsearch.open("KargvaInfoExtracion/kargva_db_v5.fasta");
    string line;
    while(std::getline(snpsearch, line)) {
        if (line.at(0) != '>')
            continue;
        vector<string> header;
        header.push_back(line.substr(0, line.find("|")));
        while(line.find("|") != -1) {
            line = line.substr(line.find("|")+1);
            header.push_back(line.substr(0, line.find("|")));
        }
        if (header[4] != "NA") {
            bool included = false;
            if (snpInfoDatabase.find(header[4]) != snpInfoDatabase.end()) {
                for (auto iter = snpInfoDatabase.at(header[4]).begin(); iter != snpInfoDatabase.at(header[4]).end(); ++iter) {
                    included = included || (dynamic_cast<KargvaModel *>(*iter))->includes(header[1]);
                    if (included) {
                        (dynamic_cast<KargvaModel*>(*iter))->addToModel(header[1]);
                        break;
                    }
                }
                if (!included){
                    Model* model = new KargvaModel(header[1]);
                    snpInfoDatabase.at(header[4]).push_back(model);
                }
            }
            else {
                Model* model = new KargvaModel(header[1]);
                list<Model*> temp;
                temp.push_back(model);
                snpInfoDatabase.emplace(header[4], temp);
            }
        }
    }
    snpsearch.close();
}