#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "KargvaModelDeletion.h"
#include "KargvaModelInsertion.h"
#include "KargvaModelNonsense.h"
#include "KargvaModelReg.h"
#include "KargvaMultipleModels.h"
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
                if (header[1].find(";") == -1 && header[1].find("STOP") != -1 && !(header[1].at(0) == '-' && isalpha(header[1].at(1)))) {
                    for (auto iter = snpInfoDatabase.at(header[4]).begin(); iter != snpInfoDatabase.at(header[4]).end(); ++iter) {
                        included = included || (dynamic_cast<KargvaModel*>(*iter))->includes(header[1]);
                        if (included) {
                            (dynamic_cast<KargvaModel*>(*iter))->addToModel(header[1]);
                            break;
                        }
                    }
                }
                if (!included){
                    Model* model;
                    if (header[1].at(0) == '-' && isalpha(header[1].at(1)))
                        model = new KargvaModelDeletion(header[1], header[2]);
                    else if (header[1].at(0) == '-')
                        model = new KargvaModelInsertion(header[1], header[2]);
                    else if (header[1].find("STOP") != -1)
                        model = new KargvaModelNonsense(header[1], header[2]);
                    else if (header[1].find(";") == -1)
                        model = new KargvaMultipleModels(header[1], header[2]);
                    else
                        model = new KargvaModelReg(header[1], header[2]);
                    snpInfoDatabase.at(header[4]).push_back(model);
                }
            }
            else {
                Model* model;
                if (header[1].at(0) == '-' && isalpha(header[1].at(1)))
                    model = new KargvaModelDeletion(header[1], header[2]);
                else if (header[1].at(0) == '-')
                    model = new KargvaModelInsertion(header[1], header[2]);
                else if (header[1].find("STOP") != -1)
                    model = new KargvaModelNonsense(header[1], header[2]);
                else if (header[1].find(";") == -1)
                    model = new KargvaMultipleModels(header[1], header[2]);
                else
                    model = new KargvaModelReg(header[1], header[2]);
                list<Model*> temp;
                temp.push_back(model);
                snpInfoDatabase.emplace(header[4], temp);
            }
        }
    }
    snpsearch.close();
}