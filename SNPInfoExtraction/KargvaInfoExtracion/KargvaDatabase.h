#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "KargvaModelDeletion.h"
#include "KargvaModelNonsense.h"
#include "KargvaModelReg.h"
#include "KargvaMultipleModels.h"
#include "../ModelDatabase.h"

using namespace std;
class KargvaDatabase : public ModelDatabase {
    private:
        CARD_database* databaseSequences;
    public:
        KargvaDatabase();
        void SNPInfo();
};

KargvaDatabase::KargvaDatabase(){
    databaseSequences = new CARD_database();
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
        if (header[2] == "ARO:3003326") {
            if (header[1] == "A314G;Y322C")
                header[1] = "A313G;Y319C";
        }
        else if (header[2] == "ARO:3003327") {
            if (header[1] == "S244T")
                continue;
        }
        else if (header[2] == "ARO:3003735") {
            if (header[1] == "F441Y" || header[1] == "T387I;E449K")
                continue;
        }
        else if (header[2] == "ARO:3004562") {
            if (header[1] == "E467V")
                header[1] = "E466V";
            else if (header[1] == "E467K")
                header[1] = "E466K";
            else if (header[1] == "K444F")
                header[1] = "L444F";
            else if (header[1] == "R388P")
                header[1] = "R389P";
        }
        else if (header[2] == "ARO:3003459") {
            if (header[1] == "N510D" || header[1] == "E501D" || header[1] == "N499T")
                continue;
        }
        else if (header[2] == "ARO:3003302") {
            if (header[1] == "G125S")
                header[1] = "G124S";
        }
        else if (header[2] == "ARO:3003392") {
            if (header[1] == "A234G" || header[1] == "A431V")
                continue;
        }
        else if (header[2] == "ARO:3003319") {
            if (header[1] == "I418N")
                continue;
        }
        else if (header[2] == "ARO:3003702") {
            if (header[1] == "S80L")
                header[1] = "S87L";
        }
        else if (header[2] == "ARO:3003394") {
            if (header[1] == "A8G")
                header[1] = "D8G";
            else if (header[1] == "Y68D")
                header[1] = "W68D";
            else if (header[1] == "A140S")
                header[1] = "R140S";
            else if (header[1] == "L72P")
                header[1] = "L172P";
            else if (header[1] == "R157W")
                header[1] = "V157W";
        }
        else if (header[2] == "ARO:3004153") {
            if (header[1] == "I211V")
                continue;
        }
        else if (header[2] == "ARO:3003970") {
            if (header[1] == "G98R")
                header[1] = "G99R";
        }
        else if (header[2] == "ARO:3003283") {
            //continue;
        }
        if (header[4] != "NA") {
            bool included = false;
            if (snpInfoDatabase.find(header[4]) != snpInfoDatabase.end()) {
                if (header[1].find(";") == -1 || header[1].find("STOP") == -1 || header[1].at(0) != '-') {
                    for (auto iter = snpInfoDatabase.at(header[4]).begin(); iter != snpInfoDatabase.at(header[4]).end(); ++iter) {
                        KargvaModelReg* tempReg = dynamic_cast<KargvaModelReg*>(*iter);
                        if (tempReg != nullptr) {
                            included = tempReg->includes(header[1]);
                            if (included) {
                                tempReg->addToModel(header[1]);
                                break;
                            }
                        }
                    }
                }
                if (!included){
                    Model* model;
                    if (header[1].at(0) == '-')
                        model = new KargvaModelDeletion(header[1], header[2], databaseSequences);
                    else if (header[1].find("STOP") != -1)
                        model = new KargvaModelNonsense(header[1], header[2], databaseSequences);
                    else if (header[1].find(";") != -1)
                        model = new KargvaMultipleModels(header[1], header[2], databaseSequences);
                    else
                        model = new KargvaModelReg(header[1], header[2], databaseSequences);
                    snpInfoDatabase.at(header[4]).push_back(model);
                }
            }
            else {
                Model* model;
                if (header[1].at(0) == '-')
                    model = new KargvaModelDeletion(header[1], header[2], databaseSequences);
                else if (header[1].find("STOP") != -1)
                    model = new KargvaModelNonsense(header[1], header[2], databaseSequences);
                else if (header[1].find(";") != -1)
                    model = new KargvaMultipleModels(header[1], header[2], databaseSequences);
                else
                    model = new KargvaModelReg(header[1], header[2], databaseSequences);
                list<Model*> temp;
                temp.push_back(model);
                snpInfoDatabase.emplace(header[4], temp);
            }
        }
    }
    snpsearch.close();
}