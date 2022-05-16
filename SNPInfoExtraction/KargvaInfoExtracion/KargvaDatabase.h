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
#include "MEG_6090MultipleModels.h"
#include "../ModelDatabase.h"

using namespace std;
class KargvaDatabase : public ModelDatabase {
    private:
        shared_ptr<CARD_database> databaseSequences;
        void addSNP(string snp, string id, string meg);
    public:
        KargvaDatabase();
        void SNPInfo();
};

KargvaDatabase::KargvaDatabase(){
    databaseSequences = make_shared<CARD_database>();
    SNPInfo();
}
void KargvaDatabase::addSNP(string snp, string id, string meg) {
    bool included = false;
    if (snpInfoDatabase.find(meg) != snpInfoDatabase.end()) {
        if (snp.find("STOP") == -1) {
            for (auto iter = snpInfoDatabase.at(meg).begin(); iter != snpInfoDatabase.at(meg).end(); ++iter) {
                KargvaModelReg* tempReg = dynamic_cast<KargvaModelReg*>(*iter);
                if (tempReg != nullptr) {
                    included = tempReg->includes(snp);
                    if (included) {
                        if (snp.find(";") == -1 && snp.find(",") == -1)
                            tempReg->addToModel(snp);
                        break;
                    }
                }
            }
        }
        if (!included) {
            Model* model;
            if (snp.find(";") != -1)
                model = new KargvaMultipleModels(snp, id, databaseSequences);
            else if (snp.find(",") != -1)
                model = new MEG_6090MultipleModels(snp, id, databaseSequences);
            else if (snp.find("-") != -1)
                throw;
            else if (snp.find("STOP") != -1)
                model = new KargvaModelNonsense(snp, id, databaseSequences);
            else
                model = new KargvaModelReg(snp, id, databaseSequences);
            snpInfoDatabase.at(meg).push_back(model->Clone());
            delete model;
        }
    }
    else {
        Model* model;
        if (snp.find(";") != -1 || snp.find(",") != -1 || snp.find("-") != -1)
            throw;
        else if (snp.find("STOP") != -1)
            model = new KargvaModelNonsense(snp, id, databaseSequences);
        else
            model = new KargvaModelReg(snp, id, databaseSequences);
        list<Model*> temp;
        temp.push_back(model->Clone());
        snpInfoDatabase.emplace(meg, temp);
        delete model;
    }
}
void KargvaDatabase::SNPInfo(){
    vector<vector<string>> mult;
    ifstream snpsearch;
    snpsearch.open("KargvaInfoExtracion/kargva_db_v5.fasta");
    string line;
    bool Q142X = false;
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
        else if (header[2] == "ARO:3003686") {
            if (Q142X)
                continue;
            else
                Q142X = true;
            if (header[1] == "Q142X")
                header[1] = "Q142STOP";
        }
        else if (header[2] == "ARO:3003042") {
            if (header[1] == "T488A")
                header[1] = "T494A";
            else if (header[1] == "E475G")
                header[1] = "E481G";
            else if (header[1] == "T445A")
                header[1] = "T451A";
        }
        else if (header[2] == "ARO:3003283") {
            continue;
        }
        if (header[4] != "NA" && header[4] != "not corresponding") {
            if (header[1].find(";") != -1 || header[1].find(",") != -1) {
                vector<string> temp = { header[1], header[2], header[4] };
                mult.push_back(temp);
            }
            else
                this->addSNP(header[1], header[2], header[4]);
        }
    }
    snpsearch.close();
    snpsearch.open("KargvaInfoExtracion/MEG_6090_SNP_list_update2.txt");
    while (std::getline(snpsearch, line)) {
        if (line.find(";") != -1 || line.find(",") != -1) {
            vector<string> temp = { line, "ARO:3003283", "MEG_6090" };
            mult.push_back(temp);
        }
        else
            this->addSNP(line, "ARO:3003283", "MEG_6090");
    }
    snpsearch.close();
    for (auto iter = mult.begin(); iter != mult.end(); ++iter) {
        this->addSNP((* iter)[0], (* iter)[1], (* iter)[2]);
    }
}