#pragma once
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include "KargvaAaDeletion.h"
#include "KargvaAaNonsense.h"
#include "KargvaAaMissense.h"
#include "KargvaAaMultiple.h"
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
        if (snp.find("FS-") != -1) throw;
        if (snp.find("STOP") == -1) {
            for (auto iter = snpInfoDatabase.at(meg).begin(); iter != snpInfoDatabase.at(meg).end(); ++iter) {
                KargvaAaMissense* tempReg = dynamic_cast<KargvaAaMissense*>(*iter);
                if (tempReg != nullptr) {
                    included = tempReg->includes(snp);
                    if (included) {
                        if (snp.find(";") == -1)
                            tempReg->addToModel(snp);
                        break;
                    }
                }
            }
        }
        if (!included) {
            Model* model;
            if (snp.find(";") != -1)
                model = new KargvaAaMultiple(snp, id, databaseSequences);
            else if (snp.find("-") != -1)
                model = new KargvaAaDeletion(snp, id, databaseSequences);
            else if (snp.find("STOP") != -1)
                model = new KargvaAaNonsense(snp, id, databaseSequences);
            else
                model = new KargvaAaMissense(snp, id, databaseSequences);
            snpInfoDatabase.at(meg).push_back(model->Clone());
            delete model;
        }
    }
    else {
        Model* model;
        if (snp.find(";") != -1 || snp.find("-") != -1)
            throw;
        else if (snp.find("STOP") != -1)
            model = new KargvaAaNonsense(snp, id, databaseSequences);
        else
            model = new KargvaAaMissense(snp, id, databaseSequences);
        list<InfoPipe*> temp;
        temp.push_back(model->Clone());
        snpInfoDatabase.emplace(meg, temp);
        delete model;
    }
}
void KargvaDatabase::SNPInfo(){
    vector<string> corrected = { "ARO:3003374", "ARO:3003373", "ARO:3004145", "ARO:3003760", "ARO:3003092", "ARO:3003761", "ARO:3003326",
                                "ARO:3003327", "ARO:3003458", "ARO:3003470", "ARO:3003901", "ARO:3003889", "ARO:3003940", "ARO:3003295",
                                "ARO:3003298", "ARO:3004631", "ARO:3003294", "ARO:3003305", "ARO:3004562", "ARO:3003459", "ARO:3003303",
                                "ARO:3003302", "ARO:3003393", "ARO:3003448", "ARO:3003392", "ARO:3003077", "ARO:3003790", "ARO:3003028",
                                "ARO:3003028", "ARO:3003574", "ARO:3003777", "ARO:3003785", "ARO:3003784", "ARO:3003461", "ARO:3003779",
                                "ARO:3003686", "ARO:3003309", "ARO:3004630", "ARO:3003941", "ARO:3003316", "ARO:3003937", "ARO:3003323",
                                "ARO:3003896", "ARO:3003394", "ARO:3003380", "ARO:3003379", "ARO:3003283", "ARO:3003288", "ARO:3004563",
                                "ARO:3004721", "ARO:3003479", "ARO:3004153", "ARO:3003445", "ARO:3003465", "ARO:3003359", "ARO:3003902",
                                "ARO:3003970", "ARO:3003325", "ARO:3003387", "ARO:3003895", "ARO:3003438"};
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
        bool needsCorrection = false;
        for (auto iter = corrected.begin(); iter != corrected.end(); ++iter) {
            if (*iter == header[2]) {
                needsCorrection = true;
                break;
            }
        }
        if (needsCorrection) {
            continue;
        }
        else if (header[2] == "ARO:3003735") {          //MEG_3065
            if (header[1] == "F441Y" || header[1] == "T387I;E449K")
                continue;
        }
        else if (header[2] == "ARO:3003928") {          //MEG_3177
            if (header[1] == "D95F")
                continue;
        }
        else if (header[2] == "ARO:3003298") {          //MEG_3181
            if (header[1] == "A91T")
                header[1] = "A91V";
        }
        else if (header[2] == "ARO:3003306") {          //MEG_3242
            if (header[1] == "S464A")
                header[1] = "S463A;S464Y";
            else
                continue;
        }
        else if (header[2] == "ARO:3003301") {          //MEG_3245
            if (header[1] == "S128L" || header[1] == "I55S")
                continue;
        }
        else if (header[2] == "ARO:3003463") {          //MEG_3445
            if (header[1] == "G312S")
                continue;
        }
        else if (header[2] == "ARO:3003079") {          //MEG_3589
            if (header[1] == "A180T" || header[1] == "H264Q")
                continue;
        }
        else if (header[2] == "ARO:3003319") {          //MEG_4057
            if (header[1] == "I418N")
                header[1] = "I422N";
        }
        else if (header[2] == "ARO:3003776") {          //MEG_4094
            if (header[1] != "L42*")
                continue;
        }
        else if (header[2] == "ARO:3003775") {          //MEG_4095
            if (header[1] == "C115S")
                continue;
        }
        else if (header[2] == "ARO:3003390") {          //MEG_4289
            if (header[1] == "R154A" || header[1] == "R154D")
                continue;
        }
        else if (header[2] == "ARO:3003702") {          //MEG_5325
            if (header[1] == "S80L")
                header[1] = "S87L";
        }
        else if (header[2] == "ARO:3003308") {          //MEG_5333
            if (header[1] == "E84G")
                continue;
        }
        else if (header[2] == "ARO:3003042") {          //MEG_5406
            if (header[1] == "T488A")
                header[1] = "T494A";
            else if (header[1] == "E475G")
                header[1] = "E481G";
            else if (header[1] == "T445A")
                header[1] = "T451A";
        }
        else if (header[2] == "ARO:3003043") {          //MEG_5407
            if (header[1] == "M400T" || header[1] == "M339F")
                continue;
        }
        else if (header[2] == "ARO:3003382") {          //MEG_6548
            if (header[1] == "G121D")
                continue;
        }
        else if (header[2] == "ARO:3003369") {          //MEG_7301
            if (header[1] == "R231V" || header[1] == "R234S" || header[1] == "R234F")
                continue;
            else if (header[1] == "T336A")
                header[1] = "T335A";
        }
        if (header[4] != "NA" && header[4] != "not corresponding") {
            if (header[1].find(";") != -1) {
                vector<string> temp = { header[1], header[2], header[4] };
                mult.push_back(temp);
            }
            else
                this->addSNP(header[1], header[2], header[4]);
        }
    }
    snpsearch.close();

    snpsearch.open("KargvaInfoExtracion/SNP_correction.csv");
    while (std::getline(snpsearch, line)) {
        string id = line.substr(0, line.find(","));
        line = line.substr(line.find(","));
        string meg = line.substr(0, line.find(","));
        vector<string> snp;
        while (line.find(",") != -1) {
            line = line.substr(line.find(","));
            if (line[0] == ',') break;
            snp.push_back(line.substr(0, line.find(",")));
        }
        for (auto iter = snp.begin(); iter != snp.end(); ++iter) {
            if ((*iter).find(";") != -1) {
                vector<string> temp = { *iter, id, meg };
                mult.push_back(temp);
            }
            else 
                this->addSNP(*iter, id, meg);
        }
    }
    for (auto iter = mult.begin(); iter != mult.end(); ++iter) {
        this->addSNP((* iter)[0], (* iter)[1], (* iter)[2]);
    }
}