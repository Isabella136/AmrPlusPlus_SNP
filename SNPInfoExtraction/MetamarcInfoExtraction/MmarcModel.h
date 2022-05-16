#pragma once
#include <string>
#include <list>
#include "../ModelReg.h"

using namespace std;

class MmarcModel: public virtual ModelReg {
    private:
        int mmarc_codon_start;
        int mmarc_codon_end;
        string context_left_aa = "";
        string context_right_aa = "";
        void makeModel(string line);
    public:
        MmarcModel(string line);
        ~MmarcModel();
        MmarcModel(const MmarcModel& other);
        Model* Clone();
        string condensedSNPinfo();
};
MmarcModel::MmarcModel(string line)
{
    makeModel(line);
}
MmarcModel::MmarcModel(const MmarcModel& other) {
    this->wt_aa = other.wt_aa;
    this->pos = other.pos;
    this->mutant_aa = other.mutant_aa;
    this->context_left_aa = other.context_left_aa;
    this->context_right_aa = other.context_right_aa;
    this->mmarc_codon_end = other.mmarc_codon_end;
    this->mmarc_codon_start = other.mmarc_codon_start;
}
Model* MmarcModel::Clone() {
    return new MmarcModel(*this);
}
MmarcModel::~MmarcModel() {}
void MmarcModel::makeModel(string line)
{
    line = line.substr(line.find(',') + 1);
    line = line.substr(line.find(',') + 1);
    line = line.substr(line.find(',') + 1);
    line = line.substr(line.find(',') + 1);
    line = line.substr(line.find(',') + 1);
    pos = stoi(line.substr(line.find('-') + 1, line.find(':') - line.find('-') - 1));
    line = line.substr(line.find(',') + 1);
    wt_aa = line[0];
    line = line.substr(line.find(',') + 1);
    string aaList = line.substr(0, line.find(','));
    for (int i = 0; i < aaList.size(); i++)
        mutant_aa.push_back(aaList[i]);
    line = line.substr(line.find(',') + 1);
    mmarc_codon_start = stoi(line.substr(0, line.find(',')));
    line = line.substr(line.find(',') + 1);
    mmarc_codon_end = stoi(line.substr(0, line.find(',')));
    line = line.substr(line.find(',') + 1);
    aaList = line.substr(0, line.find(','));
    for (int i = 0; i < aaList.size(); i++)
        context_left_aa += aaList[i];
    line = line.substr(line.find(',') + 1);
    aaList = line.substr(0, line.find(','));
    for (int i = 0; i < aaList.size(); i++)
        context_right_aa += aaList[i];
}
string MmarcModel::condensedSNPinfo()
{
    string toReturn = "Reg:";
    toReturn += wt_aa;
    toReturn += to_string(pos);
    for (int i = 0; i < mutant_aa.size(); i++)
        toReturn += mutant_aa[i];
    toReturn += '_';
    toReturn += context_left_aa;
    toReturn += '_';
    toReturn += context_right_aa;
    return toReturn;
}
