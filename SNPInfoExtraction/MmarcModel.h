#pragma once
#include <string>
#include <list>

using namespace std;

class MmarcModel
{
private:
    int pos;
    char wt_aa;
    list<char> mutant_aa;
    int mmarc_codon_start;
    int mmarc_codon_end;
    char context_left_aa[5] = {};
    char context_right_aa[5] = {};
public:
    MmarcModel(string line);
    pair<pair<char, int>, list<char>> condensedSNPinfo();
    pair<int, int> codonStartAndEnd();
    pair<char*, char*> context();
};
MmarcModel::MmarcModel(string line)
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
        context_left_aa[i] = aaList[i];
    line = line.substr(line.find(',') + 1);
    aaList = line.substr(0, line.find(','));
    for (int i = 0; i < aaList.size(); i++)
        context_right_aa[i] = aaList[i];
}
pair<pair<char, int>, list<char>> MmarcModel::condensedSNPinfo()
{
    return make_pair(make_pair(wt_aa, pos), mutant_aa);
}
pair<int, int> MmarcModel::codonStartAndEnd()
{
    return make_pair(mmarc_codon_start, mmarc_codon_end);
}
pair<char*, char*> MmarcModel::context()
{
    return make_pair(context_left_aa, context_right_aa);
}