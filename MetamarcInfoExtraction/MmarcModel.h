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
    string context_left_aa = "";
    string context_right_aa = "";
public:
    MmarcModel(string line);
    pair<pair<char, int>, list<char>> condensedSNPinfo();
    pair<int, int> codonStartAndEnd();
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
        context_left_aa += aaList[i];
    line = line.substr(line.find(',') + 1);
    aaList = line.substr(0, line.find(','));
    for (int i = 0; i < aaList.size(); i++)
        context_right_aa += aaList[i];
}
pair<pair<char, int>, list<char>> MmarcModel::condensedSNPinfo()
{
    list<char> toReturn = mutant_aa;
    toReturn.push_back('_');
    for (int i = 0; i < context_left_aa.length(); i++)
        toReturn.push_back(context_left_aa[i]);
    toReturn.push_back('_');
    for (int i = 0; i < context_right_aa.length(); i++)
        toReturn.push_back(context_right_aa[i]);
    return make_pair(make_pair(wt_aa, pos), toReturn);
}
pair<int, int> MmarcModel::codonStartAndEnd()
{
    return make_pair(mmarc_codon_start, mmarc_codon_end);
}