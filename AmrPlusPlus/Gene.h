#pragma once
#include <string>
#include <unordered_map>
#include <list>
#include "MmarcModel.h"

using namespace std;

struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2>& p) const
    {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

class Gene
{
	private:
		string geneName;
		string geneType;
		string geneClass;
		string geneMechanism;
		string geneGroup;
		string geneSequence;
		bool ofInterest = false;
		unordered_map<pair<char, int>, list<char>, hash_pair> listOfSNPs;

	public:
		Gene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
		string getSequence();
		string getName();
		string getFASTA();
		void makeOfInterest();
		void addSNP(pair<char, int> wt_pos, list<char> mutants);
		bool isOfInterest();
};

Gene::Gene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence) 
{
	geneName = _geneName;
	geneType = _geneType;
	geneClass = _geneClass;
	geneMechanism = _geneMechanism;
	geneGroup = _geneGroup;
	geneSequence = _geneSequence;
}
void Gene::addSNP(pair<char, int> wt_pos, list<char> mutants)
{
	listOfSNPs.emplace(wt_pos, mutants);
}
void Gene::makeOfInterest()
{
	ofInterest = true;
}
bool Gene::isOfInterest()
{
	return ofInterest;
}
string Gene::getSequence()
{
	return geneSequence;
}

string Gene::getName()
{
	return geneName;
}

string Gene::getFASTA()
{
	if (ofInterest) 
	{
		string SNPinfo = "";
		for (auto iter = listOfSNPs.begin(); iter != listOfSNPs.end(); ++iter)
		{
			SNPinfo = SNPinfo + "|";
			SNPinfo = SNPinfo + iter->first.first + to_string(iter->first.second);
			for (int i = 0; i < iter->second.size(); i++) 
			{
				SNPinfo = SNPinfo + iter->second.front();
				iter->second.pop_front();
			}
		}
		return ">" + geneName + "|" + geneType + "|" + geneClass + "|" + geneMechanism + "|" + geneGroup + SNPinfo + "\n" + geneSequence + "\n";
	}
	else
		return "";
}