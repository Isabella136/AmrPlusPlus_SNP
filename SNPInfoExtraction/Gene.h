#pragma once
#include <string>
#include <unordered_map>
#include <list>

using namespace std;

class Gene
{
private:
	string geneName;
	string geneType;
	string geneClass;
	string geneMechanism;
	string geneGroup;
	string geneSequence;
	vector<string> listOfSNPs;

public:
	Gene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
	string getSequence();
	string getName();
	string getFASTA();
	string getFASTA2();
	void addSNP(string SNP);
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
void Gene::addSNP(string SNP)
{
	listOfSNPs.push_back(SNP);
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
	string SNPinfo = "";
	if (listOfSNPs.empty())
		return "";
	else
	{
		for (auto iter = listOfSNPs.begin(); iter != listOfSNPs.end(); ++iter)
		{
			SNPinfo = SNPinfo + "|" + *iter;
		}
	}
	return ">" + geneName + "|" + geneType + "|" + geneClass + "|" + geneMechanism + "|" + geneGroup + SNPinfo + "\n" + geneSequence + "\n";

}

string Gene::getFASTA2()
{
	if (!(listOfSNPs.empty()))
		return "";
	return geneName + "," + geneType + "," + geneClass + "," + geneMechanism + "," + geneGroup + "," + "RequiresSNPConfirmation" + "\n";

}