#pragma once
#include <string>

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
		bool ofInterest = false;

	public:
		Gene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
		string getSequence();
		string getName();
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
string Gene::getSequence()
{
	return geneSequence;
}

string Gene::getName()
{
	return geneName;
}