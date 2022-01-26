#pragma once
#include <string>
#include <unordered_map>
#include "GeneGroup.h"

using namespace std;

class GeneMechanism
{
	private:
		GeneType* geneType;
		GeneClass* geneClass;
		string geneMechanism;
		unordered_map<string, GeneGroup*> geneGroups;
	public:
		GeneMechanism(string geneName, string geneSequence, int delimiterIndexBeforeMechanism);
		GeneGroup* getGeneGroup(string geneGroup);
		Gene* getGene(string geneName);
		void setGeneType(GeneType& gType);
		void setGeneClass(GeneClass& gClass);
		void addGeneGroup(string geneName, string geneSequence, int delimiterIndexBeforeGroup);
};

GeneMechanism::GeneMechanism(string geneName, string geneSequence, int delimiterIndexBeforeMechanism)
{

}