#pragma once
#include <string>
#include <unordered_map>
#include "Gene.h"

using namespace std;

class GeneGroup
{
	private:
		GeneType* geneType;
		GeneClass* geneClass;
		GeneMechanism* geneMechanism;
		string geneGroup;
		unordered_map<string, Gene*> genes;
	public:
		GeneGroup(string geneName, string geneSequence, int delimiterIndexBeforeGroup);
		Gene* getGene(string geneName);
		void setGeneType(GeneType& gType);
		void setGeneClass(GeneClass& gClass);
		void setGeneMechanism(GeneMechanism& gMechanisms);
		void addGene(string geneName, string geneSequence);
};