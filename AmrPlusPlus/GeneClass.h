#pragma once
#include <string>
#include <unordered_map>
#include "GeneMechanism.h"

using namespace std;

class GeneClass
{
	private:
		GeneType* geneType;
		string geneClass;
		unordered_map<string, GeneMechanism*> geneMechanisms;
	public:
		GeneClass(string geneName, string geneSequence, int delimiterIndexBeforeClass);
		GeneMechanism* getGeneMechanism(string geneMechanism); 
		Gene* getGene(string geneName);
		GeneType* getGeneType();
		void setGeneType(GeneType& gType);
		void addGeneMechanism(string geneName, string geneSequence, int delimiterIndexBeforeMechanism);
};

GeneClass::GeneClass(string geneName, string geneSequence, int delimiterIndexBeforeClass)
{
	int nextDelimiterIndex = geneName.substr(delimiterIndexBeforeClass + 1).find('|');
	geneClass = geneName.substr(delimiterIndexBeforeClass + 1, nextDelimiterIndex - delimiterIndexBeforeClass - 1);
	//geneMechanisms.push_back(new GeneMechanism(geneName, geneSequence, nextDelimiterIndex, this));
}