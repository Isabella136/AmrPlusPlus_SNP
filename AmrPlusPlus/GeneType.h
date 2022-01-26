#pragma once
#include <string>
#include <unordered_map>
#include "GeneClass.h"

using namespace std;

class GeneType
{
	private:
		string geneType;
		unordered_map<string,GeneClass*> geneClasses;
	public:
		GeneType(string geneName, string geneSequence);
		GeneClass* getGeneClass(string geneClass);
		Gene* getGene(string geneName);
		void addGeneClass(string geneName, string geneSequence, int delimiterIndexBeforeClass);
};

GeneType::GeneType(string geneName, string geneSequence)
{
	int firstDelimiterIndex = geneName.find('|');
	int nextDelimiterIndex = geneName.substr(firstDelimiterIndex + 1).find('|');
	geneType = geneName.substr(firstDelimiterIndex + 1, nextDelimiterIndex - firstDelimiterIndex - 1);
	addGeneClass(geneName, geneSequence, nextDelimiterIndex);
}

void GeneType::addGeneClass(string geneName, string geneSequence, int delimiterIndexBeforeClass)
{
	int nextDelimiterIndex = geneName.substr(delimiterIndexBeforeClass + 1).find('|');
	string geneClass = geneName.substr(delimiterIndexBeforeClass + 1, nextDelimiterIndex - delimiterIndexBeforeClass - 1);
	geneClasses.emplace(geneClass, new GeneClass(geneName, geneSequence, delimiterIndexBeforeClass));
	geneClasses[geneClass]->setGeneType(this);
}

GeneClass* GeneType::getGeneClass(string geneClass)
{
	return geneClasses[geneClass];
}

Gene* GeneType::getGene(string geneName)
{
	int firstDelimiterIndex = geneName.find('|');
	int delimiterIndexBeforeClass = geneName.substr(firstDelimiterIndex + 1).find('|');
	int nextDelimiterIndex = geneName.substr(delimiterIndexBeforeClass + 1).find('|');
	string geneClass = geneName.substr(delimiterIndexBeforeClass + 1, nextDelimiterIndex - delimiterIndexBeforeClass - 1);
	return getGeneClass(geneClass)->getGene(geneName);
}