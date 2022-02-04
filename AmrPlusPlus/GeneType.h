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
		GeneClass* getGeneClass(string geneClass);
	public:
		GeneType(string geneName, string geneSequence);
		Gene* getGene(string geneName);
		void addGene(string geneName, string geneSequence, int delimiter2);
};

GeneType::GeneType(string geneName, string geneSequence)
{
	int delimiter1 = geneName.find('|');
	int delimiter2 = geneName.substr(delimiter1 + 1).find('|') + delimiter1 + 1;
	geneType = geneName.substr(delimiter1 + 1, delimiter2 - delimiter1 - 1);
	addGene(geneName, geneSequence, delimiter2);
}

GeneClass* GeneType::getGeneClass(string geneClass)
{
	try 
	{
		return geneClasses.at(geneClass);
	}
	catch (const out_of_range& oor)
	{
		return nullptr;
	}
}

Gene* GeneType::getGene(string geneName)
{
	int delimiter1 = geneName.find('|');
	int delimiter2 = geneName.substr(delimiter1 + 1).find('|') + delimiter1 + 1;
	int delimiter3 = geneName.substr(delimiter2 + 1).find('|') + delimiter2 + 1;
	string geneClass = geneName.substr(delimiter2 + 1, delimiter3 - delimiter2 - 1);
	return getGeneClass(geneClass)->getGene(geneName, delimiter3);
}

void GeneType::addGene(string geneName, string geneSequence, int delimiter2)
{
	int delimiter3 = geneName.substr(delimiter2 + 1).find('|') + delimiter2 + 1;
	string geneClass = geneName.substr(delimiter2 + 1, delimiter3 - delimiter2 - 1);
	if (!getGeneClass(geneClass))
	{
		geneClasses.emplace(geneClass, new GeneClass(geneName, geneSequence, delimiter2, delimiter3));
	}
	else
	{
		getGeneClass(geneClass)->addGene(geneName, geneSequence, delimiter3);
	}
}