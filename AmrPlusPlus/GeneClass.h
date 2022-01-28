#pragma once
#include <string>
#include <unordered_map>
#include "GeneMechanism.h"

using namespace std;

class GeneClass
{
	private:
		string geneClass;
		unordered_map<string, GeneMechanism*> geneMechanisms;
		GeneMechanism* getGeneMechanism(string geneMechanism); 
	public:
		GeneClass(string geneName, string geneSequence, int delimiter2, int delimiter3);
		Gene* getGene(string geneName, int delimiter3);
		void addGene(string geneName, string geneSequence, int delimiter3);
};

GeneClass::GeneClass(string geneName, string geneSequence, int delimiter2, int delimiter3)
{
	geneClass = geneName.substr(delimiter2 + 1, delimiter3 - delimiter2 - 1);
	addGene(geneName, geneSequence, delimiter3);
}

GeneMechanism* GeneClass::getGeneMechanism(string geneMechanism)
{
	try 
	{
		return geneMechanisms.at(geneMechanism);
	}
	catch (const out_of_range& oor)
	{
		return nullptr;
	}
}

Gene* GeneClass::getGene(string geneName, int delimiter3)
{
	int delimiter4 = geneName.substr(delimiter3 + 1).find('|');
	string geneMechanism = geneName.substr(delimiter3 + 1, delimiter4 - delimiter3 - 1);
	return getGeneMechanism(geneMechanism)->getGene(geneName, delimiter4);
}

void GeneClass::addGene(string geneName, string geneSequence, int delimiter3)
{
	int delimiter4 = geneName.substr(delimiter3 + 1).find('|');
	string geneMechanism = geneName.substr(delimiter3 + 1, delimiter4 - delimiter3 - 1);
	if (!getGeneMechanism(geneMechanism))
	{
		geneMechanisms.emplace(geneMechanism, new GeneMechanism(geneName, geneSequence, delimiter3, delimiter4));
	}
	else
	{
		getGeneMechanism(geneMechanism)->addGene(geneName, geneSequence, delimiter4);
	}
}