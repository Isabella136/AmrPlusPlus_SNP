#pragma once
#include <string>
#include <unordered_map>
#include "GeneGroup.h"

using namespace std;

class GeneMechanism
{
	private:
		string geneMechanism;
		unordered_map<string, GeneGroup*> geneGroups;
		GeneGroup* getGeneGroup(string geneGroup);
	public:
		GeneMechanism(string geneName, string geneSequence, int delimiter3, int delimiter4);
		Gene* getGene(string geneName, int delimiter4);
		void addGene(string geneName, string geneSequence, int delimiter4);
};

GeneMechanism::GeneMechanism(string geneName, string geneSequence, int delimiter3, int delimiter4)
{
	geneMechanism = geneName.substr(delimiter3 + 1, delimiter4 - delimiter3 - 1);
	addGene(geneName, geneSequence, delimiter4);
}

GeneGroup* GeneMechanism::getGeneGroup(string geneGroup)
{
	try 
	{
		return geneGroups.at(geneGroup);
	}
	catch (const out_of_range& oor)
	{
		return nullptr;
	}
}

Gene* GeneMechanism::getGene(string geneName, int delimiter4)
{
	int delimiter5 = geneName.substr(delimiter4 + 1).find('|');
	string geneGroup= geneName.substr(delimiter4 + 1, delimiter5 - delimiter4 - 1);
	return getGeneGroup(geneGroup)->getGene(geneName);
}

void GeneMechanism::addGene(string geneName, string geneSequence, int delimiter4)
{
	int delimiter5 = geneName.substr(delimiter4 + 1).find('|');
	string geneGroup = geneName.substr(delimiter4 + 1, delimiter5 - delimiter4 - 1);
	if (!getGeneGroup(geneMechanism))
	{
		geneGroups.emplace(geneGroup, new GeneGroup(geneName, geneSequence, delimiter4, delimiter5));
	}
	else
	{
		getGeneGroup(geneGroup)->addGene(geneName, geneSequence);
	}
}