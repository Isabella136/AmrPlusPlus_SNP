#pragma once
#include <string>
#include <unordered_map>
#include "Gene.h"

using namespace std;

class GeneMechanism
{
	private:
		string geneMechanism;
		unordered_map<string,Gene*> genes;
	public:
		GeneMechanism(string _geneMechanism, Gene* gene);
		Gene* getGene(string geneName);
		void addGene(Gene* gene);
};

GeneMechanism::GeneMechanism(string _geneMechanism, Gene* gene)
{
	geneMechanism = _geneMechanism;
	addGene(gene);
}

Gene* GeneMechanism::getGene(string geneName)
{
	try 
	{
		return genes.at(geneName);
	}
	catch (const out_of_range& oor)
	{
		return nullptr;
	}
}

void GeneMechanism::addGene(Gene* gene)
{
	genes.emplace(gene->getName(), gene);
}