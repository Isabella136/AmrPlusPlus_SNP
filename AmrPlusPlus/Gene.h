#pragma once
#include <string>

using namespace std;

class Gene
{
	private:
		GeneType* geneType;
		GeneClass* geneClass;
		GeneMechanism* geneMechanism;
		GeneGroup* geneGroup;
		string geneName;
		string* geneSequence;
	public:
		Gene(string gName, string gSequence);
		string* getSequence();
		void setGeneType(GeneType& gType);
		void setGeneClass(GeneClass& gClass);
		void setGeneMechanism(GeneMechanism& gMechanisms);
		void setGeneGroup(GeneGroup& gGroup);
};