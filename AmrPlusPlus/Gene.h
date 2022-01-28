#pragma once
#include <string>

using namespace std;

class Gene
{
	private:
		string geneName;
		string* geneSequence;
	public:
		Gene(string gName, string gSequence);
		string* getSequence();
};

Gene::Gene(string gName, string gSequence)
{
	geneName = gName;
	geneSequence = &gSequence;
}
string* Gene::getSequence()
{
	return geneSequence;
}