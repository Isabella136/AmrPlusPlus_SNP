#pragma once
#include <string>
#include <vector>
#include "KargvaModelDeletion.h"
#include "KargvaModelInsertion.h"
#include "KargvaModelNonsense.h"
#include "KargvaModelReg.h"

using namespace std;

class KargvaMultipleModels : public virtual Model {
	protected:
		vector<KargvaModel*> models;
	public:
		KargvaMultipleModels(string line, string id);
		string condensedSNPinfo();

};
KargvaMultipleModels::KargvaMultipleModels(string line, string id) {

}
string KargvaMultipleModels::condensedSNPinfo()
{
	string toReturn = "Mult:";
}