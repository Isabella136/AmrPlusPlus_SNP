#pragma once
#include <string>
#include <vector>
#include "KargvaModelDeletion.h"
//#include "KargvaModelInsertion.h"
#include "KargvaModelNonsense.h"
#include "KargvaModelReg.h"

using namespace std;

class KargvaMultipleModels : public virtual KargvaModel {
	protected:
		vector<KargvaModel*> models;
	public:
		KargvaMultipleModels(string line, string id, CARD_database* dbSeq);
		void addToModel(string line);
		bool includes(string line);
		string condensedSNPinfo();

};
KargvaMultipleModels::KargvaMultipleModels(string line, string id, CARD_database* dbSeq) {
	string temp = line;
	vector<string> snp;
	snp.push_back(temp.substr(0, temp.find(";")));
	while (temp.find(";") != -1) {
		temp = temp.substr(temp.find(";") + 1);
		snp.push_back(temp.substr(0, temp.find(";")));
	}
	for (int i = 0; i < snp.size(); i++) {
		KargvaModel* model;
		if (snp[i].at(0) == '-')
			model = new KargvaModelDeletion(snp[i], id, dbSeq);
		else if (snp[i].find("STOP") != -1)
			model = new KargvaModelNonsense(snp[i], id, dbSeq);
		else
			model = new KargvaModelReg(snp[i], id, dbSeq);
		models.push_back(model);
	}
}
void KargvaMultipleModels::addToModel(string line) {
	throw std::exception("should not have been called: model type is multiple");
}
bool KargvaMultipleModels::includes(string line) {
	return false;
}
string KargvaMultipleModels::condensedSNPinfo()
{
	string toReturn = "Mult:";
	for (int i = 0; i < models.size(); i++) {
		toReturn += models[i]->condensedSNPinfo();
		if (i + 1 < models.size())
			toReturn += ";";
	}
	return toReturn;
}