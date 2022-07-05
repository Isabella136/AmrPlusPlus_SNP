#pragma once
#include <string>
#include <vector>
#include "KargvaModelDeletion.h"
#include "KargvaModelNonsense.h"
#include "KargvaModelMissense.h"

using namespace std;

class KargvaMultipleModels : public virtual KargvaModel {
	protected:
		vector<KargvaModel*> models;
	public:
		KargvaMultipleModels();
		KargvaMultipleModels(string line, string id, shared_ptr<CARD_database> dbSeq);
		~KargvaMultipleModels();
		KargvaMultipleModels(const KargvaMultipleModels& other);
		Model* Clone();
		void addToModel(string line);
		bool includes(string line);
		string condensedSNPinfo();

};
KargvaMultipleModels::KargvaMultipleModels() {}
KargvaMultipleModels::KargvaMultipleModels(string line, string id, shared_ptr<CARD_database> dbSeq) {
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
		else
			model = new KargvaModelMissense(snp[i], id, dbSeq);
		models.push_back(model);
	}
}
KargvaMultipleModels::KargvaMultipleModels(const KargvaMultipleModels& other) {
	for (auto iter = other.models.begin(); iter != other.models.end(); ++iter) {
		this->models.push_back(dynamic_cast<KargvaModel*>((*iter)->Clone()));
	}
}
Model* KargvaMultipleModels::Clone() {
	return new KargvaMultipleModels(*this);
}
KargvaMultipleModels::~KargvaMultipleModels() {
	while (!models.empty())
		models.pop_back();
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