#pragma once
#include <string>
#include <vector>
#include "LiteratureNucModelDeletion.h"
#include "LiteratureNucModelMissense.h"

using namespace std;

class LiteratureNucMultipleModels : public virtual LiteratureModel {
protected:
	vector<LiteratureModel*> models;
public:
	LiteratureNucMultipleModels();
	LiteratureNucMultipleModels(string line, string id, shared_ptr<MEGARes_database> dbSeq);
	~LiteratureNucMultipleModels();
	LiteratureNucMultipleModels(const LiteratureNucMultipleModels& other);
	Model* Clone();
	string condensedSNPinfo();
	int getPos();

};
LiteratureNucMultipleModels::LiteratureNucMultipleModels() {}
LiteratureNucMultipleModels::LiteratureNucMultipleModels(string line, string id, shared_ptr<MEGARes_database> dbSeq) {
	string temp = line;
	vector<string> snp;
	snp.push_back(temp.substr(0, temp.find(";")));
	while (temp.find(";") != -1) {
		temp = temp.substr(temp.find(";") + 1);
		snp.push_back(temp.substr(0, temp.find(";")));
	}
	for (int i = 0; i < snp.size(); i++) {
		LiteratureModel* model;
		if (snp[i].find("NucDel:") != -1)
			model = new LiteratureNucModelDeletion(snp[i].substr(7), id, dbSeq);
		else
			model = new LiteratureNucModelMissense(snp[i].substr(4), id, dbSeq);
		models.push_back(model);
	}
}
LiteratureNucMultipleModels::LiteratureNucMultipleModels(const LiteratureNucMultipleModels& other) {
	for (auto iter = other.models.begin(); iter != other.models.end(); ++iter) {
		this->models.push_back(dynamic_cast<LiteratureModel*>((*iter)->Clone()));
	}
}
Model* LiteratureNucMultipleModels::Clone() {
	return new LiteratureNucMultipleModels(*this);
}
LiteratureNucMultipleModels::~LiteratureNucMultipleModels() {
	while (!models.empty())
		models.pop_back();
}

int LiteratureNucMultipleModels::getPos() {
	throw std::exception("should not have been called: model type is multiple");
}

string LiteratureNucMultipleModels::condensedSNPinfo()
{
	string toReturn = "NucMult:";
	for (int i = 0; i < models.size(); i++) {
		toReturn += models[i]->condensedSNPinfo();
		if (i + 1 < models.size())
			toReturn += ";";
	}
	return toReturn;
}
