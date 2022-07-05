#pragma once
#include <string>
#include <vector>
#include "LiteratureModelDeletion.h"
#include "LiteratureModelMissense.h"

using namespace std;

class LiteratureMultipleModels : public virtual LiteratureModel {
protected:
	vector<LiteratureModel*> models;
public:
	LiteratureMultipleModels();
	LiteratureMultipleModels(string line, string id, shared_ptr<MEGARes_database> dbSeq);
	~LiteratureMultipleModels();
	LiteratureMultipleModels(const LiteratureMultipleModels& other);
	Model* Clone();
	string condensedSNPinfo();
	int getPos();

};
LiteratureMultipleModels::LiteratureMultipleModels() {}
LiteratureMultipleModels::LiteratureMultipleModels(string line, string id, shared_ptr<MEGARes_database> dbSeq) {
	string temp = line;
	vector<string> snp;
	snp.push_back(temp.substr(0, temp.find(";")));
	while (temp.find(";") != -1) {
		temp = temp.substr(temp.find(";") + 1);
		snp.push_back(temp.substr(0, temp.find(";")));
	}
	for (int i = 0; i < snp.size(); i++) {
		LiteratureModel* model;
		if (snp[i].find("Del:") != -1)
			model = new LiteratureModelDeletion(snp[i].substr(4), id, dbSeq);
		else
			model = new LiteratureModelMissense(snp[i].substr(4), id, dbSeq);
		models.push_back(model);
	}
}
LiteratureMultipleModels::LiteratureMultipleModels(const LiteratureMultipleModels& other) {
	for (auto iter = other.models.begin(); iter != other.models.end(); ++iter) {
		this->models.push_back(dynamic_cast<LiteratureModel*>((*iter)->Clone()));
	}
}
Model* LiteratureMultipleModels::Clone() {
	return new LiteratureMultipleModels(*this);
}
int LiteratureMultipleModels::getPos() {
	throw std::exception("should not have been called: model type is multiple");
}
LiteratureMultipleModels::~LiteratureMultipleModels() {
	while (!models.empty())
		models.pop_back();
}
string LiteratureMultipleModels::condensedSNPinfo()
{
	string toReturn = "Mult:";
	for (int i = 0; i < models.size(); i++) {
		toReturn += models[i]->condensedSNPinfo();
		if (i + 1 < models.size())
			toReturn += ";";
	}
	return toReturn;
}
