#pragma once
#include <string>
#include <vector>
#include "KargvaMultipleModels.h"

using namespace std;

class MEG_6090MultipleModels : public virtual KargvaMultipleModels {
	public:
		MEG_6090MultipleModels(string line, string id, shared_ptr<CARD_database> dbSeq);
		~MEG_6090MultipleModels();
		MEG_6090MultipleModels(const MEG_6090MultipleModels& other);
		Model* Clone();
};

MEG_6090MultipleModels::MEG_6090MultipleModels(const MEG_6090MultipleModels& other) {
	for (auto iter = other.models.begin(); iter != other.models.end(); ++iter) {
		this->models.push_back(dynamic_cast<KargvaModel*>((*iter)->Clone()));
	}
}
Model* MEG_6090MultipleModels::Clone() {
	return new MEG_6090MultipleModels(*this);
}
MEG_6090MultipleModels::~MEG_6090MultipleModels() {}
MEG_6090MultipleModels::MEG_6090MultipleModels(string line, string id, shared_ptr<CARD_database> dbSeq) {
	string temp = line;
	vector<string> snp;
	snp.push_back(temp.substr(0, temp.find(",")));
	while (temp.find(",") != -1) {
		temp = temp.substr(temp.find(",") + 1);
		snp.push_back(temp.substr(0, temp.find(",")));
	}
	for (int i = 0; i < snp.size(); i++) {
		KargvaModel* model;
		if (snp[i].at(snp[i].size() - 1) == '-')
			model = new KargvaModelDeletion(snp[i], id, dbSeq);
		else
			model = new KargvaModelMissense(snp[i], id, dbSeq);
		models.push_back(model);
	}
}