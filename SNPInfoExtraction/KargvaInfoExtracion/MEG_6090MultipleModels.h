#pragma once
#include <string>
#include <vector>
#include "KargvaMultipleModels.h"

using namespace std;

class MEG_6090MultipleModels : public virtual KargvaMultipleModels {
	public:
		MEG_6090MultipleModels(string line, string id, CARD_database* dbSeq);
};

MEG_6090MultipleModels::MEG_6090MultipleModels(string line, string id, CARD_database* dbSeq) {
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
		else if (snp[i].find("STOP") != -1)
			model = new KargvaModelNonsense(snp[i], id, dbSeq);
		else
			model = new KargvaModelReg(snp[i], id, dbSeq);
		models.push_back(model);
	}
}