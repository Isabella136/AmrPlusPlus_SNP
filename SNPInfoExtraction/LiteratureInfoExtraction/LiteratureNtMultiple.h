#pragma once
#include <string>
#include <vector>
#include "LiteratureNtDeletion.h"
#include "LiteratureNtMissense.h"

using namespace std;

class LiteratureNtMultiple : public virtual LiteratureModel {
protected:
	vector<LiteratureModel*> models;
public:
	LiteratureNtMultiple();
	LiteratureNtMultiple(string line, string id, shared_ptr<MEGARes_database> dbSeq);
	~LiteratureNtMultiple();
	LiteratureNtMultiple(const LiteratureNtMultiple& other);
	InfoPipe* Clone();
	string condensedInfo();
	int getPos();
	string getFirstPos();
	string infoType();

};
string LiteratureNtMultiple::getFirstPos() {
	return "mult" + models[0]->getFirstPos();
}
LiteratureNtMultiple::LiteratureNtMultiple() {}
LiteratureNtMultiple::LiteratureNtMultiple(string line, string id, shared_ptr<MEGARes_database> dbSeq) {
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
			model = new LiteratureNtDeletion(snp[i].substr(7), id, dbSeq);
		else
			model = new LiteratureNtMissense(snp[i].substr(4), id, dbSeq);
		models.push_back(model);
	}
}
LiteratureNtMultiple::LiteratureNtMultiple(const LiteratureNtMultiple& other) {
	for (auto iter = other.models.begin(); iter != other.models.end(); ++iter) {
		this->models.push_back(dynamic_cast<LiteratureModel*>((*iter)->Clone()));
	}
}
InfoPipe* LiteratureNtMultiple::Clone() {
	return new LiteratureNtMultiple(*this);
}
LiteratureNtMultiple::~LiteratureNtMultiple() {
	while (!models.empty())
		models.pop_back();
}

int LiteratureNtMultiple::getPos() {
	throw std::exception("should not have been called: model type is multiple");
}

string LiteratureNtMultiple::condensedInfo()
{
	string toReturn = "NucMult:";
	for (int i = 0; i < models.size(); i++) {
		toReturn += models[i]->condensedInfo();
		if (i + 1 < models.size())
			toReturn += ";";
	}
	return toReturn;
}
string LiteratureNtMultiple::infoType() {
	return "Model";
}
