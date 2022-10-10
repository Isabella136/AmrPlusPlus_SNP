#pragma once
#include <string>
#include <vector>
#include "LiteratureAaDeletion.h"
#include "LiteratureAaMissense.h"

using namespace std;

class LiteratureAaMultiple : public virtual LiteratureModel {
protected:
	vector<LiteratureModel*> models;
public:
	LiteratureAaMultiple();
	LiteratureAaMultiple(string line, string id, shared_ptr<MEGARes_database> dbSeq);
	~LiteratureAaMultiple();
	LiteratureAaMultiple(const LiteratureAaMultiple& other);
	InfoPipe* Clone();
	string condensedInfo();
	int getPos();
	string getFirstPos();
	string infoType();

};
string LiteratureAaMultiple::getFirstPos() {
	return "mult" + models[0]->getFirstPos();
}
LiteratureAaMultiple::LiteratureAaMultiple() {}
LiteratureAaMultiple::LiteratureAaMultiple(string line, string id, shared_ptr<MEGARes_database> dbSeq) {
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
			model = new LiteratureAaDeletion(snp[i].substr(4), id, dbSeq);
		else
			model = new LiteratureAaMissense(snp[i].substr(4), id, dbSeq);
		models.push_back(model);
	}
}
LiteratureAaMultiple::LiteratureAaMultiple(const LiteratureAaMultiple& other) {
	for (auto iter = other.models.begin(); iter != other.models.end(); ++iter) {
		this->models.push_back(dynamic_cast<LiteratureModel*>((*iter)->Clone()));
	}
}
InfoPipe* LiteratureAaMultiple::Clone() {
	return new LiteratureAaMultiple(*this);
}
int LiteratureAaMultiple::getPos() {
	throw std::exception("should not have been called: model type is multiple");
}
LiteratureAaMultiple::~LiteratureAaMultiple() {
	while (!models.empty())
		models.pop_back();
}
string LiteratureAaMultiple::condensedInfo()
{
	string toReturn = "Mult:";
	for (int i = 0; i < models.size(); i++) {
		toReturn += models[i]->condensedInfo();
		if (i + 1 < models.size())
			toReturn += ";";
	}
	return toReturn;
}
string LiteratureAaMultiple::infoType() {
	return "Model";
}