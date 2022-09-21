#pragma once
#include <string>
#include <vector>
#include "KargvaAaDeletion.h"
#include "KargvaAaMissense.h"
#include "KargvaAaNonstop.h"
#include "KargvaAaInsertion.h"

using namespace std;

class KargvaAaMultiple : public virtual KargvaModel {
	protected:
		vector<KargvaModel*> models;
	public:
		KargvaAaMultiple();
		KargvaAaMultiple(string line, string id, shared_ptr<CARD_database> dbSeq);
		~KargvaAaMultiple();
		KargvaAaMultiple(const KargvaAaMultiple& other);
		InfoPipe* Clone();
		void addToModel(string line);
		bool includes(string line);
		string condensedInfo();
		string infoType();

};
KargvaAaMultiple::KargvaAaMultiple() {}
KargvaAaMultiple::KargvaAaMultiple(string line, string id, shared_ptr<CARD_database> dbSeq) {
	string temp = line;
	vector<string> snp;
	snp.push_back(temp.substr(0, temp.find(";")));
	while (temp.find(";") != -1) {
		temp = temp.substr(temp.find(";") + 1);
		snp.push_back(temp.substr(0, temp.find(";")));
	}
	for (int i = 0; i < snp.size(); i++) {
		KargvaModel* model;
		if (snp[i] == "nonstop")
			model = new KargvaAaNonstop(snp[i], id, dbSeq);
		else if (snp[i].find("-") != -1)
			model = new KargvaAaDeletion(snp[i], id, dbSeq);
		else if (snp[i].find("+") != -1)
			model = new KargvaAaInsertion(snp[i], id, dbSeq);
		else
			model = new KargvaAaMissense(snp[i], id, dbSeq);
		models.push_back(model);
	}
}
KargvaAaMultiple::KargvaAaMultiple(const KargvaAaMultiple& other) {
	for (auto iter = other.models.begin(); iter != other.models.end(); ++iter) {
		this->models.push_back(dynamic_cast<KargvaModel*>((*iter)->Clone()));
	}
}
InfoPipe* KargvaAaMultiple::Clone() {
	return new KargvaAaMultiple(*this);
}
KargvaAaMultiple::~KargvaAaMultiple() {
	while (!models.empty())
		models.pop_back();
}
void KargvaAaMultiple::addToModel(string line) {
	throw std::exception("should not have been called: model type is multiple");
}
bool KargvaAaMultiple::includes(string line) {
	return false;
}
string KargvaAaMultiple::condensedInfo()
{
	string toReturn = "Mult:";
	for (int i = 0; i < models.size(); i++) {
		toReturn += models[i]->condensedInfo();
		if (i + 1 < models.size())
			toReturn += ";";
	}
	return toReturn;
}
string KargvaAaMultiple::infoType() {
	return "Model";
}