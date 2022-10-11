#pragma once
#include "KargvaAaMissense.h"

using namespace std;

class KargvaAaHypersuscetible : public virtual KargvaModel {
	protected:
		vector<KargvaModel*> models;
	public:
		KargvaAaHypersuscetible(string line, string id, shared_ptr<CARD_database> dbSeq);
		~KargvaAaHypersuscetible();
		KargvaAaHypersuscetible(const KargvaAaHypersuscetible& other);
		InfoPipe* Clone();
		void addToModel(string line);
		bool includes(string line);
		list<int> getFirstPos();
		string condensedInfo();
		string infoType();
};
list<int> KargvaAaHypersuscetible::getFirstPos() {
	return {2, *(++(models[0]->getFirstPos()).begin())};		//hyper come after snp(1) and before deletion(3)
}
KargvaAaHypersuscetible::KargvaAaHypersuscetible(string line, string id, shared_ptr<CARD_database> dbSeq) {
	string temp = line.substr(4);
	vector<string> snp;
	snp.push_back(temp.substr(0, temp.find(";")));
	while (temp.find(";") != -1) {
		temp = temp.substr(temp.find(";") + 1);
		snp.push_back(temp.substr(0, temp.find(";")));
	}
	for (int i = 0; i < snp.size(); i++) {
		KargvaModel* model;
		model = new KargvaAaMissense(snp[i], id, dbSeq);
		models.push_back(model);
	}
}
KargvaAaHypersuscetible::~KargvaAaHypersuscetible() {}
KargvaAaHypersuscetible::KargvaAaHypersuscetible(const KargvaAaHypersuscetible& other) {
	for (auto iter = other.models.begin(); iter != other.models.end(); ++iter) {
		this->models.push_back(dynamic_cast<KargvaModel*>((*iter)->Clone()));
	}
}
InfoPipe* KargvaAaHypersuscetible::Clone() {
	return new KargvaAaHypersuscetible(*this);
}
void KargvaAaHypersuscetible::addToModel(string line) {
	throw std::exception("should not have been called: model type is hypersusceptible");
}
bool KargvaAaHypersuscetible::includes(string line) {
	return false;
}
string KargvaAaHypersuscetible::condensedInfo() {
	string toReturn = "Hyper";
	for (auto iter = models.begin(); iter != models.end(); ++iter) {
		if (iter == models.begin()) {
			toReturn += ":";
		}
		else {
			toReturn += ";";
		}
		toReturn += (*iter)->condensedInfo();
	}
	return toReturn;
}
string KargvaAaHypersuscetible::infoType() {
	return "Model";
}
