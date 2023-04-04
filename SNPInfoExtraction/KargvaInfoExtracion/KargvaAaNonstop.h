#pragma once
#include "../AaNonstop.h"
#include "KargvaModel.h"

using namespace std;

class KargvaAaNonstop : public AaNonstop, public KargvaModel {
	private:
		void makeModel(string line);
	public:
		KargvaAaNonstop(string line, string id, shared_ptr<CARD_database> dbSeq);
		~KargvaAaNonstop();
		KargvaAaNonstop(const KargvaAaNonstop& other);
		InfoPipe* Clone();
		void addToModel(string line);
		bool includes(string line);
		string condensedInfo();
		string infoType();
};


KargvaAaNonstop::KargvaAaNonstop(string line, string id, shared_ptr<CARD_database> dbSeq):KargvaModel(id, dbSeq){
	makeModel(line);
}
KargvaAaNonstop::~KargvaAaNonstop() {}
KargvaAaNonstop::KargvaAaNonstop(const KargvaAaNonstop& other) {
	this->pos = other.pos;
	this->cardID = other.cardID;
	this->databaseSequences = other.databaseSequences;
}
InfoPipe* KargvaAaNonstop::Clone() {
	return new KargvaAaNonstop(*this);
}
void KargvaAaNonstop::makeModel(string line) {
	pos = findTermination() + 1;
}
void KargvaAaNonstop::addToModel(string line) {
	throw std::exception("should not have been called: model type is nonstop");
}
bool KargvaAaNonstop::includes(string line) {
	return false;
}
string KargvaAaNonstop::condensedInfo() {
	return "Nonstop:" + to_string(pos);
}
string KargvaAaNonstop::infoType() {
	return "Model";
}