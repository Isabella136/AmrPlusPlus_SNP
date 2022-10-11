#pragma once
#include <string>
#include "Model.h"

using namespace std;

class ModelSNP : public virtual Model{
	protected:
		int pos = 0;
		virtual void makeModel(string line) = 0;
		virtual ~ModelSNP();
	public:
		ModelSNP();
		list<int> getFirstPos();
		virtual string condensedInfo() = 0;
};

ModelSNP::ModelSNP() {}
ModelSNP::~ModelSNP() {}
list<int> ModelSNP::getFirstPos() {
	return { 1, pos };				//SNP come after intrinsic(0) and before hypersusceptible(2)
}