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
		virtual string condensedInfo() = 0;
};

ModelSNP::ModelSNP() {}
ModelSNP::~ModelSNP() {}