#pragma once
#include <string>
#include "Model.h"

using namespace std;

class ModelSingle : public virtual Model{
	protected:
		int pos = 0;
		virtual void makeModel(string line) = 0;
		virtual ~ModelSingle();
	public:
		ModelSingle();
		virtual string condensedSNPinfo() = 0;
};

ModelSingle::ModelSingle() {}
ModelSingle::~ModelSingle() {}