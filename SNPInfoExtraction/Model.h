#pragma once
#include <string>
#include "InfoPipe.h"

using namespace std;

class Model : public virtual InfoPipe{
	public:
		Model();
		virtual ~Model();
		virtual string condensedInfo() = 0;
};

Model::Model() {}
Model::~Model() {}
