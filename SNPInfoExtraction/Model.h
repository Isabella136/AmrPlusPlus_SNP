#pragma once
#include <string>
#include <list>
#include "InfoPipe.h"

using namespace std;

class Model : public virtual InfoPipe{
	public:
		Model();
		virtual ~Model();
		virtual list<int> getFirstPos() = 0;
		virtual string condensedInfo() = 0;
};

Model::Model() {}
Model::~Model() {}
