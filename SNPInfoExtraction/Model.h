#pragma once
#include <string>

using namespace std;

class Model {
	public:
		Model();
		virtual ~Model();
		virtual string condensedSNPinfo() = 0;
		virtual Model* Clone() = 0;
};

Model::Model() {}
Model::~Model() {}
