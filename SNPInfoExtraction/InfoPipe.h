#pragma once
#include <string>

using namespace std;

class InfoPipe {
	public:
		InfoPipe();
		virtual ~InfoPipe();
		virtual InfoPipe* Clone() = 0;
		virtual string infoType() = 0;
		virtual string condensedInfo() = 0;
};

InfoPipe::InfoPipe() {}
InfoPipe::~InfoPipe() {}
