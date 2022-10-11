#pragma once
#include <string>
#include "InfoPipe.h"

using namespace std;

class FrameshiftInfo : public InfoPipe {
	protected:
		string detail;
	public:
		FrameshiftInfo(string input, string source = "CARD");
		~FrameshiftInfo();
		FrameshiftInfo(const FrameshiftInfo& other);
		InfoPipe* Clone();
		string condensedInfo();
		string infoType();
};

FrameshiftInfo::FrameshiftInfo(string input, string source) {
	detail = input;
	this->source = source;
}
FrameshiftInfo::FrameshiftInfo(const FrameshiftInfo& other) {
	this->detail = other.detail;
	this->source = other.source;
}
InfoPipe* FrameshiftInfo::Clone() {
	return new FrameshiftInfo(*this);
}
FrameshiftInfo::~FrameshiftInfo() {}

string FrameshiftInfo::infoType() {
	return "Frameshift";
}
string FrameshiftInfo::condensedInfo() {
	return detail;
}