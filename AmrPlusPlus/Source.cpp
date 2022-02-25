#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

using namespace std;


void filterSAMfile(string fileName);

int main(){
	string fileName;
	cout << "Specify file name or type x to exit: \n";
	cin >> fileName;
	while (fileName != "x")
	{
		filterSAMfile(fileName);
		cout << "Specify file name or type x to exit: \n";
		cin >> fileName;
	}
	return 0;
}

void filterSAMfile(string fileName)
{
	ifstream samFileByCategory;
	ifstream samFileByLine;
	samFileByCategory.open(fileName);
	samFileByLine.open(fileName);
	ofstream outputSamFile;
	outputSamFile.open("Filtered_out_" + fileName, ofstream::trunc);
	string QNAME;
	string FLAG;
	string RNAME;
	string toWrite;
	string restOfLine;
	while (std::getline(samFileByCategory, QNAME, '\t'))
	{
		std::getline(samFileByLine, toWrite);
		if (QNAME.at(0) == '@')
		{
			outputSamFile << toWrite << '\n';
		}
		else
		{
			std::getline(samFileByCategory, FLAG, '\t');
			std::getline(samFileByCategory, RNAME, '\t');
			if (RNAME.find("RequiresSNPConfirmation") != string::npos)
			{
				outputSamFile << toWrite << '\n';
			}
		}
		std::getline(samFileByCategory, restOfLine);
	}
	samFileByCategory.close();
	samFileByLine.close();
	outputSamFile.close();
}