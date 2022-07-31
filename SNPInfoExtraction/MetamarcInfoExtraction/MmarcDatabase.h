#pragma once
#include <fstream>
#include <string>
#include <unordered_map>
#include <iostream>
#include "MmarcAaMissense.h"
#include "../ModelDatabase.h"

using namespace std;

class MmarcDatabase : public ModelDatabase{
    public:
        MmarcDatabase();
        void SNPInfo();
};

MmarcDatabase::MmarcDatabase()
{
    SNPInfo();
}

void MmarcDatabase::SNPInfo()
{
    ifstream snpsearch;
    snpsearch.open("MetamarcInfoExtraction/metamarc_files/mmarc_snpsearch_metadata2_modified.txt");
    unordered_map<string, list<InfoPipe*>> name_model;
    string line;
    std::getline(snpsearch, line);
    while (std::getline(snpsearch, line))
    {
        string modelName = line.substr(0, line.find(','));
        InfoPipe* model = new MmarcAaMissense(line);
        if (name_model.find(modelName) != name_model.end())
            name_model.at(modelName).push_back((model->Clone()));
        else { 
			list<InfoPipe*> temp;
			temp.push_back(model->Clone());
			name_model.emplace(modelName, temp);
		}
        delete model;
    }
    snpsearch.close();
    ifstream model_members;
    model_members.open("MetamarcInfoExtraction/metamarc_files/mmarc_model_members.csv");
    unordered_map<string, string> header_name;
    std::getline(model_members, line);
    while (std::getline(model_members, line))
    {
        string modelName = line.substr(0, line.find('\t'));
        line = line.substr(line.find('\t') + 1);
            header_name.emplace(line.substr(0, line.find(',')), modelName);
        while (line.find(',') != -1)
        {
            line = line.substr(line.find(',') + 1);
            header_name.emplace(line.substr(0, line.find(',')), modelName);
            
        }
    }
    model_members.close();
    unordered_map<string, list<InfoPipe*>> header_model;
    for (auto iter = header_name.begin(); iter != header_name.end(); ++iter)
    {
        if (name_model.find(iter->second) != name_model.end()) {
            header_model.emplace(iter->first, name_model.at(iter->second)); 
        }
    }
    ifstream v1;
    v1.open("../mapping_files/megaresv1_to_external_header_mappings_v1.01.csv");
    unordered_map<string, string> source_to_header;
    std::getline(v1, line);
    std::getline(v1, line);
    while (std::getline(v1, line))
    {
        line = line.substr(line.find(',') + 1);
        string header = line.substr(0, line.find(','));
        line = line.substr(line.find(',') + 1);
        int lastDelimiter = header.find_last_of('|');
        if (lastDelimiter != -1 && header.substr(lastDelimiter + 1) == "RequiresSNPConfirmation")
            header = header.substr(0, lastDelimiter);
        source_to_header.emplace(line, header);
    }
    v1.close();
    unordered_map<string, list<InfoPipe*>> source_model;
    for (auto iter = source_to_header.begin(); iter != source_to_header.end(); ++iter)
    {
        if (header_model.find(iter->second) != header_model.end()) {
            list<InfoPipe*> temp = header_model.at(iter->second);
            source_model.emplace(iter->first, temp); 
        }
    }
    ifstream v2;
    v2.open("../mapping_files/megaresv2_to_external_header_mappings_v2.00.csv");
    unordered_map<string, string> header2_source;
    std::getline(v2, line);
    std::getline(v2, line);
    while (std::getline(v2, line, ','))
    {
        string source;
        std::getline(v2, source, ',');
        string header2;
        std::getline(v2, header2, '|');
        header2_source.emplace(header2, source);
        std::getline(v2, line);
    }
    v2.close();
    for (auto iter = header2_source.begin(); iter != header2_source.end(); ++iter)
    {
        if (source_model.find(iter->second) != source_model.end()) {
            snpInfoDatabase.emplace(iter->first, source_model.at(iter->second)); 
        }
    }
}
