#pragma once
#include <fstream>
#include <string>
#include <unordered_map>
#include <iostream>
#include "Gene.h"
#include "MmarcModel.h"

using namespace std;

class GeneDatabase{
    private:
        unordered_map<string, Gene*>  genes;
        unordered_map<string, list<MmarcModel*>> snpInfoDatabase;
    public:
        GeneDatabase(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
        void addGene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
        void SNPInfo();
        void print(string fileName);
};

GeneDatabase::GeneDatabase(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence)
{
    SNPInfo();
    addGene(_geneName, _geneType, _geneClass, _geneMechanism, _geneGroup, _geneSequence);
}
void GeneDatabase::addGene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence)
{
    Gene* toAdd = new Gene(_geneName, _geneType, _geneClass, _geneMechanism, _geneGroup, _geneSequence);
    genes.emplace(_geneName, toAdd);

}
void GeneDatabase::SNPInfo()
{
    ifstream snpsearch;
    snpsearch.open("mmarc_snpsearch_metadata2.txt");
    unordered_map<string, list<MmarcModel*>> name_model;
    string line;
    std::getline(snpsearch, line);
    while (std::getline(snpsearch, line))
    {
        string modelName = line.substr(0, line.find(','));
        try{name_model.at(modelName).push_back(new MmarcModel(line));}
        catch (const out_of_range& oor) {name_model.emplace(modelName, (new MmarcModel(line)));}
    }
    snpsearch.close();
    ifstream model_members;
    model_members.open("mmarc_model_members.csv");
    unordered_map<string, string> header_name;
    std::getline(model_members, line);
    while (std::getline(model_members, line))
    {
        string modelName = line.substr(0, line.find(','));
        while(line.find(','))
        {
            line = line.substr(line.find(',') + 1);
            header_name.emplace(line.substr(0, line.find(',')), modelName);
        }

    }
    model_members.close();
    unordered_map<string, list<MmarcModel*>> header_model;
    for (auto iter = header_name.begin(); iter != header_name.end(); ++iter)
    {
        try{header_model.emplace(iter->first, name_model.at(iter->second));}
        catch (const out_of_range& oor) {}
    }
    ifstream v1;
    v1.open("megaresv1_to_external_header_mappings_v1.01.csv");
    unordered_map<string, string> source_to_header;
    std::getline(v1, line);
    std::getline(v1, line);
    while (std::getline(v1, line))
    {
        line = line.substr(line.find(',') + 1);
        string header = line.substr(0, line.find(','));
        while(line.find(' '))
        {
            line = line.substr(line.find(' ') + 1);
            source_to_header.emplace(line.substr(0, line.find(' ')), header);
        }
    }
    model_members.close();
    unordered_map<string, list<MmarcModel*>> source_model;
    for (auto iter = source_to_header.begin(); iter != source_to_header.end(); ++iter)
    {
        try{source_model.emplace(iter->first, header_model.at(iter->second));}
        catch (const out_of_range& oor) {}
    }
    ifstream v2;
    v2.open("megaresv1_to_external_header_mappings_v1.01.csv");
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
    model_members.close();
    for (auto iter = header2_source.begin(); iter != header2_source.end(); ++iter)
    {
        try{snpInfoDatabase.emplace(iter->first, source_model.at(iter->second));}
        catch (const out_of_range& oor) {}
    }
}
void GeneDatabase::print(string fileName)
{
    ofstream output;
    output.open(fileName);
    for (auto iter = genes.begin(); iter != genes.end(); ++iter)
    {
        output << iter->second->getFASTA();
    }
    output.close();
}