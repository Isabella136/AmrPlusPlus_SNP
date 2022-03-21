#pragma once
#include <fstream>
#include <string>
#include <unordered_map>
#include <iostream>
#include "Gene.h"
#include "MmarcModel.h"

using namespace std;

class GeneDatabase {
private:
    unordered_map<string, Gene*>  genes;
    unordered_map<string, list<MmarcModel*>> snpInfoDatabase;
public:
    GeneDatabase();
    void addGene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence);
    void SNPInfo();
    void print(string fileName);
};

GeneDatabase::GeneDatabase()
{
    SNPInfo();
}

void GeneDatabase::addGene(string _geneName, string _geneType, string _geneClass, string _geneMechanism, string _geneGroup, string _geneSequence)
{
    Gene* toAdd = new Gene(_geneName, _geneType, _geneClass, _geneMechanism, _geneGroup, _geneSequence);
    try
    {
        list<MmarcModel*> allSNPs = snpInfoDatabase.at(_geneName);
        for (auto iter = allSNPs.begin(); iter != allSNPs.end(); ++iter)
        {
            toAdd->addSNP((*iter)->condensedSNPinfo().first, (*iter)->condensedSNPinfo().second);
        }
    }
    catch (const out_of_range & oor) {}
    genes.emplace(_geneName, toAdd);
}
void GeneDatabase::SNPInfo()
{
    ifstream snpsearch;
    snpsearch.open("mmarc_snpsearch_metadata2_modified.txt");
    unordered_map<string, list<MmarcModel*>> name_model;
    string line;
    std::getline(snpsearch, line);
    while (std::getline(snpsearch, line))
    {
        string modelName = line.substr(0, line.find(','));
        try { name_model.at(modelName).push_back(new MmarcModel(line)); }
        catch (const out_of_range & oor) { 
			list<MmarcModel*> temp;
			temp.push_back(new MmarcModel(line));
			name_model.emplace(modelName, temp);
			}
    }
    snpsearch.close();
    ifstream model_members;
    model_members.open("mmarc_model_members.csv");
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
    unordered_map<string, list<MmarcModel*>> header_model;
    for (auto iter = header_name.begin(); iter != header_name.end(); ++iter)
    {
        try { 
            header_model.emplace(iter->first, name_model.at(iter->second)); 
            }
        catch (const out_of_range & oor) {
        }
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
        line = line.substr(line.find(',') + 1);
        int lastDelimiter = header.find_last_of('|');
        if (lastDelimiter != -1 && header.substr(lastDelimiter + 1) == "RequiresSNPConfirmation")
            header = header.substr(0, lastDelimiter);
        source_to_header.emplace(line.substr(0, line.find(' ')), header);
        while (line.find(' ') != -1)
        {
            line = line.substr(line.find(' ') + 1);
            source_to_header.emplace(line.substr(0, line.find(' ')), header);
        }
    }
    model_members.close();
    unordered_map<string, list<MmarcModel*>> source_model;
    ofstream test1;
    test1.open("test1.csv");
    for (auto iter = source_to_header.begin(); iter != source_to_header.end(); ++iter)
    {
        test1 << iter->first << ',' << iter->second << '\n';
    }
    test1.close();
    ofstream test2;
    test2.open("test2.csv");
    for (auto iter = header_model.begin(); iter != header_model.end(); ++iter)
    {
        test2 << iter->first << '\n';
    }
    test2.close();
    for (auto iter = source_to_header.begin(); iter != source_to_header.end(); ++iter)
    {
        try {
            list<MmarcModel*> temp = header_model.at(iter->second);
            source_model.emplace(iter->first, temp); 
            }
        catch (const out_of_range & oor) {
        }
    }
    ifstream v2;
    v2.open("megaresv2_to_external_header_mappings_v2.00.csv");
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
        try { snpInfoDatabase.emplace(iter->first, source_model.at(iter->second)); }
        catch (const out_of_range & oor) {}
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