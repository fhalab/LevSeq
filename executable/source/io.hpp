// io.cpp

#pragma once
#include <string>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <map>


struct FastaEntry {
    std::string identifier;
    std::string sequence;
};

struct FastqEntry {
    std::string identifier;
    std::string sequence;
    std::string plus_line;
    std::string quality_scores;
};

int extractNumericSuffix(const std::string& filename);

bool compareFileNames(const std::filesystem::directory_entry& a, const std::filesystem::directory_entry& b); 
bool isValidBase(char base);
bool isValidEntry(const FastqEntry& entry); 


std::vector<FastaEntry> read_fasta(const std::string& file_path);

std::vector<FastqEntry> read_fastq(const std::string& file_path);

std::string get_reverse_complement(const std::string& sequence);

std::string get_reverse_qualities(const std::string& qualities);

bool is_rbc_barcode(const std::string& s);

bool is_fbc_barcode(const std::string& s);

std::pair<std::map<std::string, std::string>, std::map<std::string, std::string>> get_barcodes(const std::string& file_path);