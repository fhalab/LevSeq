#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<filesystem>
#include"io.hpp"
#include <unordered_map>
#include <algorithm>
#include <map>
#include <zlib.h>  




std::vector<FastaEntry> read_fasta(const std::string& file_path) {
    std::vector<FastaEntry> entries;
    std::ifstream infile(file_path);

    if (!infile.is_open()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return entries;
    }

    std::string line;
    FastaEntry entry;

    while (std::getline(infile, line)) {
        // Check if the line is an identifier
        if (!line.empty() && line[0] == '>') {
            // If we're not on the first entry, save the previous one
            if (!entry.identifier.empty() && !entry.sequence.empty()) {
                entries.push_back(entry);
                entry.sequence.clear();
            }
            entry.identifier = line;
        } 
        // Sequence line
        else {
            entry.sequence += line;
        }
    }

    // Add the last entry if it's not empty
    if (!entry.identifier.empty() && !entry.sequence.empty()) {
        entries.push_back(entry);
    }
    
    return entries;
}

std::vector<FastqEntry> read_fastq(const std::string& file_path) {
    std::vector<FastqEntry> entries;
    
    // Check file extension
    bool is_gzipped = (file_path.substr(file_path.find_last_of(".")) == ".gz");

    if(is_gzipped) {
        // Reading gzipped fastq
        gzFile file = gzopen(file_path.c_str(), "rb");
        if(file == nullptr) {
            std::cerr << "Error opening gzipped file: " << file_path << std::endl;
            return entries;
        }

        char buffer[5000];
        FastqEntry entry;
        int line_num = 0;
        while(true) {
            char* res = gzgets(file, buffer, sizeof(buffer) - 1);
            if(res == Z_NULL) break;  // End of file or error

            std::string line(buffer);
            
            // Remove newline character
            if (!line.empty() && line[line.size() - 1] == '\n') {
                line.erase(line.size() - 1);
            }
            
            // Parsing fastq entry lines
            if(line_num % 4 == 0) {
                entry.identifier = line;
            } else if(line_num % 4 == 1) {
                entry.sequence = line;
            } else if(line_num % 4 == 2) {
                entry.plus_line = line;
            } else if(line_num % 4 == 3) {
                entry.quality_scores = line;
                entries.push_back(entry);
            }
            line_num++;
        }

        // std::cout << "Added entry: "
        //   << "\nID: " << entry.identifier
        //   << "\nSeq: " << entry.sequence
        //   << "\n+: " << entry.plus_line
        //   << "\nQual: " << entry.quality_scores
        //   << std::endl;

        gzclose(file);

    } else {
        // Reading non-gzipped fastq
        std::ifstream infile(file_path);
        if (!infile.is_open()) {
            std::cerr << "Error opening file: " << file_path << std::endl;
            return entries;
        }

        FastqEntry entry;
        while (std::getline(infile, entry.identifier)) {
            if (entry.identifier.empty() || entry.identifier[0] != '@') {
                std::cerr << "Malformed fastq entry - expected '@' at start of line" << std::endl;
                return {};
            }
            if (!std::getline(infile, entry.sequence) ||
                !std::getline(infile, entry.plus_line) ||
                !std::getline(infile, entry.quality_scores)) {
                std::cerr << "Malformed fastq entry - expected four lines per entry" << std::endl;
                return {};
            }
            entries.push_back(entry);
        }
    }

    return entries;
}



std::string get_reverse_complement(const std::string& sequence) {
    static const std::unordered_map<char, char> complement {
        {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'},
        {'a', 't'}, {'t', 'a'}, {'c', 'g'}, {'g', 'c'},
        {'N', 'N'}, {'n', 'n'} 
    };

    std::string rev_comp;
    rev_comp.reserve(sequence.size());

    for (auto base : sequence) {
        auto found = complement.find(base);
        if (found != complement.end()) {
            rev_comp.push_back(found->second);
        } else {
            throw std::invalid_argument("Invalid base found: " + std::string(1, base));
        }
    }

    std::reverse(rev_comp.begin(), rev_comp.end());
    return rev_comp;
}

std::string get_reverse_qualities(const std::string& qualities) {
    std::string rev_qualities;
    rev_qualities.reserve(qualities.size());

    for (auto it = qualities.rbegin(); it != qualities.rend(); ++it) {
        rev_qualities.push_back(*it);
    }

    return rev_qualities;
}




// Barcodes

bool is_rbc_barcode(const std::string& s) {
    // Check if it starts with "RB"
    if (s.substr(0, 2) != "RB") return false;
    // Check if the remaining characters are digits
    for (size_t i = 2; i < s.size(); ++i) {
        if (!std::isdigit(s[i])) return false;
    }
    return true;
}

bool is_fbc_barcode(const std::string& s) {
    // Check if it starts with "NB"
    if (s.substr(0, 2) != "NB") return false;
    // Check if the remaining characters are digits
    for (size_t i = 2; i < s.size(); ++i) {
        if (!std::isdigit(s[i])) return false;
    }
    return true;
}


std::pair<std::map<std::string, std::string>, std::map<std::string, std::string>> get_barcodes(const std::string& file_path){
    std::ofstream bar_file("barcodes.txt");
    if (!bar_file.is_open()) {
        std::cerr << "Error opening barcode file for writing\n";
        return {};
    }

    bar_file << "Barcode\tSequence\n";
    std::vector<FastaEntry> barcodes = read_fasta(file_path);
    std::map<std::string, std::string> rbc_map; // RBC
    std::map<std::string, std::string> fbc_map; // FBC

    for (const auto& entry : barcodes) {
        std::string barcode_candidate = entry.identifier.substr(1, 4);

        if (is_rbc_barcode(barcode_candidate) && barcode_candidate.size() == 4) {
            rbc_map.insert({barcode_candidate, entry.sequence});
            bar_file << barcode_candidate << "\t" << entry.sequence << "\n";
        }

        else if (is_fbc_barcode(barcode_candidate) && barcode_candidate.size() == 4) {
            fbc_map.insert({barcode_candidate, entry.sequence});
            bar_file << barcode_candidate << "\t" << entry.sequence << "\n";
        }

        else {
            std::cerr << "Barcode " << barcode_candidate << " is not a valid barcode\n";
        }
    }

    bar_file.close();

    return {fbc_map, rbc_map};
}





