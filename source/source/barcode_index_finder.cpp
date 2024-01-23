#include <iostream>
#include <utility> 
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <filesystem> 
#include "alignments.hpp" // Add alignments
#include "io.hpp" // IO fo such as read fasta files
#include <map>
#include <cctype> 
#include <chrono>
#include <tuple>
#include <sstream>
#include <algorithm>
#include <regex>
using namespace std::chrono;



#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

std::vector<FastaEntry> readFasta(const std::string& filepath) {
    std::vector<FastaEntry> entries;
    std::ifstream file(filepath);
    std::string line, seq;
    FastaEntry entry;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!entry.sequence.empty()) {
                entries.push_back(entry); // Save previous entry
                entry = FastaEntry(); // Reset entry
            }
            entry.identifier = line.substr(1); // Skip '>'
        } else {
            entry.sequence += line;
        }
    }
    if (!entry.sequence.empty()) entries.push_back(entry); // Save last entry

    return entries;
}

std::map<std::string, std::string> readBarcodes(const std::string& filepath) {
    std::map<std::string, std::string> barcodes;
    std::ifstream file(filepath);
    std::string line, name, seq;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!seq.empty()) {
                barcodes[name] = seq; // Save previous barcode
                seq.clear(); // Reset sequence
            }
            name = line.substr(1); // Skip '>'
        } else {
            seq += line;
        }
    }
    if (!seq.empty()) barcodes[name] = seq; // Save last barcode

    return barcodes;
}

void findAndWriteBarcodes(const std::vector<FastaEntry>& sequences, const std::map<std::string, std::string>& barcodes, const std::string& outputFilePath) {
    std::ofstream outfile(outputFilePath);
    for (const auto& entry : sequences) {
        for (const auto& [barcodeName, barcodeSeq] : barcodes) {
            size_t startPos = entry.sequence.find(barcodeSeq);
            if (startPos != std::string::npos) {
                size_t endPos = startPos + barcodeSeq.length() - 1;
                outfile << entry.identifier << "\t" << barcodeName << "\t" << startPos << "\t" << endPos << std::endl;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if(argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file> <barcode_file> <output_file>" << std::endl;
        return 1;
    }

    std::string fastaFilePath = argv[1];
    std::string barcodeFilePath = argv[2];
    std::string outputFilePath = argv[3];

    auto sequences = readFasta(fastaFilePath);
    auto barcodes = readBarcodes(barcodeFilePath);

    findAndWriteBarcodes(sequences, barcodes, outputFilePath);

    return 0;
}
