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


int main(){


    // Get Fastq file path
    std::filesystem::path folder_path = std::filesystem::path("/home/emre/github_repo/MinION/examples/data/tpr_fpr_testing/24_barcodes");

    
    // Check if path exists and is a directory
    if (!std::filesystem::exists(folder_path) || !std::filesystem::is_directory(folder_path)) {
        std::cerr << "Error: Path does not exist or is not a directory.\n";
        return 1; // Return with error code
        }


    std::vector<std::filesystem::directory_entry> files;
    localAlignmentResult score;
    
    std::filesystem::path demultiplex_folder = std::filesystem::path("/home/emre/github_repo/MinION/examples/data") / "Demultiplex_cpp";

    // Loop trough different thresholds
    std::vector<int> thresholds = {10};

    for (int threshold : thresholds) {

        // demultiplex_folder = std::filesystem::path("/home/emre/github_repo/MinION/examples/data") / "Demultiplex_cpp" / std::to_string(threshold);
        demultiplex_folder = std::filesystem::path("/home/emre/github_repo/MinION/examples/data") / "Demultiplex_cpp";
        // Check if Demultiplex folder exists, if not create it
    try {
        if (!std::filesystem::exists(demultiplex_folder)) {
            std::filesystem::create_directory(demultiplex_folder);
        }
    } catch (std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating main results directory: " << e.what() << std::endl;
        // Optionally: return or exit here since directory creation failed
        return 1;
    }


    // Create summary file
    std::filesystem::path summary_file_path = demultiplex_folder / "barcoding_summary_24_barcodes.txt"; // For multithread, create a file for each thread and merge them at the end
    std::filesystem::path barcode_index_file_path = demultiplex_folder / "barcode_index.txt";
    std::ofstream summary_file(summary_file_path.string());
    std::ofstream barcode_index_file(barcode_index_file_path.string());

    summary_file << "ID\tRBC\tRBC_Score\n"; // Header

    barcode_index_file << "ID\tRBC\tRBC_Score\n"; // Header


    // Get mapped barcodes
    auto [fbc_map, rbc_map] = get_barcodes("/home/emre/github_repo/MinION/minION/barcoding/minion_barcodes_sim.fasta"); //TODO, change it into an argument, so the user can use their own barcodes


    // Combine Reverse barcodes and reverse complements
    std::map<std::string, std::string> all_rbc = rbc_map;
   


    // Create a Forward - RBC map
    std::map<std::string, std::string> fbc_rev_map;
    for(const auto& [key, seq] : fbc_map) {
        std::string rev_comp = get_reverse_complement(seq);
        fbc_rev_map[key + "-Rev"] = rev_comp;
        }
    

    std::map<std::string, int> rbc_precalculated_scores;

    // Map for Bases
    std::map<char, int> base_map = {
        {'A', 96},
        {'C', 100},
        {'G', 98},
        {'T', 100}
    };


    for (const auto& [seq_id, sequence] : all_rbc) {
        int all_rbc_score = 0;
        for (char base : sequence) {
            all_rbc_score += base_map[base]; // Add the base score for each nucleotide in the sequence
        }
        rbc_precalculated_scores[seq_id] = all_rbc_score; // Store the total score with the sequence ID as the key
        }


    // Loop through all files in the folder

    for (const auto& file : std::filesystem::directory_iterator(folder_path)) {
        
        std::string file_path = file.path().string();

        if (file.path().extension() == ".fastq" || file.path().extension() == ".gz") {
            
            std::cout << "Processing file: " << file_path << "\n";
            
            std::vector<FastqEntry> entries = read_fastq(file_path);

            for (auto& entry : entries) {

                std::cout << "Processing entry: " << entry.identifier << "\n";

                // Check if entry is empty
                if (entry.identifier.empty() || entry.sequence.empty()) {
                    std::cerr << "Error: Empty entry found in file: " << file_path << std::endl;
                    continue; // Skip to the next entry
                    }

                // Define RBC variables
                std::string best_rbc_name;
                double best_percent_score = 0;

                std::string forward_subseq = entry.sequence;

                for(const auto& [ref_name, ref_seq] : all_rbc) {
                
                    double rbc_sum_score = rbc_precalculated_scores[ref_name];


                    score = perform_alignment_trim(forward_subseq, ref_seq, scoring_matrix2);

                    std::cout << "Score: " << score.score << "\n"; 

                    double percent_score = (double)score.score / rbc_sum_score * 100;


                    // If best score is already 100, skip the loop
                    if (percent_score == 100) {
                        best_rbc_name = ref_name;
                        best_percent_score = percent_score;
                        break;
                        }

                    // Update best score and associated name if current score is better
                    if (percent_score > best_percent_score) {
                        best_rbc_name = ref_name;
                        best_percent_score = percent_score;
                        }

                }

                std::cout << "Best RBC: " << best_rbc_name << "\t" << "Score: " << best_percent_score << "\n";


                // Check if best score is above thershold%, otherwise assign as unclassified
                if (best_percent_score < threshold) {
                    best_rbc_name = "unclassified";
                    } 

                // Output to file
                summary_file << entry.identifier << "\t" <<best_rbc_name << "\t" << best_percent_score << "\n";

            }
        }
        else {
            std::cout << "Skipping file: " << file_path << "\n";
            continue;
        }
    }
    }
}