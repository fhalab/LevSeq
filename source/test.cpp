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
using namespace std::chrono;


// Note: Define `get_barcodes` function as per your requirement.

int main()
{

    std::filesystem::path folder_path = std::filesystem::current_path() / ".." / "data";

    // Front & Rear Window Size for barcode matching
    int front_window_size = 100;
    int rear_window_size = 100;

    // Start Demultiplexing
    std::filesystem::path demultiplex_folder = std::filesystem::current_path()/ ".." / "data" / "Demultiplex_test"; // Create Demultiplex folder - Ideally Experiment folder


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

    // Get mapped barcodes
    auto [fbc_map, rbc_map] = get_barcodes("/home/emre/github_repo/MinION/minION/barcoding/minion_barcodes.fasta"); //TODO, change it into an argument, so the user can use their own barcodes


    // Create RBC reverse map 

    std::map<std::string, std::string> reverse_complements;
    for(const auto& [key, seq] : rbc_map) {
        std::string rev_comp = get_reverse_complement(seq);
        reverse_complements[key + "-Rev"] = rev_comp;
        }

    // Combine Reverse barcodes and reverse complements
    std::map<std::string, std::string> all_rbc = rbc_map;
    all_rbc.insert(reverse_complements.begin(), reverse_complements.end());


    // Create a Forward - RBC map
    std::map<std::string, std::string> fbc_rev_map;
    for(const auto& [key, seq] : fbc_map) {
        std::string rev_comp = get_reverse_complement(seq);
        fbc_rev_map[key + "-Rev"] = rev_comp;
        }
    

    std::map<std::string, int> rbc_precalculated_scores;
    std::map<std::string, int> fbc_precalculated_scores;

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

    for (const auto& [seq_id, sequence] : fbc_map) {
        int fbc_score = 0;
        for (char base : sequence) {
            fbc_score += base_map[base]; // Add the base score for each nucleotide in the sequence
        }
        fbc_precalculated_scores[seq_id] = fbc_score; // Store the total score with the sequence ID as the key
        }


   // Iterate through each file in the directory - Use multithreading

    // Define RBC and FBC variables in an outer scope
    localAlignmentResult score;
    localAlignmentResult fbc_score;
    int sequence_length;

    int counter = 0;
    for (const auto& file : std::filesystem::directory_iterator(folder_path)) {
        
        std::string file_path = file.path().string();
        
        if (file.path().extension() == ".fastq" || file.path().extension() == ".gz") {
            
            std::cout << "Processing file: " << file_path << "\n";
            
            std::vector<FastqEntry> entries = read_fastq(file_path);

            // Define RBC variables
            std::string best_rbc_name;
            double best_percent_score = 0;
        
            // Define FBC variables
            std::string best_fbc_name;
            double best_fbc_percent_score = 0;

            for (auto& entry : entries) {

                if (counter > 5)
                {
                    break;
                }
                

                // Check if entry is empty
                if (entry.identifier.empty() || entry.sequence.empty()) {
                    std::cerr << "Error: Empty entry found in file: " << file_path << std::endl;
                    continue; // Skip to the next entry
                    }

                if(entry.sequence.size() < front_window_size) {
                    std::cerr << "Sequence less than 100 bases" << std::endl;
                    continue; // Skip to the next entry
                    }
                // Check length of Sequence with min and max length
                else if(entry.sequence.size() > 1000) {
                    std::cerr << "Sequence more than 1000 bases" << std::endl;
                    continue; // Skip to the next entry
                    }
                else if(entry.sequence.size() < 750) {
                    std::cerr << "Sequence less than 750 bases" << std::endl;
                    continue; // Skip to the next entry
                }


                // Get Front and Rear Subsequences
                std::string forward_subseq = entry.sequence.substr(0, front_window_size);
                std::string rear_subseq = entry.sequence.substr(entry.sequence.size() - rear_window_size);


                // Reverse Barcode Demultiplexing
                for(const auto& [ref_name, ref_seq] : all_rbc) {
                    
                    double rbc_sum_score = rbc_precalculated_scores[ref_name];

                    std::string sequence_to_align;
                    if(ref_name.find("-Rev") != std::string::npos) {
                        sequence_to_align = rear_subseq;
                    } else {
                        sequence_to_align = forward_subseq;
                    }

                    //int score = perform_alignment(sequence_to_align, ref_seq); // Get
                    score = perform_alignment_trim(sequence_to_align, ref_seq, scoring_matrix2); // Get

                    
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




                if(best_rbc_name.find("-Rev") == std::string::npos) { // String does not contain -Rev

                    // Cut at end position since the sequence is forward
                    // trim from score.end to end of sequence
                    entry.sequence = entry.sequence.substr(score.end_pos);
                    entry.quality_scores = entry.quality_scores.substr(score.end_pos); // Trimming new sequence size is sequence.size() - end_pos

                    entry.sequence = get_reverse_complement(entry.sequence);
                    entry.quality_scores = get_reverse_qualities(entry.quality_scores);
                    // Update Forward and Rear subsequences
                    forward_subseq = entry.sequence.substr(0, front_window_size);
                    }

                else {
                    best_rbc_name = best_rbc_name.substr(0, best_rbc_name.size() - 4); // Remove -Rev from the end
                    // Cut at start position since the sequence is reverse
                    score.start_pos = entry.sequence.size() - rear_window_size + score.start_pos;
                    score.end_pos = entry.sequence.size() - rear_window_size + score.end_pos;
                    entry.sequence = entry.sequence.substr(0, score.start_pos);
                    entry.quality_scores = entry.quality_scores.substr(0, score.start_pos); // Trimming new sequence size is start_pos
                    }

                // Check if best score is above thershold%, otherwise assign as unclassified
                if (best_percent_score < 70) {
                    best_rbc_name = "unclassified";
                    }


                if (best_rbc_name != "unclassified") {
                    // Forward Barcode Demultiplexing
                    for (const auto& [fbc_ref_name, fbc_ref_seq] : fbc_map) {
                        double fbc_sum_score = fbc_precalculated_scores[fbc_ref_name];
                        
                        //int fbc_score = perform_alignment(forward_subseq, fbc_ref_seq);
                        fbc_score = perform_alignment_trim(forward_subseq, fbc_ref_seq, scoring_matrix2); // Get

                        
                        double fbc_percent_score = (double)fbc_score.score / fbc_sum_score * 100;
                        
                        if (fbc_percent_score == 100) {
                            best_fbc_name = fbc_ref_name;
                            best_fbc_percent_score = fbc_percent_score;
                            break;
                            }
                        
                        // Update best score and associated name if the current score is better
                        if (fbc_percent_score > best_fbc_percent_score) {
                            best_fbc_name = fbc_ref_name;
                            best_fbc_percent_score = fbc_percent_score;
                            }
                        }

                        if (best_fbc_percent_score < 70) {
                            best_fbc_name = "unclassified";
                            }

                        else{
                            // Trim
                            entry.sequence = entry.sequence.substr(fbc_score.end_pos);
                            entry.quality_scores = entry.quality_scores.substr(fbc_score.end_pos); // Trimming new sequence size is sequence.size() - end_pos
                            }
                        }

                        sequence_length = entry.sequence.size();

                    //std::cout << "Best RBC: " << best_rbc_name << " Score: " << score.score << " Start Position RBC: " << score.start_pos << " End Position RBC: " << score.end_pos << " Seq Length: " <<  sequence_length <<std::endl;
                    //std::cout << "Best FBC: " << best_fbc_name << " Score: " << fbc_score.score << " Start Position FBC: " << fbc_score.start_pos << " End Position FBC: " << fbc_score.end_pos << std::endl;
                    //std::cout << "Sequence: " << entry.sequence << std::endl;

                    entry.sequence = entry.sequence.substr(0, 200);
                    std::cout << "Sequence: " << entry.sequence << std::endl;
                    // Reverse
                    entry.sequence = get_reverse_complement(entry.sequence);
                    std::cout << "Reverse Sequence: " << entry.sequence << std::endl;

                    counter++;
                    }

                }
            
            break;

        }

    return 0;
    }