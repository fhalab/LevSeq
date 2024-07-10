
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
using namespace std::chrono;



int main(int argc, char* argv[]) {


    /* Argument Parsing 
    TODO: Add key-value pair parsing so the arguments are clear and the order does not matter
    Args:
        - Folder name, where the basecalled fastq files are located
        - Front window size
        - Rear window size

    */


    if(argc != 4) {
        std::cerr << "Usage: " << argv[0] 
                  << "Please provide all arguments required <folder_name> <front_window_size> <rear_window_size>" << std::endl;
        return 1;
        }


    // Define Folder Path - Currently only name with given path TODO: Change it to default path where you start the program
    std::string folder_name(argv[1]);
    //std::filesystem::path folder_path = std::filesystem::current_path() / ".." / folder_name;
    std::filesystem::path folder_path = std::filesystem::path("/home/emre/minION_results/MinION_RBC_0902723_sup") / folder_name;


    // Check if path exists and is a directory
    if (!std::filesystem::exists(folder_path) || !std::filesystem::is_directory(folder_path)) {
        std::cerr << "Error: Path does not exist or is not a directory.\n";
        return 1; // Return with error code
        }

    // Front & Rear Window Size for barcode matching
    int front_window_size;
    int rear_window_size;

    try {
        front_window_size = std::stoi(argv[2]);
        rear_window_size = std::stoi(argv[3]);
        } catch (std::invalid_argument& e) {
        std::cerr << "Error: Invalid window size argument(s). Must be integer(s)." << std::endl;
        return 1; // return an error code, for example 1
        } catch (std::out_of_range& e) {
        std::cerr << "Error: Window size argument(s) out of range." << std::endl;
        return 1; // return an error code, for example 1
        }

    // TODO - Include Open Gap and Extension Gap Scores

    // Start Demultiplexing
    //std::filesystem::path demultiplex_folder = std::filesystem::current_path()/ ".." / "data" / "Demultiplex"; // Create Demultiplex folder - Ideally Experiment folder

    std::filesystem::path demultiplex_folder = std::filesystem::path("/home/emre/minION_results/MinION_RBC_0902723_sup") / "Demultiplex_cpp_70";


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
    std::filesystem::path summary_file_path = demultiplex_folder / "barcode_summary.txt"; // For multithread, create a file for each thread and merge them at the end
    std::ofstream summary_file(summary_file_path.string());

    summary_file << "RBC\tRBC_Score\tFBC\tFBC_Score\n"; // Header

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


    struct BarcodeData {
        std::ofstream file_stream;
        int entry_count = 0;

        BarcodeData() = default;

        BarcodeData(std::string file_path) : file_stream(file_path, std::ios::app) {}
    };

    // Outputmap for fastq files
    std::map<std::string, BarcodeData> barcode_to_ofstream;

    std::map<std::string, int> barcode_file_counters;




    auto start = high_resolution_clock::now();



    // Iterate through each file in the directory - Use multithreading

    for (const auto& file : std::filesystem::directory_iterator(folder_path)) {
        
        std::string file_path = file.path().string();
        
        if (file.path().extension() == ".fastq" || file.path().extension() == ".gz") {
            
            std::cout << "Processing file: " << file_path << "\n";
            
            std::vector<FastqEntry> entries = read_fastq(file_path);


                for (auto& entry : entries) {

                    // Check if entry is empty
                    if (entry.identifier.empty() || entry.sequence.empty()) {
                        std::cerr << "Error: Empty entry found in file: " << file_path << std::endl;
                        continue; // Skip to the next entry
                        }

                    // Define RBC variables
                    std::string best_rbc_name;
                    double best_percent_score = 0;
                
                    // Define FBC variables
                    std::string best_fbc_name;
                    double best_fbc_percent_score = 0;

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

                        int score = perform_alignment(sequence_to_align, ref_seq);

                        double percent_score = (double)score / rbc_sum_score * 100;

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

                    if(best_rbc_name.find("-Rev") == std::string::npos) {
                    
                        entry.sequence = get_reverse_complement(entry.sequence);
                        entry.quality_scores = get_reverse_qualities(entry.quality_scores);
                        // Update Forward and Rear subsequences
                        forward_subseq = entry.sequence.substr(0, front_window_size);
                        }

                    else {
                        best_rbc_name = best_rbc_name.substr(0, best_rbc_name.size() - 4); // Remove -Rev from the end
                        }

                    // Check if best score is above thershold%, otherwise assign as unclassified
                    if (best_percent_score < 70) {
                        best_rbc_name = "unclassified";
                        }

                    std::filesystem::path barcode_folder = demultiplex_folder / "unclassified";

                    if (best_rbc_name != "unclassified"){
                        // Forward Barcode Demultipleixng

                        for(const auto& [fbc_ref_name, fbc_ref_seq] : fbc_map) {
                            
                            double fbc_sum_score = fbc_precalculated_scores[fbc_ref_name];

                            int fbc_score = perform_alignment(forward_subseq, fbc_ref_seq);

                            double fbc_percent_score = (double)fbc_score / fbc_sum_score * 100;

                            if (fbc_percent_score == 100) {
                                best_fbc_name = fbc_ref_name;
                                best_fbc_percent_score = fbc_percent_score;
                                break;
                                }

                            // Update best score and associated name if current score is better
                            if (fbc_percent_score > best_fbc_percent_score) {
                                best_fbc_name = fbc_ref_name;
                                best_fbc_percent_score = fbc_percent_score;
                                }
                            }
                        
                        // Check if best score is above 50%, otherwise assign as unclassified
                        if (best_fbc_percent_score < 70) {
                            best_fbc_name = "unclassified";
                            barcode_folder = demultiplex_folder / best_rbc_name / best_fbc_name;
                            }
                        
                        else {
                            barcode_folder = demultiplex_folder / best_rbc_name / best_fbc_name;
                            }

                        }
                    
                        if (!std::filesystem::exists(barcode_folder)) {
                            std::filesystem::create_directories(barcode_folder);
                        }
                    // Check if ofstream for barcode exists, if not create and open it
                    if(barcode_to_ofstream.find(best_fbc_name) == barcode_to_ofstream.end()) {
                        barcode_file_counters[best_fbc_name] = 1;
                        }

                    // Check if the current ofstream has written 1000 entries, if so, close and open new file
                    if(barcode_to_ofstream[best_fbc_name].entry_count >= 4000) {
                        barcode_to_ofstream[best_fbc_name].file_stream.close(); // Close the current file.
                        barcode_to_ofstream[best_fbc_name].entry_count = 0; // Reset the count.
                        barcode_file_counters[best_fbc_name]++; // Increment the file id counter.
                        }

                    // If ofstream is not open (either first time or after closing), open new file
                    if(!barcode_to_ofstream[best_fbc_name].file_stream.is_open()) {
                        std::stringstream filename;
                        filename << "demultiplexed_" 
                                << best_rbc_name << "_" << best_fbc_name << "_"
                                << std::setw(3) << std::setfill('0') 
                                << barcode_file_counters[best_fbc_name]
                                << ".fastq";
                        std::filesystem::path barcode_file_path = barcode_folder / filename.str();
                        barcode_to_ofstream[best_fbc_name].file_stream.open(barcode_file_path.string(), std::ios::app);
                    }

                    // Write fastq entry to corresponding file.
                    barcode_to_ofstream[best_fbc_name].file_stream  << entry.identifier << "\n"
                                                                    << entry.sequence << "\n"
                                                                    << "+" << "\n"
                                                                    << entry.quality_scores << "\n";

                    // Increment entry count.
                    barcode_to_ofstream[best_fbc_name].entry_count++;

                    // Output to file
                    summary_file << best_rbc_name << "\t" << best_percent_score  << "\t" << best_fbc_name << "\t" << best_fbc_percent_score << "\n";

                
                    // End of Entry Loop
                    }
                }
        }

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    std::cout << "Time taken by function: " << duration.count() << " milli seconds" << std::endl;

    return 0;
}