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
#include <vector>

using namespace std::chrono;
namespace fs = std::filesystem;


/**
 * @file main.cpp
 * @brief Main entry point for the application that processes fastq files.
 *
 * This program is designed to handle fastq files for genomic sequencing data. It takes several
 * command line arguments to specify the details of file processing, including the locations of 
 * the files, demultiplexer folder path, barcode information, and parameters for window sizes.
 * 
 * @note The program requires specific command-line arguments in a particular order. Future updates 
 *       should include key-value pair parsing for more flexible argument input.
 *
 * Usage:
 *     ./program_name <folder_name> <demultiplexer_folder_path> <barcode_fasta> <front_window_size> <rear_window_size>
 *
 * Where:
 *     - <folder_name>: String. The name of the folder where the basecalled fastq files are located.
 *     - <demultiplexer_folder_path>: String. The path to the folder containing demultiplexer tools.
 *     - <barcode_fasta>: String. File path to the barcode fasta file.
 *     - <front_window_size>: Integer. The size of the front window for processing the fastq files.
 *     - <rear_window_size>: Integer. The size of the rear window for processing the fastq files.
 *
 * @param argc The count of command-line arguments passed to the program.
 * @param argv An array of character pointers listing all the arguments.
 *
 * @return int Returns 0 on successful execution, and 1 if the program encounters an error, 
 *             specifically if the number of arguments provided is incorrect.
 */


void printHelp() {
    std::cout << "Usage: ./myprogram [OPTIONS]\n\n"
              << "Options:\n"
              << "  --folder_name, -f                 Required, Path to the fastq/fastq.gz files, string\n"
              << "  --demultiplexer_folder_path, -d   Required, Path where the demultiplexed sequences should be saved, string\n"
              << "  --barcode_fasta, -b               Required, Path to barcode fasta file, string\n"
              << "  --front_window_size, -w           Required, integer\n"
              << "  --rear_window_size, -r            Required, integer\n"
              << "  --min_length                    Optional, integer, default: 0\n"
              << "  --max_length                    Optional, integer, default: 10000\n"
              << "  --match_score                 Optional, integer, default: 1\n"
              << "  --mismatch_score              Optional, integer, default: -1\n"
              << "  -h, --help                   Show this help message and exit\n";
}


std::vector<std::string> checkRequiredArgs(const std::map<std::string, std::string>& args, const std::vector<std::string>& requiredArgs) {
    std::vector<std::string> missingArgs;
    for (const auto& arg : requiredArgs) {
        if (args.find(arg) == args.end()) {
            missingArgs.push_back(arg);
        }
    }
    return missingArgs;
}

int main(int argc, char* argv[]) {

    // Helper function output
    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
    printHelp();
    return 0;
    }

    std::map<std::string, std::string> args;
    std::map<std::string, std::string> keyMap = {
        {"-f", "--folder_name"},
        {"-d", "--demultiplexer_folder_path"},
        {"-b", "--barcode_fasta"},
        {"-w", "--front_window_size"},
        {"-r", "--rear_window_size"},
        {"-m", "--min_length"},
        {"-x", "--max_length"},
        {"-s", "--match_score"},
        {"-t", "--mismatch_score"}
    };

    for (int i = 1; i < argc; i += 2) {
        std::string key(argv[i]);

        if (key[0] == '-' && i + 1 < argc) { 
            if (keyMap.find(key) != keyMap.end()) {
                // Store both short and long forms in the map
                args[key] = argv[i + 1];
                args[keyMap[key]] = argv[i + 1];
                std::cout << "Processed argument: " << key << " with value: " << argv[i + 1] << std::endl;

            } else {
                std::cerr << "Invalid argument: " << key << std::endl;
                return 1;
            }
        } else {
            std::cerr << "Invalid argument format: " << key << std::endl;
            return 1;
        }
    }

    std::vector<std::string> requiredArgs = {"-f", "--folder_name", "-d", "--demultiplexer_folder_path", "-b", "--barcode_fasta", "-w", "--front_window_size", "-r", "--rear_window_size"};
    std::vector<std::string> missingArgs = checkRequiredArgs(args, requiredArgs);

    if (!missingArgs.empty()) {
        std::cerr << "Error: The following arguments were not provided:" << std::endl;
        for (const auto& arg : missingArgs) {
            std::cerr << arg.substr(1) << std::endl;
        }
        return 1;
    }

    std::string folderName = args["--folder_name"];
    std::string demultiplexerFolderPath = args["--demultiplexer_folder_path"];
    std::string barcodeFasta = args["--barcode_fasta"];
    int frontWindowSize;
    int rearWindowSize;
    try {
        frontWindowSize = std::stoi(args["--front_window_size"]);
        rearWindowSize = std::stoi(args["--rear_window_size"]);
    } catch (const std::exception& e) {
        std::cerr << "Error: Invalid window size provided." << std::endl;
        return 1;
    }

    // Min and Max Length
    int min_length = 0;
    int max_length = 10000;

    // Check if arguments are provided
    if (args.count("--min_length") > 0) {
        try {
            min_length = std::stoi(args["--min_length"]);
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid min_length provided." << std::endl;
            return 1;
        }
    }
    if (args.count("--max_length") > 0) {
        try {
            max_length = std::stoi(args["--max_length"]);
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid max_length provided." << std::endl;
            return 1;
        }
    }

    // Check if max_length is smaller than min_length
    if (max_length < min_length) {
        std::cerr << "Error: max_length cannot be smaller than min_length." << std::endl;
        return 1;
    }

    // Folder Paths
    fs::path folderPath = fs::path(folderName);
    fs::path demultiplex_folder = fs::path(demultiplexerFolderPath);
    fs::path barcode_fasta_path = fs::path(barcodeFasta);


    // Check if folder exists
    if (!fs::exists(folderPath) || !fs::is_directory(folderPath)) {
    fs::create_directory(folderPath);
    }
    if (!fs::exists(demultiplex_folder) || !fs::is_directory(demultiplex_folder)) {
        fs::create_directory(demultiplex_folder);
    }

    // Check if barcode_fasta_path exists
    if (!fs::exists(barcode_fasta_path)) {
        std::cerr << "Error: Barcode fasta file path does not exist." << std::endl;
        return 1;
    }
    

    // Sort fastq files, First sequences sequences are usually better
    std::vector<fs::directory_entry> files;
    for (const auto& file : fs::directory_iterator(folderPath)) {
        files.push_back(file);
    }

    try {
        std::sort(files.begin(), files.end(), compareFileNames);
    } catch (const std::exception& e) {
        std::cerr << "Warning: An error occurred while sorting the files. Skipping sorting.\n";
    }


    // Create summary file
    fs::path summary_file_path = demultiplex_folder / "barcoding_summary.txt"; // TDDO: For multithread, create a file for each thread and merge them at the end
    std::ofstream summary_file(summary_file_path.string());
    summary_file << "RBC\tRBC_Score\tFBC\tFBC_Score\n"; // TODO: Automatically detect if the file exists and append to it instead of overwriting it


    // Barcodes
    auto [fbc_map, rbc_map] = get_barcodes(barcodeFasta);

    std::map<std::string, std::string> all_rbc = rbc_map;
    for(const auto& [key, seq] : rbc_map) {
        all_rbc[key + "-Rev"] = get_reverse_complement(seq);
    }

    // Create a Forward - RBC map
    std::map<std::string, std::string> fbc_rev_map;
    for(const auto& [key, seq] : fbc_map) {
        fbc_rev_map[key + "-Rev"] = get_reverse_complement(seq);
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

    auto calculate_scores = [&](const auto& map, auto& scores) {
        for (const auto& [seq_id, sequence] : map) {
            int score = 0;
            for (char base : sequence) {
                score += base_map[base]; // Add the base score for each nucleotide in the sequence
            }
            scores[seq_id] = score; // Store the total score with the sequence ID as the key
        }
    };

    calculate_scores(all_rbc, rbc_precalculated_scores);
    calculate_scores(fbc_map, fbc_precalculated_scores);

    struct BarcodeData {
        std::ofstream file_stream;
        int entry_count = 0; // Initialize entry_count to 0
        std::string subdirectory_name; // Store the subdirectory name

        BarcodeData() = default;

        BarcodeData(const std::string& file_path, const std::string& subdirectory)
            : file_stream(file_path, std::ios::app), subdirectory_name(subdirectory) {}
    };

    // Outputmap for fastq files
    std::map<std::string, std::map<std::string, BarcodeData>> barcode_to_ofstream;
    std::map<std::string, int> barcode_file_counters;



    // Alignment Score Values
    localAlignmentResult score;
    localAlignmentResult fbc_score;
    int file_count = 0;


    // Get number of files
    int n_files = std::distance(std::filesystem::directory_iterator(folderPath), std::filesystem::directory_iterator{});
    std::cout << "Number of files: " << n_files << "\n";

    // Process files
    int processed_files = 0;


    // Main Loop, that processes each file
    for (const auto& file : std::filesystem::directory_iterator(folderPath)) {

        std::string file_path = file.path().string();
        
        if (file.path().extension() == ".fastq" || file.path().extension() == ".gz") {
                        
            std::vector<FastqEntry> entries = read_fastq(file_path);


            for (auto& entry : entries) {

                // Check if entry is empty
                if (entry.identifier.empty() || entry.sequence.empty()) {
                    std::cerr << "Error: Empty entry found in file: " << file_path << std::endl;
                    continue; // Skip to the next entry
                    }

                if (!isValidEntry(entry)) {
                    //std::cerr << "Error: Invalid base found in sequence in file: " << file_path << std::endl; # TODO: Count reads with invalid bases
                    continue; 
                    }

                // Define RBC variables
                std::string best_rbc_name;
                double best_rbc_percent_score = 0;
            
                // Define FBC variables
                std::string best_fbc_name;
                double best_fbc_percent_score = 0;

                if(entry.sequence.size() < frontWindowSize || entry.sequence.size() < rearWindowSize ||
                entry.sequence.size() > max_length || 
                entry.sequence.size() < min_length ||
                (entry.quality_scores.size() != entry.sequence.size())){
                continue;
                }


                // Get Front and Rear Subsequences
                std::string forward_subseq = entry.sequence.substr(0, frontWindowSize);
                std::string rear_subseq = entry.sequence.substr(entry.sequence.size() - rearWindowSize);

                // Check if sizes of front and rear subsequences are correct
                if (forward_subseq.size() != frontWindowSize || rear_subseq.size() != rearWindowSize) {
                    std::cerr << "Error: Subsequence sizes are incorrect." << std::endl;
                    continue;
                }



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
                    score = perform_alignment_trim(sequence_to_align, ref_seq, scoring_matrix2); 


                    double percent_score = (double)score.score / rbc_sum_score * 100;

                    // If best score is already 100, skip the loop
                    if (percent_score == 100) {
                        best_rbc_name = ref_name;
                        best_rbc_percent_score = percent_score;
                        break;
                        }

                    // Update best score and associated name if current score is better
                    if (percent_score > best_rbc_percent_score) {
                        best_rbc_name = ref_name;
                        best_rbc_percent_score = percent_score;
                        }
                    }

                int trim_front = 0;
                int trim_rear = 0;


                if(best_rbc_name.find("-Rev") == std::string::npos) {
                    entry.sequence = entry.sequence.substr(score.end_pos + trim_rear);
                    entry.quality_scores = entry.quality_scores.substr(score.end_pos + trim_rear); // Trimming new sequence size is sequence.size() - end_pos
                    entry.sequence = get_reverse_complement(entry.sequence);
                    entry.quality_scores = get_reverse_qualities(entry.quality_scores);
                    // Update Forward and Rear subsequences
                    forward_subseq = entry.sequence.substr(0, frontWindowSize);
                    }

                else {
                    best_rbc_name = best_rbc_name.substr(0, best_rbc_name.size() - 4); // Remove -Rev from the end
                    score.start_pos = entry.sequence.size() - rearWindowSize + score.start_pos;
                    score.end_pos = entry.sequence.size() - rearWindowSize + score.end_pos;
                    entry.sequence = entry.sequence.substr(0, score.start_pos - trim_rear);
                    entry.quality_scores = entry.quality_scores.substr(0, score.start_pos - trim_rear); // Trimming new sequence size is start_pos
                    }

                // Check if best score is above thershold%, otherwise assign as unclassified
                if (best_rbc_percent_score < 70) {
                    best_rbc_name = "unclassified";
                    }

                if (best_rbc_name != "unclassified") {
                    // Forward Barcode Demultiplexing
                    for (const auto& [fbc_ref_name, fbc_ref_seq] : fbc_map) {
                        double fbc_sum_score = fbc_precalculated_scores[fbc_ref_name];
                        
                        //int fbc_score = perform_alignment(forward_subseq, fbc_ref_seq);
                        fbc_score = perform_alignment_trim(forward_subseq, fbc_ref_seq, scoring_matrix2);


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
                        
                        if (best_fbc_name != "unclassified") {
                            entry.sequence = entry.sequence.substr(fbc_score.end_pos + trim_front);
                            entry.quality_scores = entry.quality_scores.substr(fbc_score.end_pos + trim_front);
                            }


                        }
                    
                    fs::path barcode_subdirectory = best_rbc_name + "/" + best_fbc_name;
                    fs::path barcode_folder = demultiplexerFolderPath / barcode_subdirectory;

                    if (!fs::exists(barcode_folder)) {
                            fs::create_directories(barcode_folder);
                        }


                    // Check if the current entry count for the specific barcode exceeds the limit (e.g., 4000)
                    if (barcode_to_ofstream[best_rbc_name][best_fbc_name].entry_count >= 4000) {
                        barcode_to_ofstream[best_rbc_name][best_fbc_name].entry_count = 0; // Reset entry count
                        barcode_to_ofstream[best_rbc_name][best_fbc_name].file_stream.close(); // Close the current file
                        barcode_file_counters[best_fbc_name]++; // Increment the file counter
                    }

                    // Create and open the ofstream for the current entry
                    if (!barcode_to_ofstream[best_rbc_name][best_fbc_name].file_stream.is_open()) {
                        std::stringstream filename;
                        filename << "demultiplexed_" << best_rbc_name << "_" << best_fbc_name << "_"
                                << std::setw(3) << std::setfill('0') << barcode_file_counters[best_fbc_name] << ".fastq";
                        std::filesystem::path barcode_file_path = barcode_folder / filename.str();

                        barcode_to_ofstream[best_rbc_name][best_fbc_name].file_stream.open(barcode_file_path.string(), std::ios::app);
                    }

                    // Increment the entry count for the current barcode
                    barcode_to_ofstream[best_rbc_name][best_fbc_name].entry_count++;

                    // Write fastq entry to the corresponding file
                    barcode_to_ofstream[best_rbc_name][best_fbc_name].file_stream << entry.identifier << "\n"
                                                                                << entry.sequence << "\n"
                                                                                << "+" << "\n"
                                                                                << entry.quality_scores << "\n";

                    // Output to file
                    summary_file << best_rbc_name << "\t" << best_rbc_percent_score << "\t" << best_fbc_name << "\t" << best_fbc_percent_score << "\n";

            }

        ++processed_files;
        std::cout << "\rProcessing files: [";
        int progress = (processed_files * 50 / n_files);
        for (int i = 0; i < 50; ++i) {
            if (i < progress) {
                std::cout << '#';
            } else {
                std::cout << ' ';
            }
        }
        std::cout << "] " << processed_files * 100 / n_files << "%" << std::flush;
    }

    }
    std::cout << "\n";
    return 0;
}

