#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/configuration/align_config_parallel.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/convert.hpp>
#include "alignments.hpp"
#include <map>
#include <string>
#include <thread>
#include <execution>
#include <mutex>
#include <future>
#include <cmath>



using namespace seqan3;


int perform_alignment(std::string const & seq1_str, std::string const & seq2_str)
{   

    auto seq1 = seq1_str | std::views::transform([](char c) { return seqan3::dna15{}.assign_char(c); });
    auto seq2 = seq2_str | std::views::transform([](char c) { return seqan3::dna15{}.assign_char(c); });


    
    
    // Define a custom 15x15 scoring matrix for dna15.
    // TODO: Make a cofiguration function 
    constexpr nucleotide_scoring_scheme<int>::matrix_type scoring_matrix
    {{
        //A   _   C   G  U  AAN  ACN  AGN  ATN  UN  CNN  GNN  TNN  UNN  N
        {{ 96, 0, -316, 0, -192, 0, 0, 0, 0, 0, 0, -396, 0, 0, 0 }}, // A_ind = 0, C_ind = 2, G_ind = 4, T_ind = 12
        {{0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{-316, 0,  100, 0, -352, 0, 0, 0, 0, 0, 0, -295, 0, 0, 0 }},
        {{0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{-192, 0, -352, 0, 98, 0, 0, 0, 0, 0, 0, -329, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{-396, 0, -295, 0, -329, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }}
    }};
    
    nucleotide_scoring_scheme<int> my_scoring_scheme{scoring_matrix};

    // Configure the output:
    auto output_config =    seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_begin_position{} |
                            seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_alignment{};

    // Apply the local alignment configuration.
    auto min_cfg =          seqan3::align_cfg::method_local{} | 
                            seqan3::align_cfg::scoring_scheme{my_scoring_scheme} | output_config |
                            seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-80}, seqan3::align_cfg::extension_score{-40}}; // Set gap open and extension scores
                            


    // Score
    int score = 0;
    const double BARCODE_LENGTH = 24.0;
    try
    {
        for (auto res : seqan3::align_pairwise(std::tie(seq1, seq2), min_cfg))
        {
            
            // seqan3::debug_stream << "Score: " << res.score() << '\n' << "Alignment: \n" << res.alignment() << '\n';
            // seqan3::debug_stream << "Begin: (" << res.sequence1_begin_position() << "," << res.sequence1_end_position()<< ")\n";
            // Get edit distance
            double edit_distance = 0;
            edit_distance = std::abs(static_cast<double>(res.sequence1_end_position() - res.sequence1_begin_position()) - BARCODE_LENGTH);
            //seqan3::debug_stream << "Edit Distance: " << edit_distance << '\n';
            if (edit_distance < 5){
                score = res.score();
            }
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << "An error occurred: " << e.what() << std::endl;
    }

    return score;
}



localAlignmentResult perform_alignment_trim(std::string const & seq1_str, std::string const & seq2_str, const nucleotide_scoring_scheme<int>::matrix_type& scoring_matrix2)
{   

    auto seq1 = seq1_str | std::views::transform([](char c) { return seqan3::dna15{}.assign_char(c); });
    auto seq2 = seq2_str | std::views::transform([](char c) { return seqan3::dna15{}.assign_char(c); });


    nucleotide_scoring_scheme<int> my_scoring_scheme{scoring_matrix2};

    // Configure the output:
    auto output_config =    seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_begin_position{} |
                            seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_alignment{};

    // Apply the local alignment configuration.
    auto min_cfg =          seqan3::align_cfg::method_local{} | 
                            seqan3::align_cfg::scoring_scheme{my_scoring_scheme} | output_config |
                            seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-80}, seqan3::align_cfg::extension_score{-40}}; // Set gap open and extension scores
                            



    // Variables to store the result
    int score = 0;
    int start_pos = 0;
    int end_pos = 0;
    const double BARCODE_LENGTH = 24.0;

    try
    {
        for (auto res : seqan3::align_pairwise(std::tie(seq1, seq2), min_cfg))
        {
            // Get edit distance
            double edit_distance = std::abs(static_cast<double>(res.sequence1_end_position() - res.sequence1_begin_position()) - BARCODE_LENGTH);
            if (edit_distance < 24){
                score = res.score();
                start_pos = res.sequence1_begin_position();
                end_pos = res.sequence1_end_position();
            }
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << "An error occurred: " << e.what() << std::endl;
    }

    return localAlignmentResult{score, start_pos, end_pos};
}








int perform_semiglobal(std::string const & seq1_str, std::string const & seq2_str)
{   

    auto seq1 = seq1_str | std::views::transform([](char c) { return seqan3::dna15{}.assign_char(c); });
    auto seq2 = seq2_str | std::views::transform([](char c) { return seqan3::dna15{}.assign_char(c); });

    return 0;


}


std::vector<sequence_pair_t> build_sequence_pairs(seqan3::dna15_vector const& common_seq, std::map<std::string, std::string> const& barcodes)
{
    std::vector<sequence_pair_t> sequences;
    sequences.reserve(barcodes.size());

    for (auto const& [id, seq_str] : barcodes)
    {
        seqan3::dna15_vector barcode_seq; 
        barcode_seq.reserve(seq_str.size());
        
        for (char c : seq_str) 
        {
            barcode_seq.push_back(seqan3::dna15{}.assign_char(c));
        }

        sequences.emplace_back(id, std::make_pair(common_seq, barcode_seq));
    }

    return sequences;
}



AlignmentResult perform_multi_alignment(std::vector<sequence_pair_t> sequences)
{   

    // Define a custom 15x15 scoring matrix for dna15.
    // TODO: Make a cofiguration function 
    constexpr nucleotide_scoring_scheme<int>::matrix_type scoring_matrix
    {{
        //A   _   C   G  U  AAN  ACN  AGN  ATN  UN  CNN  GNN  TNN  UNN  N
        {{ 96, 0, -316, 0, -192, 0, 0, 0, 0, 0, 0, -396, 0, 0, 0 }}, // A_ind = 0, C_ind = 2, G_ind = 4, T_ind = 12
        {{0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{-316, 0,  100, 0, -352, 0, 0, 0, 0, 0, 0, -295, 0, 0, 0 }},
        {{0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{-192, 0, -352, 0, 98, 0, 0, 0, 0, 0, 0, -329, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{-396, 0, -295, 0, -329, 0, 0, 0, 0, 0, 0, 100, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }},
        {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }}
    }};
    
    nucleotide_scoring_scheme<int> my_scoring_scheme{scoring_matrix};

    // Configure the output:
    auto output_config =    seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_begin_position{} |
                            seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_alignment{};

    // Apply the local alignment configuration.
    auto min_cfg =          seqan3::align_cfg::method_local{} |
                            seqan3::align_cfg::scoring_scheme{my_scoring_scheme} | output_config |
                            seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-80}, seqan3::align_cfg::extension_score{-40}}|
                            seqan3::align_cfg::parallel{32}; // Set gap open and extension scores
                            

    // Score
    int max_score = std::numeric_limits<int>::min();
    size_t max_scoring_ind = 0; 


    try
    {
        for (size_t i = 0; i < sequences.size(); ++i)
        {
        auto& [barcode_name, seq_pair] = sequences[i];
        auto& [seq_pair1, seq_pair2] = seq_pair;

            // Assuming that seq_pair1 and seq_pair2 are vectors, 
            // use them directly without '.second'
            for (auto res : seqan3::align_pairwise(std::tie(seq_pair1, seq_pair2), min_cfg))
            {
                // seqan3::debug_stream << "Score: " << res.score() << '\n'
                //                     << "Alignment: \n" << res.alignment() << '\n'
                //                     << "Second Sequence: \n" << seq_pair2 << '\n';


                if (res.score() > max_score)
                {
                    max_score = res.score();
                    max_scoring_ind = i;
                }
            }
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << "An error occurred: " << e.what() << std::endl;
    }

    return {max_score, max_scoring_ind};
}