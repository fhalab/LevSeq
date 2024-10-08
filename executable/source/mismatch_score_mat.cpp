#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/convert.hpp>
using namespace seqan3;

void perform_alignment(std::string const & seq1_str, std::string const & seq2_str)
{
    auto seq1 = seq1_str | std::views::transform([](char c) { return seqan3::dna15{}.assign_char(c); });
    auto seq2 = seq2_str | std::views::transform([](char c) { return seqan3::dna15{}.assign_char(c); });
    
    // Define a custom 15x15 scoring matrix for dna15.
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
    auto output_config = seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_begin_position{}
                       | seqan3::align_cfg::output_end_position{} | seqan3::align_cfg::output_alignment{};

    // Apply the local alignment configuration.
    auto min_cfg = align_cfg::method_local{} |
                   align_cfg::scoring_scheme{my_scoring_scheme} | output_config;

    // Compute the pairwise alignment and output the score and the alignment.
    for (auto res : seqan3::align_pairwise(std::tie(seq1, seq2), min_cfg))
     {
        seqan3::debug_stream << "Score: " << res.score() << '\n';
        seqan3::debug_stream << "Begin: (" << res.sequence1_begin_position() << "," << res.sequence2_begin_position()
                             << ")\n";
        seqan3::debug_stream << "End: (" << res.sequence1_end_position() << "," << res.sequence2_end_position()
                             << ")\n";
        seqan3::debug_stream << "Alignment: \n" << res.alignment() << '\n';
    }
}


int main();