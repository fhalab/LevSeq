// alignment.hpp

#pragma once
#include <string>
#include <iostream>
#include <map>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>

//using sequence_pair_t = std::pair<seqan3::dna15_vector, seqan3::dna15_vector>;
using sequence_pair_t = std::pair<std::string, std::pair<seqan3::dna15_vector, seqan3::dna15_vector>>;
using namespace seqan3;

struct AlignmentResult {
    int max_score;
    size_t max_scoring_ind;
    std::string max_scoring_seq;
};

// Define a custom 15x15 scoring matrix for dna15.
constexpr nucleotide_scoring_scheme<int>::matrix_type scoring_matrix2
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


struct localAlignmentResult
{
    int score;
    int start_pos;
    int end_pos;
};

int perform_alignment(const std::string& seq1, const std::string& seq2);

localAlignmentResult perform_alignment_trim(std::string const & seq1_str, std::string const & seq2_str, const nucleotide_scoring_scheme<int>::matrix_type& scoring_matrix2);


std::vector<sequence_pair_t> build_sequence_pairs(seqan3::dna15_vector const& common_seq, std::map<std::string, std::string> const& barcodes);
AlignmentResult perform_multi_alignment(std::vector<sequence_pair_t> sequences);

