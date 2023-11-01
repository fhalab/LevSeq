cmake_minimum_required (VERSION 3.14)
 
project (my_app LANGUAGES CXX VERSION 1.0.0)
 
set (seqan3_git_tag "master") # adapt as needed, e.g. "3.2.0" or "master"
 
message (STATUS "Fetching SeqAn3 ${seqan3_git_tag}:")
 
include (FetchContent)
FetchContent_Declare (
    seqan3_fetch_content
    GIT_REPOSITORY "https://github.com/seqan/seqan3.git"
    GIT_TAG "${seqan3_git_tag}"
)
 
# Download and make SeqAn3 available.
FetchContent_MakeAvailable (seqan3_fetch_content)

list(APPEND CMAKE_PREFIX_PATH "/home/emre/bin")

# add_executable(seq_io seq_io.cpp)
# target_link_libraries(seq_io PUBLIC seqan3::seqan3)

# add_executable(SW_align seq_io.cpp io.cpp alignments.cpp)
# target_link_libraries(SW_align PUBLIC seqan3::seqan3)

# add_executable(test io.cpp alignments.cpp test.cpp)
# target_link_libraries(test PUBLIC seqan3::seqan3)

# add_executable(test_main io.cpp alignments.cpp test_main.cpp)
# target_link_libraries(test_main PUBLIC seqan3::seqan3)


add_executable(test_local_alignment io.cpp alignments.cpp test_local_alignment.cpp)
target_link_libraries(test_local_alignment PUBLIC seqan3::seqan3)