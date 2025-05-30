cmake_minimum_required (VERSION 3.14)
 
project (my_app LANGUAGES CXX VERSION 1.0.0)
 
set (seqan3_git_tag "3.3.0") # adapt as needed, e.g. "3.2.0" or "master"

# Place binaries in a bin/ folder inside the build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -flto")
message (STATUS "Fetching SeqAn3 ${seqan3_git_tag}:")
 
include (FetchContent)
FetchContent_Declare (
    seqan3_fetch_content
    GIT_REPOSITORY "https://github.com/seqan/seqan3.git"
    GIT_TAG "${seqan3_git_tag}"
)
 
# Download and make SeqAn3 available.
FetchContent_MakeAvailable (seqan3_fetch_content)

list(APPEND CMAKE_PREFIX_PATH "/bin")

# add_executable(seq_io seq_io.cpp)
# target_link_libraries(seq_io PUBLIC seqan3::seqan3)

# add_executable(SW_align seq_io.cpp io.cpp alignments.cpp)
# target_link_libraries(SW_align PUBLIC seqan3::seqan3)

# add_executable(test_main io.cpp alignments.cpp test_main.cpp)
# target_link_libraries(test_main PUBLIC seqan3::seqan3)

add_executable(demultiplex io.cpp alignments.cpp demultiplex.cpp)
target_link_libraries(demultiplex PUBLIC seqan3::seqan3)



