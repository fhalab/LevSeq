
if(NOT "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src/seqan3_fetch_content-populate-stamp/seqan3_fetch_content-populate-gitinfo.txt" IS_NEWER_THAN "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src/seqan3_fetch_content-populate-stamp/seqan3_fetch_content-populate-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: '/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src/seqan3_fetch_content-populate-stamp/seqan3_fetch_content-populate-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E remove_directory "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: '/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "/usr/bin/git"  clone --no-checkout "https://github.com/seqan/seqan3.git" "seqan3_fetch_content-src"
    WORKING_DIRECTORY "/home/emre/github_repo/MinION/source/build/_deps"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/seqan/seqan3.git'")
endif()

execute_process(
  COMMAND "/usr/bin/git"  checkout master --
  WORKING_DIRECTORY "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'master'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "/usr/bin/git"  submodule update --recursive --init 
    WORKING_DIRECTORY "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: '/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src/seqan3_fetch_content-populate-stamp/seqan3_fetch_content-populate-gitinfo.txt"
    "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src/seqan3_fetch_content-populate-stamp/seqan3_fetch_content-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: '/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src/seqan3_fetch_content-populate-stamp/seqan3_fetch_content-populate-gitclone-lastrun.txt'")
endif()

