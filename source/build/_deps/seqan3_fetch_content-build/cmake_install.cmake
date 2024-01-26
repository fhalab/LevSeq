# Install script for directory: /home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/doc/seqan3" TYPE FILE FILES
    "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src/CHANGELOG.md"
    "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src/CODE_OF_CONDUCT.md"
    "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src/CONTRIBUTING.md"
    "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src/LICENSE.md"
    "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src/README.md"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/cmake/seqan3" TYPE FILE FILES
    "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src/build_system/seqan3-config.cmake"
    "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src/build_system/seqan3-config-version.cmake"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src/include/seqan3")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/seqan3/submodules/sdsl-lite" TYPE DIRECTORY FILES "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src/submodules/sdsl-lite/include")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/seqan3/submodules/cereal" TYPE DIRECTORY FILES "/home/emre/github_repo/MinION/source/build/_deps/seqan3_fetch_content-src/submodules/cereal/include")
endif()

