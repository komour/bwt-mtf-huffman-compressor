# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = "/Users/asap/Library/Application Support/JetBrains/Toolbox/apps/CLion/ch-0/213.5744.254/CLion.app/Contents/bin/cmake/mac/bin/cmake"

# The command to remove a file.
RM = "/Users/asap/Library/Application Support/JetBrains/Toolbox/apps/CLion/ch-0/213.5744.254/CLion.app/Contents/bin/cmake/mac/bin/cmake" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/asap/CLionProjects/BWT-MTF-Huffman

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/asap/CLionProjects/BWT-MTF-Huffman/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/BWT_MTF_Huffman.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/BWT_MTF_Huffman.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/BWT_MTF_Huffman.dir/flags.make

CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.o: CMakeFiles/BWT_MTF_Huffman.dir/flags.make
CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/asap/CLionProjects/BWT-MTF-Huffman/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.o -c /Users/asap/CLionProjects/BWT-MTF-Huffman/main.cpp

CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/asap/CLionProjects/BWT-MTF-Huffman/main.cpp > CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.i

CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/asap/CLionProjects/BWT-MTF-Huffman/main.cpp -o CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.s

# Object files for target BWT_MTF_Huffman
BWT_MTF_Huffman_OBJECTS = \
"CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.o"

# External object files for target BWT_MTF_Huffman
BWT_MTF_Huffman_EXTERNAL_OBJECTS =

BWT_MTF_Huffman: CMakeFiles/BWT_MTF_Huffman.dir/main.cpp.o
BWT_MTF_Huffman: CMakeFiles/BWT_MTF_Huffman.dir/build.make
BWT_MTF_Huffman: CMakeFiles/BWT_MTF_Huffman.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/asap/CLionProjects/BWT-MTF-Huffman/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable BWT_MTF_Huffman"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BWT_MTF_Huffman.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/BWT_MTF_Huffman.dir/build: BWT_MTF_Huffman
.PHONY : CMakeFiles/BWT_MTF_Huffman.dir/build

CMakeFiles/BWT_MTF_Huffman.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/BWT_MTF_Huffman.dir/cmake_clean.cmake
.PHONY : CMakeFiles/BWT_MTF_Huffman.dir/clean

CMakeFiles/BWT_MTF_Huffman.dir/depend:
	cd /Users/asap/CLionProjects/BWT-MTF-Huffman/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/asap/CLionProjects/BWT-MTF-Huffman /Users/asap/CLionProjects/BWT-MTF-Huffman /Users/asap/CLionProjects/BWT-MTF-Huffman/cmake-build-debug /Users/asap/CLionProjects/BWT-MTF-Huffman/cmake-build-debug /Users/asap/CLionProjects/BWT-MTF-Huffman/cmake-build-debug/CMakeFiles/BWT_MTF_Huffman.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/BWT_MTF_Huffman.dir/depend

