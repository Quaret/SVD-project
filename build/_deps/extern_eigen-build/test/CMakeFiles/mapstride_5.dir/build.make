# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/quaret/Documents/C++/SVD-project

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/quaret/Documents/C++/SVD-project/build

# Include any dependencies generated for this target.
include _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/compiler_depend.make

# Include the progress variables for this target.
include _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/progress.make

# Include the compile flags for this target's objects.
include _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/flags.make

_deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/mapstride.cpp.o: _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/flags.make
_deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/mapstride.cpp.o: _deps/extern_eigen-src/test/mapstride.cpp
_deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/mapstride.cpp.o: _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/quaret/Documents/C++/SVD-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/mapstride.cpp.o"
	cd /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/mapstride.cpp.o -MF CMakeFiles/mapstride_5.dir/mapstride.cpp.o.d -o CMakeFiles/mapstride_5.dir/mapstride.cpp.o -c /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-src/test/mapstride.cpp

_deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/mapstride.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/mapstride_5.dir/mapstride.cpp.i"
	cd /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-src/test/mapstride.cpp > CMakeFiles/mapstride_5.dir/mapstride.cpp.i

_deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/mapstride.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/mapstride_5.dir/mapstride.cpp.s"
	cd /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-src/test/mapstride.cpp -o CMakeFiles/mapstride_5.dir/mapstride.cpp.s

# Object files for target mapstride_5
mapstride_5_OBJECTS = \
"CMakeFiles/mapstride_5.dir/mapstride.cpp.o"

# External object files for target mapstride_5
mapstride_5_EXTERNAL_OBJECTS =

_deps/extern_eigen-build/test/mapstride_5: _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/mapstride.cpp.o
_deps/extern_eigen-build/test/mapstride_5: _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/build.make
_deps/extern_eigen-build/test/mapstride_5: _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/quaret/Documents/C++/SVD-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mapstride_5"
	cd /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mapstride_5.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
_deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/build: _deps/extern_eigen-build/test/mapstride_5
.PHONY : _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/build

_deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/clean:
	cd /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/test && $(CMAKE_COMMAND) -P CMakeFiles/mapstride_5.dir/cmake_clean.cmake
.PHONY : _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/clean

_deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/depend:
	cd /home/quaret/Documents/C++/SVD-project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/quaret/Documents/C++/SVD-project /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-src/test /home/quaret/Documents/C++/SVD-project/build /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/test /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : _deps/extern_eigen-build/test/CMakeFiles/mapstride_5.dir/depend

