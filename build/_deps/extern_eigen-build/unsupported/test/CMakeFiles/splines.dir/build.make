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
include _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/compiler_depend.make

# Include the progress variables for this target.
include _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/progress.make

# Include the compile flags for this target's objects.
include _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/flags.make

_deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/splines.cpp.o: _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/flags.make
_deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/splines.cpp.o: _deps/extern_eigen-src/unsupported/test/splines.cpp
_deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/splines.cpp.o: _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/quaret/Documents/C++/SVD-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/splines.cpp.o"
	cd /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/unsupported/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/splines.cpp.o -MF CMakeFiles/splines.dir/splines.cpp.o.d -o CMakeFiles/splines.dir/splines.cpp.o -c /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-src/unsupported/test/splines.cpp

_deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/splines.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/splines.dir/splines.cpp.i"
	cd /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/unsupported/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-src/unsupported/test/splines.cpp > CMakeFiles/splines.dir/splines.cpp.i

_deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/splines.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/splines.dir/splines.cpp.s"
	cd /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/unsupported/test && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-src/unsupported/test/splines.cpp -o CMakeFiles/splines.dir/splines.cpp.s

# Object files for target splines
splines_OBJECTS = \
"CMakeFiles/splines.dir/splines.cpp.o"

# External object files for target splines
splines_EXTERNAL_OBJECTS =

_deps/extern_eigen-build/unsupported/test/splines: _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/splines.cpp.o
_deps/extern_eigen-build/unsupported/test/splines: _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/build.make
_deps/extern_eigen-build/unsupported/test/splines: _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/quaret/Documents/C++/SVD-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable splines"
	cd /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/unsupported/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/splines.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
_deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/build: _deps/extern_eigen-build/unsupported/test/splines
.PHONY : _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/build

_deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/clean:
	cd /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/unsupported/test && $(CMAKE_COMMAND) -P CMakeFiles/splines.dir/cmake_clean.cmake
.PHONY : _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/clean

_deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/depend:
	cd /home/quaret/Documents/C++/SVD-project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/quaret/Documents/C++/SVD-project /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-src/unsupported/test /home/quaret/Documents/C++/SVD-project/build /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/unsupported/test /home/quaret/Documents/C++/SVD-project/build/_deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : _deps/extern_eigen-build/unsupported/test/CMakeFiles/splines.dir/depend
