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
include CMakeFiles/svd.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/svd.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/svd.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/svd.dir/flags.make

CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.o: CMakeFiles/svd.dir/flags.make
CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.o: /home/quaret/Documents/C++/SVD-project/src/GolubKahanSVD.cpp
CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.o: CMakeFiles/svd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/quaret/Documents/C++/SVD-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.o -MF CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.o.d -o CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.o -c /home/quaret/Documents/C++/SVD-project/src/GolubKahanSVD.cpp

CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quaret/Documents/C++/SVD-project/src/GolubKahanSVD.cpp > CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.i

CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quaret/Documents/C++/SVD-project/src/GolubKahanSVD.cpp -o CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.s

CMakeFiles/svd.dir/src/iterative_refinement.cpp.o: CMakeFiles/svd.dir/flags.make
CMakeFiles/svd.dir/src/iterative_refinement.cpp.o: /home/quaret/Documents/C++/SVD-project/src/iterative_refinement.cpp
CMakeFiles/svd.dir/src/iterative_refinement.cpp.o: CMakeFiles/svd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/quaret/Documents/C++/SVD-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/svd.dir/src/iterative_refinement.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/svd.dir/src/iterative_refinement.cpp.o -MF CMakeFiles/svd.dir/src/iterative_refinement.cpp.o.d -o CMakeFiles/svd.dir/src/iterative_refinement.cpp.o -c /home/quaret/Documents/C++/SVD-project/src/iterative_refinement.cpp

CMakeFiles/svd.dir/src/iterative_refinement.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/svd.dir/src/iterative_refinement.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quaret/Documents/C++/SVD-project/src/iterative_refinement.cpp > CMakeFiles/svd.dir/src/iterative_refinement.cpp.i

CMakeFiles/svd.dir/src/iterative_refinement.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/svd.dir/src/iterative_refinement.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quaret/Documents/C++/SVD-project/src/iterative_refinement.cpp -o CMakeFiles/svd.dir/src/iterative_refinement.cpp.s

CMakeFiles/svd.dir/src/jacobi.cpp.o: CMakeFiles/svd.dir/flags.make
CMakeFiles/svd.dir/src/jacobi.cpp.o: /home/quaret/Documents/C++/SVD-project/src/jacobi.cpp
CMakeFiles/svd.dir/src/jacobi.cpp.o: CMakeFiles/svd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/quaret/Documents/C++/SVD-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/svd.dir/src/jacobi.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/svd.dir/src/jacobi.cpp.o -MF CMakeFiles/svd.dir/src/jacobi.cpp.o.d -o CMakeFiles/svd.dir/src/jacobi.cpp.o -c /home/quaret/Documents/C++/SVD-project/src/jacobi.cpp

CMakeFiles/svd.dir/src/jacobi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/svd.dir/src/jacobi.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quaret/Documents/C++/SVD-project/src/jacobi.cpp > CMakeFiles/svd.dir/src/jacobi.cpp.i

CMakeFiles/svd.dir/src/jacobi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/svd.dir/src/jacobi.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quaret/Documents/C++/SVD-project/src/jacobi.cpp -o CMakeFiles/svd.dir/src/jacobi.cpp.s

CMakeFiles/svd.dir/src/generate_svd.cpp.o: CMakeFiles/svd.dir/flags.make
CMakeFiles/svd.dir/src/generate_svd.cpp.o: /home/quaret/Documents/C++/SVD-project/src/generate_svd.cpp
CMakeFiles/svd.dir/src/generate_svd.cpp.o: CMakeFiles/svd.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/quaret/Documents/C++/SVD-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/svd.dir/src/generate_svd.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/svd.dir/src/generate_svd.cpp.o -MF CMakeFiles/svd.dir/src/generate_svd.cpp.o.d -o CMakeFiles/svd.dir/src/generate_svd.cpp.o -c /home/quaret/Documents/C++/SVD-project/src/generate_svd.cpp

CMakeFiles/svd.dir/src/generate_svd.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/svd.dir/src/generate_svd.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quaret/Documents/C++/SVD-project/src/generate_svd.cpp > CMakeFiles/svd.dir/src/generate_svd.cpp.i

CMakeFiles/svd.dir/src/generate_svd.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/svd.dir/src/generate_svd.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quaret/Documents/C++/SVD-project/src/generate_svd.cpp -o CMakeFiles/svd.dir/src/generate_svd.cpp.s

# Object files for target svd
svd_OBJECTS = \
"CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.o" \
"CMakeFiles/svd.dir/src/iterative_refinement.cpp.o" \
"CMakeFiles/svd.dir/src/jacobi.cpp.o" \
"CMakeFiles/svd.dir/src/generate_svd.cpp.o"

# External object files for target svd
svd_EXTERNAL_OBJECTS =

libsvd.a: CMakeFiles/svd.dir/src/GolubKahanSVD.cpp.o
libsvd.a: CMakeFiles/svd.dir/src/iterative_refinement.cpp.o
libsvd.a: CMakeFiles/svd.dir/src/jacobi.cpp.o
libsvd.a: CMakeFiles/svd.dir/src/generate_svd.cpp.o
libsvd.a: CMakeFiles/svd.dir/build.make
libsvd.a: CMakeFiles/svd.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/quaret/Documents/C++/SVD-project/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libsvd.a"
	$(CMAKE_COMMAND) -P CMakeFiles/svd.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/svd.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/svd.dir/build: libsvd.a
.PHONY : CMakeFiles/svd.dir/build

CMakeFiles/svd.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/svd.dir/cmake_clean.cmake
.PHONY : CMakeFiles/svd.dir/clean

CMakeFiles/svd.dir/depend:
	cd /home/quaret/Documents/C++/SVD-project/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/quaret/Documents/C++/SVD-project /home/quaret/Documents/C++/SVD-project /home/quaret/Documents/C++/SVD-project/build /home/quaret/Documents/C++/SVD-project/build /home/quaret/Documents/C++/SVD-project/build/CMakeFiles/svd.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/svd.dir/depend

