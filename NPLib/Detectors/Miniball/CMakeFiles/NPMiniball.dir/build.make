# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.24.2/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.24.2/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/johannesheines/nptool/NPLib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/johannesheines/nptool/NPLib

# Include any dependencies generated for this target.
include Detectors/Miniball/CMakeFiles/NPMiniball.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include Detectors/Miniball/CMakeFiles/NPMiniball.dir/compiler_depend.make

# Include the progress variables for this target.
include Detectors/Miniball/CMakeFiles/NPMiniball.dir/progress.make

# Include the compile flags for this target's objects.
include Detectors/Miniball/CMakeFiles/NPMiniball.dir/flags.make

Detectors/Miniball/TMiniballDataDict.cxx: Detectors/Miniball/TMiniballData.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/johannesheines/nptool/NPLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating TMiniballDataDict.cxx"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && ../../scripts/build_dict.sh TMiniballData.h TMiniballDataDict.cxx TMiniballData.rootmap libNPMiniball.dylib

Detectors/Miniball/TMiniballPhysicsDict.cxx: Detectors/Miniball/TMiniballPhysics.h
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/johannesheines/nptool/NPLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating TMiniballPhysicsDict.cxx"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && ../../scripts/build_dict.sh TMiniballPhysics.h TMiniballPhysicsDict.cxx TMiniballPhysics.rootmap libNPMiniball.dylib

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.o: Detectors/Miniball/CMakeFiles/NPMiniball.dir/flags.make
Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.o: Detectors/Miniball/TMiniballSpectra.cxx
Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.o: Detectors/Miniball/CMakeFiles/NPMiniball.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/johannesheines/nptool/NPLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.o"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.o -MF CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.o.d -o CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.o -c /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballSpectra.cxx

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.i"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballSpectra.cxx > CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.i

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.s"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballSpectra.cxx -o CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.s

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballData.cxx.o: Detectors/Miniball/CMakeFiles/NPMiniball.dir/flags.make
Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballData.cxx.o: Detectors/Miniball/TMiniballData.cxx
Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballData.cxx.o: Detectors/Miniball/CMakeFiles/NPMiniball.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/johannesheines/nptool/NPLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballData.cxx.o"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballData.cxx.o -MF CMakeFiles/NPMiniball.dir/TMiniballData.cxx.o.d -o CMakeFiles/NPMiniball.dir/TMiniballData.cxx.o -c /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballData.cxx

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballData.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NPMiniball.dir/TMiniballData.cxx.i"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballData.cxx > CMakeFiles/NPMiniball.dir/TMiniballData.cxx.i

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballData.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NPMiniball.dir/TMiniballData.cxx.s"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballData.cxx -o CMakeFiles/NPMiniball.dir/TMiniballData.cxx.s

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.o: Detectors/Miniball/CMakeFiles/NPMiniball.dir/flags.make
Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.o: Detectors/Miniball/TMiniballPhysics.cxx
Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.o: Detectors/Miniball/CMakeFiles/NPMiniball.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/johannesheines/nptool/NPLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.o"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.o -MF CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.o.d -o CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.o -c /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballPhysics.cxx

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.i"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballPhysics.cxx > CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.i

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.s"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballPhysics.cxx -o CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.s

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.o: Detectors/Miniball/CMakeFiles/NPMiniball.dir/flags.make
Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.o: Detectors/Miniball/TMiniballDataDict.cxx
Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.o: Detectors/Miniball/CMakeFiles/NPMiniball.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/johannesheines/nptool/NPLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.o"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.o -MF CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.o.d -o CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.o -c /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballDataDict.cxx

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.i"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballDataDict.cxx > CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.i

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.s"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballDataDict.cxx -o CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.s

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.o: Detectors/Miniball/CMakeFiles/NPMiniball.dir/flags.make
Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.o: Detectors/Miniball/TMiniballPhysicsDict.cxx
Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.o: Detectors/Miniball/CMakeFiles/NPMiniball.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/johannesheines/nptool/NPLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.o"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.o -MF CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.o.d -o CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.o -c /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballPhysicsDict.cxx

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.i"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballPhysicsDict.cxx > CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.i

Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.s"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && /Users/johannesheines/opt/miniconda3/envs/CERN/bin/x86_64-apple-darwin13.4.0-clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/johannesheines/nptool/NPLib/Detectors/Miniball/TMiniballPhysicsDict.cxx -o CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.s

# Object files for target NPMiniball
NPMiniball_OBJECTS = \
"CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.o" \
"CMakeFiles/NPMiniball.dir/TMiniballData.cxx.o" \
"CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.o" \
"CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.o" \
"CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.o"

# External object files for target NPMiniball
NPMiniball_EXTERNAL_OBJECTS =

lib/libNPMiniball.dylib: Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballSpectra.cxx.o
lib/libNPMiniball.dylib: Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballData.cxx.o
lib/libNPMiniball.dylib: Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysics.cxx.o
lib/libNPMiniball.dylib: Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballDataDict.cxx.o
lib/libNPMiniball.dylib: Detectors/Miniball/CMakeFiles/NPMiniball.dir/TMiniballPhysicsDict.cxx.o
lib/libNPMiniball.dylib: Detectors/Miniball/CMakeFiles/NPMiniball.dir/build.make
lib/libNPMiniball.dylib: lib/libNPCore.dylib
lib/libNPMiniball.dylib: Detectors/Miniball/CMakeFiles/NPMiniball.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/johannesheines/nptool/NPLib/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX shared library ../../lib/libNPMiniball.dylib"
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NPMiniball.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Detectors/Miniball/CMakeFiles/NPMiniball.dir/build: lib/libNPMiniball.dylib
.PHONY : Detectors/Miniball/CMakeFiles/NPMiniball.dir/build

Detectors/Miniball/CMakeFiles/NPMiniball.dir/clean:
	cd /Users/johannesheines/nptool/NPLib/Detectors/Miniball && $(CMAKE_COMMAND) -P CMakeFiles/NPMiniball.dir/cmake_clean.cmake
.PHONY : Detectors/Miniball/CMakeFiles/NPMiniball.dir/clean

Detectors/Miniball/CMakeFiles/NPMiniball.dir/depend: Detectors/Miniball/TMiniballDataDict.cxx
Detectors/Miniball/CMakeFiles/NPMiniball.dir/depend: Detectors/Miniball/TMiniballPhysicsDict.cxx
	cd /Users/johannesheines/nptool/NPLib && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/johannesheines/nptool/NPLib /Users/johannesheines/nptool/NPLib/Detectors/Miniball /Users/johannesheines/nptool/NPLib /Users/johannesheines/nptool/NPLib/Detectors/Miniball /Users/johannesheines/nptool/NPLib/Detectors/Miniball/CMakeFiles/NPMiniball.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Detectors/Miniball/CMakeFiles/NPMiniball.dir/depend

