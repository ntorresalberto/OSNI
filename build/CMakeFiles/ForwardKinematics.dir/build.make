# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/andrea/Desktop/code/OSNI

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/andrea/Desktop/code/OSNI/build

# Include any dependencies generated for this target.
include CMakeFiles/ForwardKinematics.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/ForwardKinematics.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/ForwardKinematics.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ForwardKinematics.dir/flags.make

CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.o: CMakeFiles/ForwardKinematics.dir/flags.make
CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.o: ../src/quaternion_integrator.cpp
CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.o: CMakeFiles/ForwardKinematics.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/andrea/Desktop/code/OSNI/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.o -MF CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.o.d -o CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.o -c /home/andrea/Desktop/code/OSNI/src/quaternion_integrator.cpp

CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/andrea/Desktop/code/OSNI/src/quaternion_integrator.cpp > CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.i

CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/andrea/Desktop/code/OSNI/src/quaternion_integrator.cpp -o CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.s

# Object files for target ForwardKinematics
ForwardKinematics_OBJECTS = \
"CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.o"

# External object files for target ForwardKinematics
ForwardKinematics_EXTERNAL_OBJECTS =

libForwardKinematics.a: CMakeFiles/ForwardKinematics.dir/src/quaternion_integrator.cpp.o
libForwardKinematics.a: CMakeFiles/ForwardKinematics.dir/build.make
libForwardKinematics.a: CMakeFiles/ForwardKinematics.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/andrea/Desktop/code/OSNI/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libForwardKinematics.a"
	$(CMAKE_COMMAND) -P CMakeFiles/ForwardKinematics.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ForwardKinematics.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ForwardKinematics.dir/build: libForwardKinematics.a
.PHONY : CMakeFiles/ForwardKinematics.dir/build

CMakeFiles/ForwardKinematics.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ForwardKinematics.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ForwardKinematics.dir/clean

CMakeFiles/ForwardKinematics.dir/depend:
	cd /home/andrea/Desktop/code/OSNI/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/andrea/Desktop/code/OSNI /home/andrea/Desktop/code/OSNI /home/andrea/Desktop/code/OSNI/build /home/andrea/Desktop/code/OSNI/build /home/andrea/Desktop/code/OSNI/build/CMakeFiles/ForwardKinematics.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ForwardKinematics.dir/depend
