# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/cmake-3.16.7/bin/cmake

# The command to remove a file.
RM = /opt/cmake-3.16.7/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/peter/Desktop/UNIPD/MuLifetime/MC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/peter/Desktop/UNIPD/MuLifetime/MC

# Include any dependencies generated for this target.
include CMakeFiles/main.exe.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main.exe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.exe.dir/flags.make

CMakeFiles/main.exe.dir/montecarlosim.C.o: CMakeFiles/main.exe.dir/flags.make
CMakeFiles/main.exe.dir/montecarlosim.C.o: montecarlosim.C
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/peter/Desktop/UNIPD/MuLifetime/MC/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.exe.dir/montecarlosim.C.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/main.exe.dir/montecarlosim.C.o -c /home/peter/Desktop/UNIPD/MuLifetime/MC/montecarlosim.C

CMakeFiles/main.exe.dir/montecarlosim.C.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.exe.dir/montecarlosim.C.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/peter/Desktop/UNIPD/MuLifetime/MC/montecarlosim.C > CMakeFiles/main.exe.dir/montecarlosim.C.i

CMakeFiles/main.exe.dir/montecarlosim.C.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.exe.dir/montecarlosim.C.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/peter/Desktop/UNIPD/MuLifetime/MC/montecarlosim.C -o CMakeFiles/main.exe.dir/montecarlosim.C.s

# Object files for target main.exe
main_exe_OBJECTS = \
"CMakeFiles/main.exe.dir/montecarlosim.C.o"

# External object files for target main.exe
main_exe_EXTERNAL_OBJECTS =

main.exe: CMakeFiles/main.exe.dir/montecarlosim.C.o
main.exe: CMakeFiles/main.exe.dir/build.make
main.exe: libgenerator.so
main.exe: /home/peter/CERN_ROOT/root/lib/libGpad.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libGraf.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libHist.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libMatrix.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libTree.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libMathCore.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libImt.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libMultiProc.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libNet.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libRIO.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libThread.so.6.24.08
main.exe: /home/peter/CERN_ROOT/root/lib/libCore.so.6.24.08
main.exe: CMakeFiles/main.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/peter/Desktop/UNIPD/MuLifetime/MC/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable main.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.exe.dir/build: main.exe

.PHONY : CMakeFiles/main.exe.dir/build

CMakeFiles/main.exe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.exe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.exe.dir/clean

CMakeFiles/main.exe.dir/depend:
	cd /home/peter/Desktop/UNIPD/MuLifetime/MC && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/peter/Desktop/UNIPD/MuLifetime/MC /home/peter/Desktop/UNIPD/MuLifetime/MC /home/peter/Desktop/UNIPD/MuLifetime/MC /home/peter/Desktop/UNIPD/MuLifetime/MC /home/peter/Desktop/UNIPD/MuLifetime/MC/CMakeFiles/main.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.exe.dir/depend

