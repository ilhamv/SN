# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ilham/Work/MC

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ilham/Work/MC/cmake

# Include any dependencies generated for this target.
include CMakeFiles/MC.exe.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/MC.exe.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MC.exe.dir/flags.make

CMakeFiles/MC.exe.dir/Main.cpp.o: CMakeFiles/MC.exe.dir/flags.make
CMakeFiles/MC.exe.dir/Main.cpp.o: ../Main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ilham/Work/MC/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/MC.exe.dir/Main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/MC.exe.dir/Main.cpp.o -c /home/ilham/Work/MC/Main.cpp

CMakeFiles/MC.exe.dir/Main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/MC.exe.dir/Main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ilham/Work/MC/Main.cpp > CMakeFiles/MC.exe.dir/Main.cpp.i

CMakeFiles/MC.exe.dir/Main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/MC.exe.dir/Main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ilham/Work/MC/Main.cpp -o CMakeFiles/MC.exe.dir/Main.cpp.s

CMakeFiles/MC.exe.dir/Main.cpp.o.requires:

.PHONY : CMakeFiles/MC.exe.dir/Main.cpp.o.requires

CMakeFiles/MC.exe.dir/Main.cpp.o.provides: CMakeFiles/MC.exe.dir/Main.cpp.o.requires
	$(MAKE) -f CMakeFiles/MC.exe.dir/build.make CMakeFiles/MC.exe.dir/Main.cpp.o.provides.build
.PHONY : CMakeFiles/MC.exe.dir/Main.cpp.o.provides

CMakeFiles/MC.exe.dir/Main.cpp.o.provides.build: CMakeFiles/MC.exe.dir/Main.cpp.o


# Object files for target MC.exe
MC_exe_OBJECTS = \
"CMakeFiles/MC.exe.dir/Main.cpp.o"

# External object files for target MC.exe
MC_exe_EXTERNAL_OBJECTS =

MC.exe: CMakeFiles/MC.exe.dir/Main.cpp.o
MC.exe: CMakeFiles/MC.exe.dir/build.make
MC.exe: liblibMC.a
MC.exe: /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5_cpp.so
MC.exe: /usr/lib/x86_64-linux-gnu/hdf5/serial/lib/libhdf5.so
MC.exe: libpugixml.a
MC.exe: CMakeFiles/MC.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ilham/Work/MC/cmake/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable MC.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MC.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MC.exe.dir/build: MC.exe

.PHONY : CMakeFiles/MC.exe.dir/build

CMakeFiles/MC.exe.dir/requires: CMakeFiles/MC.exe.dir/Main.cpp.o.requires

.PHONY : CMakeFiles/MC.exe.dir/requires

CMakeFiles/MC.exe.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MC.exe.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MC.exe.dir/clean

CMakeFiles/MC.exe.dir/depend:
	cd /home/ilham/Work/MC/cmake && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ilham/Work/MC /home/ilham/Work/MC /home/ilham/Work/MC/cmake /home/ilham/Work/MC/cmake /home/ilham/Work/MC/cmake/CMakeFiles/MC.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MC.exe.dir/depend

