# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.16

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "D:\ProgramFiles\JetBrains\CLion 2020.1\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "D:\ProgramFiles\JetBrains\CLion 2020.1\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\Tima\CompressFiles

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\Tima\CompressFiles\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/CompressFiles.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CompressFiles.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CompressFiles.dir/flags.make

CMakeFiles/CompressFiles.dir/main.cpp.obj: CMakeFiles/CompressFiles.dir/flags.make
CMakeFiles/CompressFiles.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\Tima\CompressFiles\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/CompressFiles.dir/main.cpp.obj"
	C:\PROGRA~1\mingw-w64\winlibs-x86_64-posix-seh-gcc-11.2.0-llvm-13.0.0-mingw-w64ucrt-9.0.0-r2\mingw64\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\CompressFiles.dir\main.cpp.obj -c D:\Tima\CompressFiles\main.cpp

CMakeFiles/CompressFiles.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CompressFiles.dir/main.cpp.i"
	C:\PROGRA~1\mingw-w64\winlibs-x86_64-posix-seh-gcc-11.2.0-llvm-13.0.0-mingw-w64ucrt-9.0.0-r2\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\Tima\CompressFiles\main.cpp > CMakeFiles\CompressFiles.dir\main.cpp.i

CMakeFiles/CompressFiles.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CompressFiles.dir/main.cpp.s"
	C:\PROGRA~1\mingw-w64\winlibs-x86_64-posix-seh-gcc-11.2.0-llvm-13.0.0-mingw-w64ucrt-9.0.0-r2\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\Tima\CompressFiles\main.cpp -o CMakeFiles\CompressFiles.dir\main.cpp.s

# Object files for target CompressFiles
CompressFiles_OBJECTS = \
"CMakeFiles/CompressFiles.dir/main.cpp.obj"

# External object files for target CompressFiles
CompressFiles_EXTERNAL_OBJECTS =

CompressFiles.exe: CMakeFiles/CompressFiles.dir/main.cpp.obj
CompressFiles.exe: CMakeFiles/CompressFiles.dir/build.make
CompressFiles.exe: CMakeFiles/CompressFiles.dir/linklibs.rsp
CompressFiles.exe: CMakeFiles/CompressFiles.dir/objects1.rsp
CompressFiles.exe: CMakeFiles/CompressFiles.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\Tima\CompressFiles\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable CompressFiles.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\CompressFiles.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CompressFiles.dir/build: CompressFiles.exe

.PHONY : CMakeFiles/CompressFiles.dir/build

CMakeFiles/CompressFiles.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\CompressFiles.dir\cmake_clean.cmake
.PHONY : CMakeFiles/CompressFiles.dir/clean

CMakeFiles/CompressFiles.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\Tima\CompressFiles D:\Tima\CompressFiles D:\Tima\CompressFiles\cmake-build-debug D:\Tima\CompressFiles\cmake-build-debug D:\Tima\CompressFiles\cmake-build-debug\CMakeFiles\CompressFiles.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CompressFiles.dir/depend

