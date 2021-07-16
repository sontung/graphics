# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /home/sontung/bin/cmake

# The command to remove a file.
RM = /home/sontung/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sontung/graphics/DrawSVG

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sontung/graphics/DrawSVG/debug

# Include any dependencies generated for this target.
include src/CMakeFiles/drawsvg.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include src/CMakeFiles/drawsvg.dir/compiler_depend.make

# Include the progress variables for this target.
include src/CMakeFiles/drawsvg.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/drawsvg.dir/flags.make

src/CMakeFiles/drawsvg.dir/svg.cpp.o: src/CMakeFiles/drawsvg.dir/flags.make
src/CMakeFiles/drawsvg.dir/svg.cpp.o: ../src/svg.cpp
src/CMakeFiles/drawsvg.dir/svg.cpp.o: src/CMakeFiles/drawsvg.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sontung/graphics/DrawSVG/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/drawsvg.dir/svg.cpp.o"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/drawsvg.dir/svg.cpp.o -MF CMakeFiles/drawsvg.dir/svg.cpp.o.d -o CMakeFiles/drawsvg.dir/svg.cpp.o -c /home/sontung/graphics/DrawSVG/src/svg.cpp

src/CMakeFiles/drawsvg.dir/svg.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/drawsvg.dir/svg.cpp.i"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sontung/graphics/DrawSVG/src/svg.cpp > CMakeFiles/drawsvg.dir/svg.cpp.i

src/CMakeFiles/drawsvg.dir/svg.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/drawsvg.dir/svg.cpp.s"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sontung/graphics/DrawSVG/src/svg.cpp -o CMakeFiles/drawsvg.dir/svg.cpp.s

src/CMakeFiles/drawsvg.dir/png.cpp.o: src/CMakeFiles/drawsvg.dir/flags.make
src/CMakeFiles/drawsvg.dir/png.cpp.o: ../src/png.cpp
src/CMakeFiles/drawsvg.dir/png.cpp.o: src/CMakeFiles/drawsvg.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sontung/graphics/DrawSVG/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/CMakeFiles/drawsvg.dir/png.cpp.o"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/drawsvg.dir/png.cpp.o -MF CMakeFiles/drawsvg.dir/png.cpp.o.d -o CMakeFiles/drawsvg.dir/png.cpp.o -c /home/sontung/graphics/DrawSVG/src/png.cpp

src/CMakeFiles/drawsvg.dir/png.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/drawsvg.dir/png.cpp.i"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sontung/graphics/DrawSVG/src/png.cpp > CMakeFiles/drawsvg.dir/png.cpp.i

src/CMakeFiles/drawsvg.dir/png.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/drawsvg.dir/png.cpp.s"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sontung/graphics/DrawSVG/src/png.cpp -o CMakeFiles/drawsvg.dir/png.cpp.s

src/CMakeFiles/drawsvg.dir/texture.cpp.o: src/CMakeFiles/drawsvg.dir/flags.make
src/CMakeFiles/drawsvg.dir/texture.cpp.o: ../src/texture.cpp
src/CMakeFiles/drawsvg.dir/texture.cpp.o: src/CMakeFiles/drawsvg.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sontung/graphics/DrawSVG/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/CMakeFiles/drawsvg.dir/texture.cpp.o"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/drawsvg.dir/texture.cpp.o -MF CMakeFiles/drawsvg.dir/texture.cpp.o.d -o CMakeFiles/drawsvg.dir/texture.cpp.o -c /home/sontung/graphics/DrawSVG/src/texture.cpp

src/CMakeFiles/drawsvg.dir/texture.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/drawsvg.dir/texture.cpp.i"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sontung/graphics/DrawSVG/src/texture.cpp > CMakeFiles/drawsvg.dir/texture.cpp.i

src/CMakeFiles/drawsvg.dir/texture.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/drawsvg.dir/texture.cpp.s"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sontung/graphics/DrawSVG/src/texture.cpp -o CMakeFiles/drawsvg.dir/texture.cpp.s

src/CMakeFiles/drawsvg.dir/viewport.cpp.o: src/CMakeFiles/drawsvg.dir/flags.make
src/CMakeFiles/drawsvg.dir/viewport.cpp.o: ../src/viewport.cpp
src/CMakeFiles/drawsvg.dir/viewport.cpp.o: src/CMakeFiles/drawsvg.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sontung/graphics/DrawSVG/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/CMakeFiles/drawsvg.dir/viewport.cpp.o"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/drawsvg.dir/viewport.cpp.o -MF CMakeFiles/drawsvg.dir/viewport.cpp.o.d -o CMakeFiles/drawsvg.dir/viewport.cpp.o -c /home/sontung/graphics/DrawSVG/src/viewport.cpp

src/CMakeFiles/drawsvg.dir/viewport.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/drawsvg.dir/viewport.cpp.i"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sontung/graphics/DrawSVG/src/viewport.cpp > CMakeFiles/drawsvg.dir/viewport.cpp.i

src/CMakeFiles/drawsvg.dir/viewport.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/drawsvg.dir/viewport.cpp.s"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sontung/graphics/DrawSVG/src/viewport.cpp -o CMakeFiles/drawsvg.dir/viewport.cpp.s

src/CMakeFiles/drawsvg.dir/triangulation.cpp.o: src/CMakeFiles/drawsvg.dir/flags.make
src/CMakeFiles/drawsvg.dir/triangulation.cpp.o: ../src/triangulation.cpp
src/CMakeFiles/drawsvg.dir/triangulation.cpp.o: src/CMakeFiles/drawsvg.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sontung/graphics/DrawSVG/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object src/CMakeFiles/drawsvg.dir/triangulation.cpp.o"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/drawsvg.dir/triangulation.cpp.o -MF CMakeFiles/drawsvg.dir/triangulation.cpp.o.d -o CMakeFiles/drawsvg.dir/triangulation.cpp.o -c /home/sontung/graphics/DrawSVG/src/triangulation.cpp

src/CMakeFiles/drawsvg.dir/triangulation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/drawsvg.dir/triangulation.cpp.i"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sontung/graphics/DrawSVG/src/triangulation.cpp > CMakeFiles/drawsvg.dir/triangulation.cpp.i

src/CMakeFiles/drawsvg.dir/triangulation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/drawsvg.dir/triangulation.cpp.s"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sontung/graphics/DrawSVG/src/triangulation.cpp -o CMakeFiles/drawsvg.dir/triangulation.cpp.s

src/CMakeFiles/drawsvg.dir/software_renderer.cpp.o: src/CMakeFiles/drawsvg.dir/flags.make
src/CMakeFiles/drawsvg.dir/software_renderer.cpp.o: ../src/software_renderer.cpp
src/CMakeFiles/drawsvg.dir/software_renderer.cpp.o: src/CMakeFiles/drawsvg.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sontung/graphics/DrawSVG/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object src/CMakeFiles/drawsvg.dir/software_renderer.cpp.o"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/drawsvg.dir/software_renderer.cpp.o -MF CMakeFiles/drawsvg.dir/software_renderer.cpp.o.d -o CMakeFiles/drawsvg.dir/software_renderer.cpp.o -c /home/sontung/graphics/DrawSVG/src/software_renderer.cpp

src/CMakeFiles/drawsvg.dir/software_renderer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/drawsvg.dir/software_renderer.cpp.i"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sontung/graphics/DrawSVG/src/software_renderer.cpp > CMakeFiles/drawsvg.dir/software_renderer.cpp.i

src/CMakeFiles/drawsvg.dir/software_renderer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/drawsvg.dir/software_renderer.cpp.s"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sontung/graphics/DrawSVG/src/software_renderer.cpp -o CMakeFiles/drawsvg.dir/software_renderer.cpp.s

src/CMakeFiles/drawsvg.dir/drawsvg.cpp.o: src/CMakeFiles/drawsvg.dir/flags.make
src/CMakeFiles/drawsvg.dir/drawsvg.cpp.o: ../src/drawsvg.cpp
src/CMakeFiles/drawsvg.dir/drawsvg.cpp.o: src/CMakeFiles/drawsvg.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sontung/graphics/DrawSVG/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object src/CMakeFiles/drawsvg.dir/drawsvg.cpp.o"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/drawsvg.dir/drawsvg.cpp.o -MF CMakeFiles/drawsvg.dir/drawsvg.cpp.o.d -o CMakeFiles/drawsvg.dir/drawsvg.cpp.o -c /home/sontung/graphics/DrawSVG/src/drawsvg.cpp

src/CMakeFiles/drawsvg.dir/drawsvg.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/drawsvg.dir/drawsvg.cpp.i"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sontung/graphics/DrawSVG/src/drawsvg.cpp > CMakeFiles/drawsvg.dir/drawsvg.cpp.i

src/CMakeFiles/drawsvg.dir/drawsvg.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/drawsvg.dir/drawsvg.cpp.s"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sontung/graphics/DrawSVG/src/drawsvg.cpp -o CMakeFiles/drawsvg.dir/drawsvg.cpp.s

src/CMakeFiles/drawsvg.dir/main.cpp.o: src/CMakeFiles/drawsvg.dir/flags.make
src/CMakeFiles/drawsvg.dir/main.cpp.o: ../src/main.cpp
src/CMakeFiles/drawsvg.dir/main.cpp.o: src/CMakeFiles/drawsvg.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sontung/graphics/DrawSVG/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object src/CMakeFiles/drawsvg.dir/main.cpp.o"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT src/CMakeFiles/drawsvg.dir/main.cpp.o -MF CMakeFiles/drawsvg.dir/main.cpp.o.d -o CMakeFiles/drawsvg.dir/main.cpp.o -c /home/sontung/graphics/DrawSVG/src/main.cpp

src/CMakeFiles/drawsvg.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/drawsvg.dir/main.cpp.i"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sontung/graphics/DrawSVG/src/main.cpp > CMakeFiles/drawsvg.dir/main.cpp.i

src/CMakeFiles/drawsvg.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/drawsvg.dir/main.cpp.s"
	cd /home/sontung/graphics/DrawSVG/debug/src && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sontung/graphics/DrawSVG/src/main.cpp -o CMakeFiles/drawsvg.dir/main.cpp.s

# Object files for target drawsvg
drawsvg_OBJECTS = \
"CMakeFiles/drawsvg.dir/svg.cpp.o" \
"CMakeFiles/drawsvg.dir/png.cpp.o" \
"CMakeFiles/drawsvg.dir/texture.cpp.o" \
"CMakeFiles/drawsvg.dir/viewport.cpp.o" \
"CMakeFiles/drawsvg.dir/triangulation.cpp.o" \
"CMakeFiles/drawsvg.dir/software_renderer.cpp.o" \
"CMakeFiles/drawsvg.dir/drawsvg.cpp.o" \
"CMakeFiles/drawsvg.dir/main.cpp.o"

# External object files for target drawsvg
drawsvg_EXTERNAL_OBJECTS =

drawsvg: src/CMakeFiles/drawsvg.dir/svg.cpp.o
drawsvg: src/CMakeFiles/drawsvg.dir/png.cpp.o
drawsvg: src/CMakeFiles/drawsvg.dir/texture.cpp.o
drawsvg: src/CMakeFiles/drawsvg.dir/viewport.cpp.o
drawsvg: src/CMakeFiles/drawsvg.dir/triangulation.cpp.o
drawsvg: src/CMakeFiles/drawsvg.dir/software_renderer.cpp.o
drawsvg: src/CMakeFiles/drawsvg.dir/drawsvg.cpp.o
drawsvg: src/CMakeFiles/drawsvg.dir/main.cpp.o
drawsvg: src/CMakeFiles/drawsvg.dir/build.make
drawsvg: src/libdrawsvghdwr.a
drawsvg: ../src/reference/libdrawsvgref.a
drawsvg: /usr/lib/x86_64-linux-gnu/libfreetype.so
drawsvg: /usr/lib/x86_64-linux-gnu/libGL.so
drawsvg: /usr/lib/x86_64-linux-gnu/libGLU.so
drawsvg: CMU462/src/libCMU462.a
drawsvg: CMU462/deps/glew/libglew.a
drawsvg: CMU462/deps/glfw/src/libglfw3.a
drawsvg: src/libdrawsvghdwr.a
drawsvg: ../src/reference/libdrawsvgref.a
drawsvg: /usr/lib/x86_64-linux-gnu/libfreetype.so
drawsvg: /usr/lib/x86_64-linux-gnu/libGL.so
drawsvg: /usr/lib/x86_64-linux-gnu/libGLU.so
drawsvg: /usr/lib/x86_64-linux-gnu/libfreetype.so
drawsvg: /usr/lib/x86_64-linux-gnu/librt.so
drawsvg: /usr/lib/x86_64-linux-gnu/libm.so
drawsvg: /usr/lib/x86_64-linux-gnu/libX11.so
drawsvg: /usr/lib/x86_64-linux-gnu/libXrandr.so
drawsvg: /usr/lib/x86_64-linux-gnu/libXinerama.so
drawsvg: /usr/lib/x86_64-linux-gnu/libXxf86vm.so
drawsvg: /usr/lib/x86_64-linux-gnu/libXcursor.so
drawsvg: src/CMakeFiles/drawsvg.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sontung/graphics/DrawSVG/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Linking CXX executable ../drawsvg"
	cd /home/sontung/graphics/DrawSVG/debug/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/drawsvg.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/drawsvg.dir/build: drawsvg
.PHONY : src/CMakeFiles/drawsvg.dir/build

src/CMakeFiles/drawsvg.dir/clean:
	cd /home/sontung/graphics/DrawSVG/debug/src && $(CMAKE_COMMAND) -P CMakeFiles/drawsvg.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/drawsvg.dir/clean

src/CMakeFiles/drawsvg.dir/depend:
	cd /home/sontung/graphics/DrawSVG/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sontung/graphics/DrawSVG /home/sontung/graphics/DrawSVG/src /home/sontung/graphics/DrawSVG/debug /home/sontung/graphics/DrawSVG/debug/src /home/sontung/graphics/DrawSVG/debug/src/CMakeFiles/drawsvg.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/drawsvg.dir/depend

