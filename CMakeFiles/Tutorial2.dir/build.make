# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/geant4/Geant4/Scint_SensitiveDetector

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/geant4/Geant4/Scint_SensitiveDetector

# Include any dependencies generated for this target.
include CMakeFiles/Tutorial2.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Tutorial2.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Tutorial2.dir/flags.make

CMakeFiles/Tutorial2.dir/Tutorial2.cc.o: CMakeFiles/Tutorial2.dir/flags.make
CMakeFiles/Tutorial2.dir/Tutorial2.cc.o: Tutorial2.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/geant4/Geant4/Scint_SensitiveDetector/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/Tutorial2.dir/Tutorial2.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Tutorial2.dir/Tutorial2.cc.o -c /home/geant4/Geant4/Scint_SensitiveDetector/Tutorial2.cc

CMakeFiles/Tutorial2.dir/Tutorial2.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tutorial2.dir/Tutorial2.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/geant4/Geant4/Scint_SensitiveDetector/Tutorial2.cc > CMakeFiles/Tutorial2.dir/Tutorial2.cc.i

CMakeFiles/Tutorial2.dir/Tutorial2.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tutorial2.dir/Tutorial2.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/geant4/Geant4/Scint_SensitiveDetector/Tutorial2.cc -o CMakeFiles/Tutorial2.dir/Tutorial2.cc.s

CMakeFiles/Tutorial2.dir/Tutorial2.cc.o.requires:
.PHONY : CMakeFiles/Tutorial2.dir/Tutorial2.cc.o.requires

CMakeFiles/Tutorial2.dir/Tutorial2.cc.o.provides: CMakeFiles/Tutorial2.dir/Tutorial2.cc.o.requires
	$(MAKE) -f CMakeFiles/Tutorial2.dir/build.make CMakeFiles/Tutorial2.dir/Tutorial2.cc.o.provides.build
.PHONY : CMakeFiles/Tutorial2.dir/Tutorial2.cc.o.provides

CMakeFiles/Tutorial2.dir/Tutorial2.cc.o.provides.build: CMakeFiles/Tutorial2.dir/Tutorial2.cc.o

CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o: CMakeFiles/Tutorial2.dir/flags.make
CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o: src/SensitiveDetector.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/geant4/Geant4/Scint_SensitiveDetector/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o -c /home/geant4/Geant4/Scint_SensitiveDetector/src/SensitiveDetector.cc

CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/geant4/Geant4/Scint_SensitiveDetector/src/SensitiveDetector.cc > CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.i

CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/geant4/Geant4/Scint_SensitiveDetector/src/SensitiveDetector.cc -o CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.s

CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o.requires:
.PHONY : CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o.requires

CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o.provides: CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o.requires
	$(MAKE) -f CMakeFiles/Tutorial2.dir/build.make CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o.provides.build
.PHONY : CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o.provides

CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o.provides.build: CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o

CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o: CMakeFiles/Tutorial2.dir/flags.make
CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o: src/PhysicsList.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/geant4/Geant4/Scint_SensitiveDetector/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o -c /home/geant4/Geant4/Scint_SensitiveDetector/src/PhysicsList.cc

CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/geant4/Geant4/Scint_SensitiveDetector/src/PhysicsList.cc > CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.i

CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/geant4/Geant4/Scint_SensitiveDetector/src/PhysicsList.cc -o CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.s

CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o.requires:
.PHONY : CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o.requires

CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o.provides: CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o.requires
	$(MAKE) -f CMakeFiles/Tutorial2.dir/build.make CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o.provides.build
.PHONY : CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o.provides

CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o.provides.build: CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o

CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o: CMakeFiles/Tutorial2.dir/flags.make
CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o: src/ActionInitialization.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/geant4/Geant4/Scint_SensitiveDetector/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o -c /home/geant4/Geant4/Scint_SensitiveDetector/src/ActionInitialization.cc

CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/geant4/Geant4/Scint_SensitiveDetector/src/ActionInitialization.cc > CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.i

CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/geant4/Geant4/Scint_SensitiveDetector/src/ActionInitialization.cc -o CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.s

CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o.requires:
.PHONY : CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o.requires

CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o.provides: CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o.requires
	$(MAKE) -f CMakeFiles/Tutorial2.dir/build.make CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o.provides.build
.PHONY : CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o.provides

CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o.provides.build: CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o

CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/Tutorial2.dir/flags.make
CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o: src/PrimaryGeneratorAction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/geant4/Geant4/Scint_SensitiveDetector/CMakeFiles $(CMAKE_PROGRESS_5)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o -c /home/geant4/Geant4/Scint_SensitiveDetector/src/PrimaryGeneratorAction.cc

CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/geant4/Geant4/Scint_SensitiveDetector/src/PrimaryGeneratorAction.cc > CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.i

CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/geant4/Geant4/Scint_SensitiveDetector/src/PrimaryGeneratorAction.cc -o CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.s

CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o.requires:
.PHONY : CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o.requires

CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o.provides: CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o.requires
	$(MAKE) -f CMakeFiles/Tutorial2.dir/build.make CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o.provides.build
.PHONY : CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o.provides

CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o.provides.build: CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o

CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o: CMakeFiles/Tutorial2.dir/flags.make
CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o: src/DetectorConstruction.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/geant4/Geant4/Scint_SensitiveDetector/CMakeFiles $(CMAKE_PROGRESS_6)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o -c /home/geant4/Geant4/Scint_SensitiveDetector/src/DetectorConstruction.cc

CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/geant4/Geant4/Scint_SensitiveDetector/src/DetectorConstruction.cc > CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.i

CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/geant4/Geant4/Scint_SensitiveDetector/src/DetectorConstruction.cc -o CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.s

CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o.requires:
.PHONY : CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o.requires

CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o.provides: CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o.requires
	$(MAKE) -f CMakeFiles/Tutorial2.dir/build.make CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o.provides.build
.PHONY : CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o.provides

CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o.provides.build: CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o

# Object files for target Tutorial2
Tutorial2_OBJECTS = \
"CMakeFiles/Tutorial2.dir/Tutorial2.cc.o" \
"CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o" \
"CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o" \
"CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o" \
"CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o" \
"CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o"

# External object files for target Tutorial2
Tutorial2_EXTERNAL_OBJECTS =

Tutorial2: CMakeFiles/Tutorial2.dir/Tutorial2.cc.o
Tutorial2: CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o
Tutorial2: CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o
Tutorial2: CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o
Tutorial2: CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o
Tutorial2: CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o
Tutorial2: CMakeFiles/Tutorial2.dir/build.make
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4Tree.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4FR.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4GMocren.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4visHepRep.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4RayTracer.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4VRML.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4OpenGL.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4gl2ps.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4OpenInventor.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4vis_management.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4modeling.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4interfaces.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4persistency.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4analysis.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4error_propagation.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4readout.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4physicslists.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4run.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4event.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4tracking.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4parmodels.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4processes.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4digits_hits.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4track.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4particles.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4geometry.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4materials.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4graphics_reps.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4intercoms.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4global.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4FR.so
Tutorial2: /usr/lib/i386-linux-gnu/libXmu.so
Tutorial2: /usr/lib/i386-linux-gnu/libQtOpenGL.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4gl2ps.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4vis_management.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4modeling.so
Tutorial2: /usr/lib/i386-linux-gnu/libQtGui.so
Tutorial2: /usr/lib/i386-linux-gnu/libQtCore.so
Tutorial2: /usr/lib/i386-linux-gnu/libCoin.so
Tutorial2: /usr/lib/i386-linux-gnu/libSM.so
Tutorial2: /usr/lib/i386-linux-gnu/libICE.so
Tutorial2: /usr/lib/i386-linux-gnu/libX11.so
Tutorial2: /usr/lib/i386-linux-gnu/libXext.so
Tutorial2: /usr/lib/i386-linux-gnu/libGLU.so
Tutorial2: /usr/lib/i386-linux-gnu/libGL.so
Tutorial2: /usr/local/lib/libSoXt.so
Tutorial2: /usr/lib/i386-linux-gnu/libXm.so
Tutorial2: /usr/lib/i386-linux-gnu/libXpm.so
Tutorial2: /usr/lib/i386-linux-gnu/libxerces-c.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4run.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4event.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4tracking.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4processes.so
Tutorial2: /usr/lib/i386-linux-gnu/libexpat.so
Tutorial2: /usr/lib/i386-linux-gnu/libz.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4digits_hits.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4track.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4particles.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4geometry.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4materials.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4graphics_reps.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4intercoms.so
Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4global.so
Tutorial2: /usr/local/lib/libCLHEP.so
Tutorial2: CMakeFiles/Tutorial2.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable Tutorial2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Tutorial2.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Tutorial2.dir/build: Tutorial2
.PHONY : CMakeFiles/Tutorial2.dir/build

# Object files for target Tutorial2
Tutorial2_OBJECTS = \
"CMakeFiles/Tutorial2.dir/Tutorial2.cc.o" \
"CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o" \
"CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o" \
"CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o" \
"CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o" \
"CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o"

# External object files for target Tutorial2
Tutorial2_EXTERNAL_OBJECTS =

CMakeFiles/CMakeRelink.dir/Tutorial2: CMakeFiles/Tutorial2.dir/Tutorial2.cc.o
CMakeFiles/CMakeRelink.dir/Tutorial2: CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o
CMakeFiles/CMakeRelink.dir/Tutorial2: CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o
CMakeFiles/CMakeRelink.dir/Tutorial2: CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o
CMakeFiles/CMakeRelink.dir/Tutorial2: CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o
CMakeFiles/CMakeRelink.dir/Tutorial2: CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o
CMakeFiles/CMakeRelink.dir/Tutorial2: CMakeFiles/Tutorial2.dir/build.make
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4Tree.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4FR.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4GMocren.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4visHepRep.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4RayTracer.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4VRML.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4OpenGL.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4gl2ps.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4OpenInventor.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4vis_management.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4modeling.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4interfaces.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4persistency.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4analysis.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4error_propagation.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4readout.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4physicslists.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4run.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4event.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4tracking.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4parmodels.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4processes.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4digits_hits.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4track.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4particles.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4geometry.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4materials.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4graphics_reps.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4intercoms.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4global.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4FR.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libXmu.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libQtOpenGL.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4gl2ps.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4vis_management.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4modeling.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libQtGui.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libQtCore.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libCoin.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libSM.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libICE.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libX11.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libXext.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libGLU.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libGL.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/lib/libSoXt.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libXm.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libXpm.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libxerces-c.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4run.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4event.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4tracking.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4processes.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libexpat.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/lib/i386-linux-gnu/libz.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4digits_hits.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4track.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4particles.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4geometry.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4materials.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4graphics_reps.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4intercoms.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/Geant4/Geant4.10.00.p02/lib/libG4global.so
CMakeFiles/CMakeRelink.dir/Tutorial2: /usr/local/lib/libCLHEP.so
CMakeFiles/CMakeRelink.dir/Tutorial2: CMakeFiles/Tutorial2.dir/relink.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable CMakeFiles/CMakeRelink.dir/Tutorial2"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Tutorial2.dir/relink.txt --verbose=$(VERBOSE)

# Rule to relink during preinstall.
CMakeFiles/Tutorial2.dir/preinstall: CMakeFiles/CMakeRelink.dir/Tutorial2
.PHONY : CMakeFiles/Tutorial2.dir/preinstall

CMakeFiles/Tutorial2.dir/requires: CMakeFiles/Tutorial2.dir/Tutorial2.cc.o.requires
CMakeFiles/Tutorial2.dir/requires: CMakeFiles/Tutorial2.dir/src/SensitiveDetector.cc.o.requires
CMakeFiles/Tutorial2.dir/requires: CMakeFiles/Tutorial2.dir/src/PhysicsList.cc.o.requires
CMakeFiles/Tutorial2.dir/requires: CMakeFiles/Tutorial2.dir/src/ActionInitialization.cc.o.requires
CMakeFiles/Tutorial2.dir/requires: CMakeFiles/Tutorial2.dir/src/PrimaryGeneratorAction.cc.o.requires
CMakeFiles/Tutorial2.dir/requires: CMakeFiles/Tutorial2.dir/src/DetectorConstruction.cc.o.requires
.PHONY : CMakeFiles/Tutorial2.dir/requires

CMakeFiles/Tutorial2.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Tutorial2.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Tutorial2.dir/clean

CMakeFiles/Tutorial2.dir/depend:
	cd /home/geant4/Geant4/Scint_SensitiveDetector && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/geant4/Geant4/Scint_SensitiveDetector /home/geant4/Geant4/Scint_SensitiveDetector /home/geant4/Geant4/Scint_SensitiveDetector /home/geant4/Geant4/Scint_SensitiveDetector /home/geant4/Geant4/Scint_SensitiveDetector/CMakeFiles/Tutorial2.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Tutorial2.dir/depend
