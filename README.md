README docx's in folders will contain example output images for all relevant scripts





README: Dipole Fitting

Central goal: Get coefficients (x2) to describe our cylindrical magnet as a dipole, allowing for fast computation.  

File descriptions:

Magdata.fld
Raw data output from ANSYS Maxwell FEM simulation giving a grid of magnetic field vectors over an area

Non-Normalized_Vectors+RSME_Calculation.py
Fits magnetic field vector data to the point-dipole equation
Gives a moderately accurate model for fitting a magnet to the equation of a perfect dipole
Omits a certain area of data closest to the magnetic field.  See df variable in script.
Requires: Vector data for magnetic field at various sample points
Formatting can be seen in MagData.fld
This is the automatic formatting for exporting a grid of sample points on a pane in ANSYS Maxwell
Should consider the distances in meters but that doesn’t make all that much sense.  Please verify this.

Normalized vectors.py 
Does basically the same thing as Non-Normalized_Vectors_RSME, just graphs them with normalized vectors so it’s a little easier to see field directions

There are also 2 ansys maxwell files for the FEM simulation of the cylindrical magnet.





README: Potential Energy Scripts
Descriptions and example outputs for files:

ElasticEnergyTester3.m
Finds total elastic energy in a bent OMCR
Notable parameters:
Num_segments: number of segments
EI: bending stiffness per segment
Lengths of segments
Internal angles for each segment
Assumes uniform bending and no complex internal forces
Example output



MagneticEnergyTester3.m
Gets total magnetic potential energy in the system
Starting parameters:
Dipole moment of the external magnitude
Lengths of tubes




TotalEnergyTest.m
Combines magnetic and elastic energy testing procedures



LocalMinsTester.m
Reimplements procedures from MagneticEnergyTester3 and EleasticEnergyTester3, now using them to find total energy, and the most intuitive local minimum thereof.
Starting parameters:
Number of segments (I’ve only ever tried 2 segments)
Bending stiffness per segment
Internal magnet dipole moment (was calculated from magnet data sheet)
External magnet dipole moment (calculated in dipole script, see other folder)
Initial robot stem angle
Initial positions for internal magnet angles
Lengths of segments

Uses fminunc to find the natural equilibrium from a certain start position, assuming no change in stem position or angle.

Outputs 

CallableMinsTester.m
Does the same thing as LocalMinsTester, but makes it callable from other functions
Inputs:
start_angle: Stem angle
ptart_pos Stem position
theta1_0: starting internal angle 1
theta2_0: start internal angle 2
do_plot: bool for if we want to draw the current position
Doesn’t close the drawing, meaning over multiple calls we get multiple overlaid drawings
iteration: current iteration
max_iteration: checks with current iteration to know if we should close the drawing


PathFinder3.m
Attempts to find pathing from an initial starting position, to an end starting Shape
Notable start parameters:
Start_angle: angle of stem
Start_pos: stem position
theta_current: angles for internal magnets at start
target_theta: goal end theta
This program repeatedly makes small increments (can be changed) in all directions and in both rotations, and calls callableminstester.m to check which brings us closest to our goal.  Then repeat.  If it runs into a dead end, we do a random perturbation and keep going.
Output example:


EnergyVisualizer.m
Generates a 3D Plot of the total potential energy over the angles of internal magnets
Makes a grid, repeatedly calls callableminstester.m, then plots it
Great for seeing stable configurations from a specific stem position and angle
Notable parameters:
Stem position (start_x, start_y) and angle (start_angle)
Bending stiffness per section
Length per section
Magnetic dipoles for big magnet an internal magnets
Range of angles to include 

Still in progress:
Visualization of bifurcation points.
Goal: show points between start positions that will lead to different stable configurations
Makes a grid, checks the end position for its stable configuration; checks if each position has neighbors that go to a different stable configuration.
Idea: we can use this to find minimum perturbation to permanently change the current configuration
Doesn’t really work all that well





