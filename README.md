# Jacques Moss s5616052 CFGAA Assigment
------------------------------------------------------------------------------------------------------------------------------------------------------
## Report for my Assignment:

### Please use FluidSimQT as that is the final save with added QT. However FluidSim is quicker efficiency wise.
------------------------------------------------------------------------------------------------------------------------------------------------------
## Initial Idea:

I wanted to create an SPH fluid simulation that applies the physics off of water and how it works in real life and allows a particle simulation to replicate and mimic those movements and algorithms.

------------------------------------------------------------------------------------------------------------------------------------------------------
##  Video:

The first video shows how it simulates a fluid simulation quite accurately at the right recommended intialization of floats. The second video shows how the variables can be changed such as bounding box dimensions, all physics values and even starting intialization before runing with manual changing of the variables eg. m_pressureMultiplier, m_viscosityStrength and m_particleSpacing in physics.h. Particle spacing involves the amount of particles initially distributed within the square and therefore directly affects how many particles are in the scene. Also the dimensions of the bounding box can be changed in BoundingBox.h with the measurements affecting the bounds collisions as well as there VAO.

https://youtu.be/GjH93HdGTUc

https://youtu.be/XdgDH1SLI6Q 

Messing with the working variables in QT. PressureMultiplier, ViscosityStrength and yOffset. As well as the 2 buttons.

https://youtu.be/DP7xY7tK5vg 

FluidSim (Non-QT version; with faster and more accurate fluid movement but no UI controls; all changes to variables have to be manual.)

https://youtu.be/Kcz7_Tj5SKw 


------------------------------------------------------------------------------------------------------------------------------------------------------
## How to use:

- Make sure you have all of the files in FluidSimQT and make sure they all organised as it is in the github repository.
- Run "clion ." in the terminal within FluidSimQT.
- Within Clion make sure to go to your build in the top right corner and link your copyshaders to the build using cmake target.
- Press ctrl + f5/ build and run using the buttons in the top right corner.
------------------------------------------------------------------------------------------------------------------------------------------------------
## What my code does:

- Uses ParticleSOA by Jon Macey to create a particle system.
- Calculates density, pressure and viscosity and applys to the particles based off acceleration and closeness of other particles.
- Creates and draws a boundary box that repels particles when they collide and dampens the force.
- Uses a vertex and fragment shader in order to determine and normalize the density in certain areas and colour those areas more red where density is higher.
- Created a QT GUI using BlankQTNgl and qt creqator. Applying yOffset, Pressure Multiplier and Viscosity Strength based off of input values from the user.
------------------------------------------------------------------------------------------------------------------------------------------------------
## Variables and Usage:
- **pressureMultiplier** in **physics.h** effects how particles are seperated based of the amount of particles (density) within a certain radius(smoothing kernel.) Can be changed to make the particles more repulsive agaisnt eachother.
- **viscosityStrength** in **physics.h** effects the movement of the particles and how they flow in a more similiar direction when the viscosity is higher.
- **particleSpacing** effects the distant between particles when initialized and therefore the amount of particles in the intialized cube at the start. The cube will stay the same size having the same dimensions unless changed manually which it can be in the **initializeParticles** function. However the spacing between particles can be changed. The particle spacing currently directly effects the **smoothingKernel**. **Unfortunately the values from qt dont get passed through properly for this variable yet. Only can be changed manually in code**
- **m_width, m_height and m_depth** are the variables in **BoundingBox.h** that can change the physical dimensions of the m_boxVAO and how it is drawn and resolves collisions.**Unfortunately the values from qt dont get passed through properly for this variable yet. Only can be changed manually in code**
- **yOffset** found in **initializeParticles** effects how high the cube is in the y-axis when initialzed.

- There is alot of other variables but for user access these are the best variables to mess around with if you don't have any knowledge about SPH fluid simulations.

### ALL variables work properly and can be changed manually in the code before running to effect the conditions. Unfortunately not all of the qt values are able to be transferred to these variables as of yet. The variables are all the variables that can be changed to effect this simulation.


![Screenshot from 2025-05-16 10-04-41](https://github.com/user-attachments/assets/2de94072-1abc-4a23-8a89-a69aec072acd)
![Screenshot from 2025-05-16 10-05-45](https://github.com/user-attachments/assets/11c19ebb-6d49-48bb-9ee9-d5ca7eccbb3f)
![Screenshot from 2025-05-16 10-06-47](https://github.com/user-attachments/assets/37bb2fc2-045c-4a57-a989-8d4a59063a08)

-------------------------------------------------------------
## Structure of Directory:

- src: main.cpp, NGLScene.cpp, NGLSceneMouseControls.cpp, Physics.cpp, Emitter.cpp, BoundingBox.cpp, MainWindow.cpp
- include:MainWindow.h, NGLScene.h, NGLSceneMouseControls.h, Physics.h, Emitter.h, BoundingBox.h
- CMakeList.txt
- ui: MainWindow.ui
- shaders: ParticleFragment.glsl, ParticleVertex.glsl
------------------------------------------------------------------------------------------------------------------------------------------------------
## Potential improvements:

- I would've liked to include a spatial hashing system in order to speed up calculations as less not all densities would be calculated and only the densities within the smoothing kernel.
- I would've liked if my values from QT were properly passed and assessable within my functions to allow for better user interface for changing ther variables.
- I wanted to allow the user to effect the fluid with repulsion and attraction mechanics around the mouse cursor.

------------------------------------------------------------------------------------------------------------------------------------------------------
## Annotated Bibliography and Research:

Reference list
- Lague, S., 2023. Coding Adventure: Simulating Fluids [online]. www.youtube.com. Available from: https://www.youtube.com/watch?v=rSKMYc1CQHE [Accessed 15 May 2025].
- Macey, J., 2025. Jon Macey’s WebPages [online]. Jon Macey’s WebPages. Available from: https://nccastaff.bournemouth.ac.uk/jmacey/ [Accessed 15 May 2025].
- NCCA, 2024. GitHub - NCCA/BlankQtNGL [online]. GitHub. Available from: https://github.com/NCCA/BlankQtNGL [Accessed 15 May 2025].
- NCCA, 2025. GitHub - NCCA/labcode-jmacey-2: labcode-jmacey-2 created by GitHub Classroom [online]. GitHub. Available from: https://github.com/NCCA/labcode-jmacey-2 [Accessed 20 Apr 2025].
- Weaver, T. and Xiao, Z., 2016. Fluid Simulation by the Smoothed Particle Hydrodynamics Method: A Survey [online]. Available from: https://eprints.bournemouth.ac.uk/23384/1/2016%20Fluid%20simulation.pdf [Accessed 9 Apr 2025].

------------------------------------------------------------------------------------------------------------------------------------------------------



