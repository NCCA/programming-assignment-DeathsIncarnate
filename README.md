Jacques Moss s5616052 CFGAA Assigment
------------------------------------------------------------------------------------------------------------------------------------------------------
## Report for my Assignment:

### Please use FluidSimQT as that is the final save. However FluidSim is quicker efficiency wise.
------------------------------------------------------------------------------------------------------------------------------------------------------
## Initial Idea:

I wanted to create an SPH fluid simulation that applies the physics off of water and how it works in real life and allows a particle simulation to replicate and mimic those movements and algorithms.

------------------------------------------------------------------------------------------------------------------------------------------------------
##  Video:

https://youtu.be/GjH93HdGTUc

https://www.youtube.com/watch?v=4upU4l5SSSU

https://youtu.be/8ev3agKAYug 

FluidSim (Non-QT version; with faster and more accurate fluid movement but no UI controls; all changes to variables have to be manual.)

https://www.youtube.com/watch?v=REkTUbwJWNQ

- The first video shows how it simulates a fluid simulation quite accurately at the right recommended intialization of floats. The second video shows how the variables can be changed such as bounding box dimentions, all physics values and even starting intialization before runing with manual changing of the variables eg. m_pressureMultiplier, m_viscosityStrength and m_particleSpacing in physics.h. Particle spacing involves the amount of particles initially distributed within the square and therefore directly affects how many particles are in the scene. Also the dimensions of the bounding box can be changed in BoundingBox.h with the m,easurements affecting the bounds collisions as well as there VAO.


------------------------------------------------------------------------------------------------------------------------------------------------------
## What my code does:

- Uses ParticleSOA by Jon Macey to create a partcle system.
- Calculates density, pressure and viscosity and applys to the particles based off acceleration and closeness of other particles.
- Creates and draws a boundary box that repels particles when they collide and dampens the force.
- Uses a vertex and fragment shader in order to determine and normalize the density in certain areas and colour those areas more red where density is higher.
- Created a QT GUI that applys yOffset, Pressure Multiplier and Viscosity Strength based off of input values from the user.

------------------------------------------------------------------------------------------------------------------------------------------------------
## Potential improvements:

- I would've liked to include a spatial hashing system in order to speed up calculations as less not all densities would be calculated and only the densities within the smoothing kernel.
- I would've liked if my values from QT were properly passed and assessable within my functions to allow for better user interface for changing ther variables.

------------------------------------------------------------------------------------------------------------------------------------------------------
## Structure of Directory:

- src: main.cpp, NGLScene.cpp, NGLSceneMouseControls.cpp, Physics.cpp, Emitter.cpp, BoundingBox.cpp, MainWindow.cpp
- include:MainWindow.h, NGLScene.h, NGLSceneMouseControls.h, Physics.h, Emitter.h, BoundingBox.h
- CMakeList.txt
- ui: MainWindow.ui
- shaders: ParticleFragment.glsl, ParticleVertex.glsl

------------------------------------------------------------------------------------------------------------------------------------------------------
## Annotated Bibliography and Research:

Reference list
- Lague, S., 2023. Coding Adventure: Simulating Fluids [online]. www.youtube.com. Available from: https://www.youtube.com/watch?v=rSKMYc1CQHE [Accessed 5 Apr 2023].
- Weaver, T. and Xiao, Z., 2016. Fluid Simulation by the Smoothed Particle Hydrodynamics Method: A Survey [online]. Available from: https://eprints.bournemouth.ac.uk/23384/1/2016%20Fluid%20simulation.pdf [Accessed 9 Apr 2023].

------------------------------------------------------------------------------------------------------------------------------------------------------



