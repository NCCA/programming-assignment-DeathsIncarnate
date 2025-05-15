Jacques Moss s5616052 CFGAA Assigment
------------------------------------------------------------------------------------------------------------------------------------------------------
## Report for my Assignement:


------------------------------------------------------------------------------------------------------------------------------------------------------
## Initial Idea:
I wanted to create an SPH fluid simulation that applies the physics off of water and how it works in real life and allows a particle simulation to replicate and mimic those movements and algorithms.

------------------------------------------------------------------------------------------------------------------------------------------------------
##  Video:


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
- 

------------------------------------------------------------------------------------------------------------------------------------------------------
## Structure of Directory:


------------------------------------------------------------------------------------------------------------------------------------------------------
## Annotated Bibliography and Research:

------------------------------------------------------------------------------------------------------------------------------------------------------



