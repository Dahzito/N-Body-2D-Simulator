# N-Body-2D-Simulator
A simple N-Body simulator made by @Dahzito (me) during my free time.
Uses Newtonian Mechanics, and runs on the Verlet Integration method to compute the velocity and position of every body that is added to the program.

To add a new body to the program, add the following line of code, and change whatever you want:

bodies.append(Body(1.867*10**29, [0,0], [0,31320], [1.07e12, 0], "Star", 6.371e6))
#Mass, acceleration[x,y], velocity [x,y], coordinates[x,y], name, radii;

Or simply change the parameters of the already existing bodies.

This needs to be added after the piece of code that contains: class Body ; and the array bodies = [], btw

Runs on Python, and I do not plan to make it run in any other coding language.
