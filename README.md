# N-Body-2D-Simulator
A simple N-Body simulator made by @Dahzito (me) during my free time.
Uses Newtonian Mechanics, and for now it runs on the Euler's Method to compute the velocity and position of every body that is added to the program.

To add a new body to the program, add the following line of code, and change whatever you want:

bodies.append(Body(5.68*10**20, [0,0], [0,2380], [3.0e12,0], "Comet"))
                   #Mass, acceleration, velocity, coordinates, name;

This needs to be added after the piece of code that contains: class Body ; and the array bodies = []

Runs on Python, and I do not plan to make it run in any other coding language.
However, I plan to make it run on other numerical methods in the close future.
