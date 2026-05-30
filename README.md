# N-Body-2D-Simulator
A simple N-Body simulator made by @Dahzito (me) during my free time.
Uses Newtonian Mechanics, and runs on the Verlet Integration method to compute the velocity and position of every body that is added to the program.

To this program, I've also added, automatic star's luminosity and temperature calculation using the Mass-Luminosity relation and the Stefan-Boltzmann's law, force felt by every body in the N-Body system, a simple mass loss calculation for planets that pass a certain temperature threshold, flux received by every body in the system, and then calculation of the temperature of every planet in the N-Body system using also the Stefan-Boltzmann's law and the planet's albedo, and gravity at the surface of every body to eventually see the Roche limit.

To add a new body to the program, add the following line of code, and change whatever you want:

bodies.append(Body(1.99*10**30, [0,0], [0, -25080], [-1.5e11, 0], "Star 1", 6.96e8, 5878, 0.0, 1000, 0, "Star"))

bodies.append(Body(2.54*10**24, [0,0], [0, 45080], [2.5e11, 0], "Planet 1", 6.3e6, 280, 0.15, 2500, 0, "Planet"))
#Mass, acceleration, velocity, coordinates, name, radii, temperature, albedo, temp_threshold, flux, type;

Or simply change the parameters of the already existing bodies.

This needs to be added after the piece of code that contains: class Body ; and the array bodies = [], btw

Runs on Python, and I do not plan to make it run in any other coding language.
