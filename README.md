# N-Body-2D-Simulator
A simple N-Body simulator made by @Dahzito (me) during my free time.
Uses Newtonian Mechanics, and runs on the Verlet Integration method to compute the velocity and position of every body that is added to the program.

To this program, I've also added, automatic star's luminosity and temperature calculation using the Mass-Luminosity relation and the Stefan-Boltzmann's law, force felt by every body in the N-Body system, a simple body vaporation system, flux received by every body in the n-body system, and then calculation of the temperature of every planet in it using the Stefan-Boltzmann's law, the body's albedo, surface emissivity and greenhouse factor, and gravity at the surface of every body to eventually see the Roche limit.

Stars have an automatic RGB value system calculation depending on their effective temperature, all the other bodies need a manual input for their RGBA values.

To add a new body to the program, add the following line of code, and change whatever you want:

bodies.append(Body(1.989e30, [0,0], [0, 0], [0, 0], "Sun", 6.96e8, 5878, 0.0, 1000, "Star", [1, 1, 0.7, 1], 0, 1))
#The RGB values of the Star are automatically determined depending on their surface temperature, also calculated using the Mass-Luminosity Relation, and Stefan-Boltzmann's law.

bodies.append(Body(5.97*10**24, [0,0], [0, 30290], [1.471e11, 0], "Earth", 6.371e6, 280, 0.31, 2000, "Planet", [15/255, 59/255, 141/255, 0.80], 0.46, 0.9))
bodies.append(Body(7.35*10**22, [0,0], [0, 31312], [1.471e11 + 3.633e8, 0], "Moon", 1.737e6, 280, 0.11, 2300, "Planet", [158/255, 158/255, 158/255, 0.80], 0.0, 0.95))

bodies.append(Body(1.898e27, [0,0], [0, 13070], [7.785e11, 0], "Jupiter", 6.9911e7, 120, 0.50, 2000, "Planet", [205/255, 133/255, 63/255, 0.8], 0, 0.9))
#Mass, acceleration, velocity, coordinates, name, radii, temperature, albedo, temp_threshold, type, RGB colors (0-1), greenhouse effect (0-1), ε - Surface Emissitivity (0-1);

(Or simply change the parameters of the already existing bodies)

This needs to be added after the piece of code that contains: class Body ; and the array bodies = [], btw

Runs on Python, and I do not plan to make it run in any other coding language.
