# N-Body-2D-Simulator
A simple N-Body simulator made by @Dahzito (me) during my free time.
Uses Newtonian Mechanics, and runs on the Verlet Integration method to compute the velocity and position of every body that is added to the program.

To this program, I've also added, automatic star's luminosity and temperature calculation using the Mass-Luminosity relation and the Stefan-Boltzmann's law, force felt by every body in the N-Body system, a simple mass loss calculation for planets that pass a certain temperature threshold, flux received by every body in the system, and then calculation of the temperature of every planet in the N-Body system using also the Stefan-Boltzmann's law and the planet's albedo, and gravity at the surface of every body to eventually see the Roche limit.

Stars also have on the simulation an automatic RGB value calculation depending on their temperature, and planets need a manual input for their RGBA values.

To add a new body to the program, add the following line of code, and change whatever you want:

bodies.append(Body(1.99*10**30, [0,0], [0, 0], [0, 0], "Star", 6.96e8, 0, 0.0, 1000, 0, "Star", [0.5, 0.5, 0.5, 1])) 
#The RGB values of the Star are automatically determined depending on their surface temperature, also calculated using the Mass-Luminosity Relation, and Stefan-Boltzmann's law.

bodies.append(Body(6.54*10**14, [0,0], [0, 22580], [2.5e10, 0], "Comet 1", 3.5e3, 280, 0.7, 100, 0, "Asteroid", [0, 0.5007843137, 0.6521568627, 0.8]))

bodies.append(Body(1.58*10**27, [0,0], [0, 27800], [5.5e10, 0], "Planet 2", 3.5e6, 280, 0.4, 900, 0, "Planet", [0.768627451, 0.4980392157, 0.3176470588, 0.8]))
#Mass, acceleration, velocity, coordinates, name, radii, temperature, albedo, temp_threshold, flux, type, RGB colors;

Or simply change the parameters of the already existing bodies.

This needs to be added after the piece of code that contains: class Body ; and the array bodies = [], btw

Runs on Python, and I do not plan to make it run in any other coding language.
