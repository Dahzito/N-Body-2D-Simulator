from math import*
from sys import*
from os import*
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

δt:float = 3 #days
δt_default:float = 2 #default timestep for normal conditions
Δt:float = 6000.0 #days

ε:float = 0 #Gravitational softening
force_temp:float = 0 #Temporary variable for storing gravitational force

G = 6.67430 * 10**-11 #Gravitational constant
σ = 5.67037 * 10**-8  #Stefan-Boltzmann constant

_t = 24*3600 #Conversion factor from days to seconds

MAX_TRAIL = 1000

class Body:
        def __init__(self, mass, acc, vel, coord, name, radii, temperature, albedo, temp_threshold, flux, type):
            
            #Dynamic Properties
            self.acceleration = list(acc)   # List as [acc_x, acc_y]
            self.velocity = list(vel)       # List as [vel_x, vel_y]
            self.coordinates = list(coord)  # List as [coord_x, coord_y]
            self.prev_acceleration = list(acc)  # Store previous acceleration for velocity Verlet

            #Physical Properties 
            self.name = name
            self.radii = radii
            self.mass = mass

            self.temperature = temperature  # Kelvin
            self.temp_threshold = temp_threshold  # Arbitrary threshold for temperature-based effects (e.g., radiation pressure)

            self.gravity = G * self.mass / self.radii**2  # Surface gravity

            self.luminosity = σ * 4 * pi * self.radii**2 * self.temperature**4
            self.albedo = albedo  # Albedo for energy absorption calculations
            self.flux = flux      # Flux received from other bodies (for energy balance calculations)

            #Other properties
            self.type = type  # Type of body (e.g., "Star", "Planet", "Asteroid")

            self.status = "Stable"  # Status can be "Stable", "Torn Apart", for Roche's limit effects

bodies = []

bodies.append(Body(1.99*10**30, [0,0], [0, -25080], [-1.5e11, 0], "Star 1", 6.96e8, 5878, 0.0, 1000, 0, "Star"))

bodies.append(Body(5.97*10**30, [0,0], [0, 22080], [1.5e11, 0], "Star 2", 6.37e6, 6578, 0.0, 1000, 0, "Star"))

bodies.append(Body(1.898*10**28, [0,0], [0, 73270], [3.0e11, 0], "Star 3", 7.15e7, 4578, 0.0, 1000, 0, "Star"))

bodies.append(Body(2.54*10**24, [0,0], [0, 50080], [2.5e11, 0], "Planet 1", 6.3e6, 280, 0.15, 250, 0, "Planet"))
#Mass, acceleration, velocity, coordinates, name, radii, temperature, albedo, temp_threshold, flux, type;

def Update_():
    global δt, ε, force_temp
    radius = []
    dt_seconds = δt * _t  # Convert to seconds

    force_temp_dx = 0; force_temp_dy = 0;
    
    ε = min(b.radii for b in bodies) * 1 # Softening length based on smallest radius

    #Calculate the distance of every object on the class Body
    for b in range(len(bodies)):
        bodies[b].acceleration = [0 , 0]
        bodies[b].flux = 0
        
    for i in range(len(bodies)):
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].coordinates[0] - bodies[i].coordinates[0]
                dy = bodies[j].coordinates[1] - bodies[i].coordinates[1]

                r = sqrt(dx**2 + dy**2)
                radius.append(r)

                r = sqrt(dx**2 + dy**2 + ε**2)
                bodies[i].acceleration[0] += (G * bodies[j].mass * dx) / r**3
                bodies[i].acceleration[1] += (G * bodies[j].mass * dy) / r**3



                # Accumulate flux from body j
                r_true = sqrt(dx**2 + dy**2)
                bodies[i].flux += bodies[j].luminosity / (4 * pi * r_true**2)

                force_temp_dx += G*bodies[i].mass*bodies[j].mass*dx/r**3
                force_temp_dy += G*bodies[i].mass*bodies[j].mass*dy/r**3
                force_temp = sqrt(force_temp_dx**2 + force_temp_dy**2)  # Magnitude of the total force

                if force_temp > bodies[i].gravity * bodies[i].mass:
                    bodies[i].status = "Torn Apart"



    #adaptative timestep based on minimum distance to improve accuracy during close encounters
    if min(radius) < 5e8: δt = 0.005
    if min(radius) < 5e9: δt = 0.05
    elif min(radius) < 4e10: δt = 0.25
    elif min(radius) < 2e11: δt = 0.5
    elif min(radius) >  3e11: δt = δt_default 
    
    dt_seconds = δt * _t

    # Calculate temperature and mass loss after all flux has been accumulated
    for i in range(len(bodies)):
        temp_new = (bodies[i].flux * ((1 - bodies[i].albedo) / σ))**0.25
        
        k0 = 1e5; T0 = 1000
        k = k0 * (temp_new / T0)**2
        
        if bodies[i].type == "Planet":
            bodies[i].temperature = temp_new
            if temp_new > bodies[i].temp_threshold:
                bodies[i].mass -= k * dt_seconds
        
        elif bodies[i].type == "Star":
            bodies[i].luminosity = 3.828*10**26*((bodies[i].mass/(1.989*10**30))**3.5)
            bodies[i].temperature = (bodies[i].luminosity / (4 * pi * bodies[i].radii**2 * σ))**0.25
        
        else:
            if temp_new > bodies[i].temperature:
                bodies[i].temperature = temp_new

        

    # Velocity Verlet integration: more accurate and better energy conservation
    # Step 1: Update positions: x_new = x + v*dt + 0.5*a*dt²
    for n in range(len(bodies)):
        b = bodies[n]
        b.coordinates[0] += b.velocity[0] * dt_seconds + 0.5 * b.acceleration[0] * dt_seconds**2
        b.coordinates[1] += b.velocity[1] * dt_seconds + 0.5 * b.acceleration[1] * dt_seconds**2
    
    # Step 2: Compute new accelerations at new positions
    # (already done at start of next Update_() call, but we need to do it once here first)
    # We'll store old acceleration and compute new one
    for b in range(len(bodies)):
        bodies[b].prev_acceleration = bodies[b].acceleration.copy()
    
    for b in range(len(bodies)):
        bodies[b].acceleration = [0, 0]
    
    for i in range(len(bodies)):
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].coordinates[0] - bodies[i].coordinates[0]
                dy = bodies[j].coordinates[1] - bodies[i].coordinates[1]

                r = sqrt(dx**2 + dy**2 + ε**2)

                bodies[i].acceleration[0] += (G * bodies[j].mass * dx) / r**3
                bodies[i].acceleration[1] += (G * bodies[j].mass * dy) / r**3

    # Step 3: Update velocities using average acceleration: v_new = v + 0.5*(a_old + a_new)*dt
    for n in range(len(bodies)):
        b = bodies[n]
        b.velocity[0] += 0.5 * (b.prev_acceleration[0] + b.acceleration[0]) * dt_seconds
        b.velocity[1] += 0.5 * (b.prev_acceleration[1] + b.acceleration[1]) * dt_seconds


#------ Visualization of the N Body system using Matplotlib ------

def calculate_energy():
    energies = [0] 
    #First value is for the total energy of the system, then each body will have its own energy value

    for i in range(len(bodies)):
        b = bodies[i]
        kinetic_energy = 0.5 * b.mass * (b.velocity[0]**2 + b.velocity[1]**2)

        # Potential energy: -G * m1 * m2 / r
        potential_energy = 0
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].coordinates[0] - b.coordinates[0]
                dy = bodies[j].coordinates[1] - b.coordinates[1]
                r = sqrt(dx**2 + dy**2 + ε**2)
                potential_energy -= 0.5 * (G * b.mass * bodies[j].mass) / r

        # Total energy
        total_energy = kinetic_energy + potential_energy
        energies[0] += total_energy
        energies.append(total_energy)

    return energies

def calculate_forces():
    forces = [] 
    #Each body will have its own force value

    for i in range(len(bodies)):
        force_temp_dx = 0
        force_temp_dy = 0
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].coordinates[0] - bodies[i].coordinates[0]
                dy = bodies[j].coordinates[1] - bodies[i].coordinates[1]
                r = sqrt(dx**2 + dy**2 + ε**2)

                force_temp_dx += G*bodies[i].mass*bodies[j].mass*dx/r**3
                force_temp_dy += G*bodies[i].mass*bodies[j].mass*dy/r**3

        # Total energy
        total_force = sqrt(force_temp_dx**2 + force_temp_dy**2)  # Magnitude of the total force
        forces.append(total_force)

    return forces



def Visualize(energy_reset_period=6000):
    """Visualize the n-body system with energy tracking.
    
    Args:
        energy_reset_period: Reset the energy graph every N days (default: 365)
        force_reset_period: Reset the force graph every N days (default: 365)
    """
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1], width_ratios=[1.5, 1, 1], hspace=0.35, wspace=0.35)
    ax = fig.add_subplot(gs[:, 0])  # Orbital plot spans all 3 rows on left
    ax_energy = fig.add_subplot(gs[0, 1])  # Energy plot
    ax_force = fig.add_subplot(gs[1, 1])   # Force plot
    ax_mass = fig.add_subplot(gs[0, 2])    # Mass plot
    ax_flux = fig.add_subplot(gs[1, 2])    # Flux plot
    ax_temperature = fig.add_subplot(gs[2, 2])  # Temperature plot
    ax.set_aspect('equal')

    # auto limits based on initial system
    max_r = max(max(abs(x.coordinates[0]), abs(x.coordinates[1])) for x in bodies)
    limit = max(max_r * 2, 1e10)

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("x coordinates (m)")
    ax.set_ylabel("y coordinates (m)")

    colors = ["#FFBB00", "#424242F8", "#CA9B00", "#00329E", "#C92F00", "#8B7200", "#FFE054", "#00A1FF"]

    points = []
    trails = []

    trails_x = [[] for _ in bodies]
    trails_y = [[] for _ in bodies]

    for i, b in enumerate(bodies):
        p, = ax.plot([], [], 'o', color=colors[i % len(colors)], markersize=6)
        t, = ax.plot([], [], color=colors[i % len(colors)], linewidth=1)

        points.append(p)
        trails.append(t)

    info_text = ax.text(
        0.02,
        0.98,
        '',
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='top',
        bbox=dict(facecolor='black', alpha=0.6, boxstyle='round,pad=0.4'),
        color='white'
    )

    # Energy plot setup
    initial_energies = calculate_energy()
    energy_lines = []
    energy_data = [[e] for e in initial_energies]
    time_data = [0.0]

    limit = max(abs(e) for e in initial_energies) * 1.5 if initial_energies else 1.0

    for i in range(len(initial_energies)):
        label = "Entire System" if i == 0 else bodies[i-1].name
        line, = ax_energy.plot(time_data, energy_data[i], label=label, color=colors[i-1 % len(colors)])
        energy_lines.append(line)

    ax_energy.set_xlim(0, energy_reset_period)
    ax_energy.set_ylim(-limit, limit)  # Adjust limits as needed
    ax_energy.set_xlabel("Time (days)")
    ax_energy.set_ylabel("Total Mechanical Energy (J)")
    ax_energy.legend()


    # Force plot setup
    initial_forces = calculate_forces()
    force_lines = []
    force_data = [[f] for f in initial_forces]
    time_data_f = [0.0]  # Separate time data for force plot

    limit = max(abs(f) for f in initial_forces) * 1.1 if initial_forces else 1.0

    for i in range(len(initial_forces)):
        label = bodies[i].name
        line, = ax_force.plot(time_data_f, force_data[i], label=label, color=colors[i % len(colors)])
        force_lines.append(line)

    ax_force.set_xlim(0, energy_reset_period)
    ax_force.set_ylim(0, limit)  # Adjust limits as needed
    ax_force.set_xlabel("Time (days)")
    ax_force.set_ylabel("Total Force felt by each body (N)")
    ax_force.legend()

    # Mass plot setup
    mass_lines = []
    mass_data = [[b.mass] for b in bodies]
    time_data_mass = [0.0]

    initial_masses = [b.mass for b in bodies]
    limit = max(initial_masses) * 1.5 if initial_masses else 1.0

    for i, b in enumerate(bodies):
        line, = ax_mass.plot(time_data_mass, mass_data[i], label=b.name, color=colors[i % len(colors)])
        mass_lines.append(line)

    ax_mass.set_xlim(0, energy_reset_period)
    ax_mass.set_ylim(0, limit)
    ax_mass.set_xlabel("Time (days)")
    ax_mass.set_ylabel("Mass (kg)")
    ax_mass.legend()

    # Flux plot setup
    flux_lines = []
    flux_data = [[b.flux] for b in bodies]
    time_data_flux = [0.0]

    initial_fluxes = [b.flux for b in bodies]
    max_flux = max(initial_fluxes) if initial_fluxes else 1.0
    limit = max(max_flux * 1.5, 1.0)  # Set a minimum limit to avoid zero range

    for i, b in enumerate(bodies):
        line, = ax_flux.plot(time_data_flux, flux_data[i], color=colors[i % len(colors)]) #, label=b.name
        flux_lines.append(line)

    ax_flux.set_xlim(0, energy_reset_period)
    ax_flux.set_ylim(0, limit)
    ax_flux.set_xlabel("Time (days)")
    ax_flux.set_ylabel("Flux received by each body (W/m²)")
    ax_flux.legend()

    # Temperature plot setup
    temperature_lines = []
    temperature_data = [[b.temperature] for b in bodies]
    time_data_temperature = [0.0]

    initial_temperatures = [b.temperature for b in bodies]
    limit = max(initial_temperatures) * 1.5 if initial_temperatures else 1.0

    for i, b in enumerate(bodies):
        line, = ax_temperature.plot(time_data_temperature, temperature_data[i], color=colors[i % len(colors)]) #, label=b.name
        temperature_lines.append(line)

    ax_temperature.set_xlim(0, energy_reset_period)
    ax_temperature.set_ylim(0, limit)
    ax_temperature.set_xlabel("Time (days)")
    ax_temperature.set_ylabel("Temperature of the body (K)")
    ax_temperature.legend()

    # Track loop iterations for continuous time display
    loop_count = [0]
    previous_frame = [0]
    last_reset_time_energy = [0]
    last_reset_time_force = [0]
    last_reset_time_mass = [0]
    last_reset_time_flux = [0]
    last_reset_time_temperature = [0]
    collisions = []  # Track collision events: (frame, time, body1, body2)

    def init():
        return points + trails + energy_lines + force_lines + mass_lines + flux_lines + temperature_lines

    def update(frame):
        # Detect when animation loops back to start
        if frame < previous_frame[0]:
            loop_count[0] += 1
        previous_frame[0] = frame

        Update_()

        for i, b in enumerate(bodies):

            trails_x[i].append(b.coordinates[0])
            trails_y[i].append(b.coordinates[1])

            if len(trails_x[i]) > MAX_TRAIL:
                trails_x[i].pop(0)
                trails_y[i].pop(0)

            points[i].set_data([b.coordinates[0]], [b.coordinates[1]])
            trails[i].set_data(trails_x[i], trails_y[i])

        frame_limit_reached_n = 0
        if frame*δt >= Δt - δt:
            frame_limit_reached_n += 1

        # Update energy data
        energies = calculate_energy()
        forces = calculate_forces()
        current_time = frame + loop_count[0] * Δt

        # Reset energy graph if it exceeds the reset period
        time_since_reset_energy = current_time - last_reset_time_energy[0]
        if time_since_reset_energy >= energy_reset_period:
            last_reset_time_energy[0] = current_time
            energy_data[:] = [[] for _ in energies]
            time_data[:] = []
            time_since_reset_energy = 0

        # Reset force graph if it exceeds the reset period
        time_since_reset_force = current_time - last_reset_time_force[0]
        if time_since_reset_force >= energy_reset_period:
            last_reset_time_force[0] = current_time
            force_data[:] = [[] for _ in forces]
            time_data_f[:] = []
            time_since_reset_force = 0

        time_data.append(time_since_reset_energy); time_data_f.append(time_since_reset_force)
        
        for energy_list, line, energy in zip(energy_data, energy_lines, energies):
            energy_list.append(energy)
            line.set_data(time_data, energy_list)
        
        for force_list, line, force in zip(force_data, force_lines, forces):
            force_list.append(force)
            line.set_data(time_data_f, force_list)

        # Handle mass, flux, and temperature resets and updates
        time_since_reset_mass = current_time - last_reset_time_mass[0]
        if time_since_reset_mass >= energy_reset_period:
            last_reset_time_mass[0] = current_time
            mass_data[:] = [[] for _ in bodies]
            time_data_mass[:] = []
            time_since_reset_mass = 0

        time_since_reset_flux = current_time - last_reset_time_flux[0]
        if time_since_reset_flux >= energy_reset_period:
            last_reset_time_flux[0] = current_time
            flux_data[:] = [[] for _ in bodies]
            time_data_flux[:] = []
            time_since_reset_flux = 0

        time_since_reset_temperature = current_time - last_reset_time_temperature[0]
        if time_since_reset_temperature >= energy_reset_period:
            last_reset_time_temperature[0] = current_time
            temperature_data[:] = [[] for _ in bodies]
            time_data_temperature[:] = []
            time_since_reset_temperature = 0

        time_data_mass.append(time_since_reset_mass)
        time_data_flux.append(time_since_reset_flux)
        time_data_temperature.append(time_since_reset_temperature)

        for mass_list, line, b in zip(mass_data, mass_lines, bodies):
            mass_list.append(b.mass)
            line.set_data(time_data_mass, mass_list)

        for flux_list, line, b in zip(flux_data, flux_lines, bodies):
            flux_list.append(b.flux)
            line.set_data(time_data_flux, flux_list)

        for temp_list, line, b in zip(temperature_data, temperature_lines, bodies):
            temp_list.append(b.temperature)
            line.set_data(time_data_temperature, temp_list)

        # Build info text with collision history
        info_lines = [f"Time passed: {current_time:.1f} days | {current_time / 365.25:.1f} years", f"Timestep: {δt} days"]
        
        # Show last 3 collisions
        recent_collisions = collisions[-3:] if collisions else []
        if recent_collisions:
            info_lines.append("\nRecent collisions:")
            for _, col_time, col_msg in recent_collisions:
                info_lines.append(f"  {col_msg} @ {col_time:.1f}d")
        
        info_text.set_text("\n".join(info_lines))

        ax_energy.set_xlim(0, energy_reset_period)
        ax_force.set_xlim(0, energy_reset_period)
        ax_mass.set_xlim(0, energy_reset_period)
        ax_flux.set_xlim(0, energy_reset_period)
        ax_temperature.set_xlim(0, energy_reset_period)

        # Dynamic y-limit updates based on current data        
        """
        if mass_data and any(mass_data):
            max_mass = max(max(m) for m in mass_data if m)
            ax_mass.set_ylim(0, max_mass * 1.2)
        """
        
        if flux_data and any(flux_data):
            max_flux = max(max(f) for f in flux_data if f)
            ax_flux.set_ylim(0, max_flux * 1.2)
        
        
        if temperature_data and any(temperature_data):
            max_temp = max(max(t) for t in temperature_data if t)
            ax_temperature.set_ylim(0, max_temp * 1.2)
        
        # Check for collisions between all body pairs
        for i in range(len(bodies)):
            for j in range(i + 1, len(bodies)):
                dx = bodies[j].coordinates[0] - bodies[i].coordinates[0]
                dy = bodies[j].coordinates[1] - bodies[i].coordinates[1]
                r = sqrt(dx**2 + dy**2)
                
                if r < (bodies[i].radii + bodies[j].radii):
                    collision_msg = f"Collision: {bodies[i].name} & {bodies[j].name}"
                    collisions.append((frame, current_time, collision_msg))
                    print(f"Collision detected between {bodies[i].name} and {bodies[j].name} at time {frame*δt:.2f} days!")
                    # Simple collision response: merge bodies (conservation of mass and momentum)
                    total_mass = bodies[i].mass + bodies[j].mass
                    new_velocity_x = (bodies[i].velocity[0] * bodies[i].mass + bodies[j].velocity[0] * bodies[j].mass) / total_mass
                    new_velocity_y = (bodies[i].velocity[1] * bodies[i].mass + bodies[j].velocity[1] * bodies[j].mass) / total_mass
                    new_coordinates_x = (bodies[i].coordinates[0] * bodies[i].mass + bodies[j].coordinates[0] * bodies[j].mass) / total_mass
                    new_coordinates_y = (bodies[i].coordinates[1] * bodies[i].mass + bodies[j].coordinates[1] * bodies[j].mass) / total_mass
                            
                    # Update body i to be the merged body, and remove body j
                    bodies[i].mass = total_mass
                    bodies[i].velocity = [new_velocity_x, new_velocity_y]
                    bodies[i].coordinates = [new_coordinates_x, new_coordinates_y]
                    bodies.pop(j)
                    # Exit inner loop to avoid index issues after pop
                    break
        
        return points + trails + energy_lines + force_lines + mass_lines + flux_lines + temperature_lines + [info_text]
    

    ani = FuncAnimation(
        fig,
        update,
        frames=np.arange(0, Δt, δt),
        init_func=init,
        interval=2,
        blit=True,
        repeat=False
    )

    plt.show()

# ---------------- RUN PROGRAM ----------------

Visualize()#Mass, acceleration, velocity, coordinates, name, radii;

def Update_():
    global δt, ε
    radius = []
    dt_seconds = δt * _t  # Convert to seconds
    
    ε = min(b.radii for b in bodies) * 1 # Softening length based on smallest radius

    #Calculate the distance of every object on the class Body
    for b in range(len(bodies)):
        bodies[b].acceleration = [0 , 0]
        
    for i in range(len(bodies)):
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].coordinates[0] - bodies[i].coordinates[0]
                dy = bodies[j].coordinates[1] - bodies[i].coordinates[1]

                r = sqrt(dx**2 + dy**2)
                radius.append(r)

                r = sqrt(dx**2 + dy**2 + ε**2)
                bodies[i].acceleration[0] += (G * bodies[j].mass * dx) / r**3
                bodies[i].acceleration[1] += (G * bodies[j].mass * dy) / r**3

    #adaptative timestep based on minimum distance to improve accuracy during close encounters
    if min(radius) < 5e8: δt = 0.005
    if min(radius) < 5e9: δt = 0.05
    elif min(radius) < 4e10: δt = 0.25
    elif min(radius) < 2e11: δt = 0.5
    elif min(radius) >  3e11: δt = δt_default 
    
    dt_seconds = δt * _t

    # Velocity Verlet integration: more accurate and better energy conservation
    # Step 1: Update positions: x_new = x + v*dt + 0.5*a*dt²
    for n in range(len(bodies)):
        b = bodies[n]
        b.coordinates[0] += b.velocity[0] * dt_seconds + 0.5 * b.acceleration[0] * dt_seconds**2
        b.coordinates[1] += b.velocity[1] * dt_seconds + 0.5 * b.acceleration[1] * dt_seconds**2
    
    # Step 2: Compute new accelerations at new positions
    # (already done at start of next Update_() call, but we need to do it once here first)
    # We'll store old acceleration and compute new one
    for b in range(len(bodies)):
        bodies[b].prev_acceleration = bodies[b].acceleration.copy()
    
    for b in range(len(bodies)):
        bodies[b].acceleration = [0, 0]
    
    for i in range(len(bodies)):
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].coordinates[0] - bodies[i].coordinates[0]
                dy = bodies[j].coordinates[1] - bodies[i].coordinates[1]

                r = sqrt(dx**2 + dy**2 + ε**2)

                bodies[i].acceleration[0] += (G * bodies[j].mass * dx) / r**3
                bodies[i].acceleration[1] += (G * bodies[j].mass * dy) / r**3

    # Step 3: Update velocities using average acceleration: v_new = v + 0.5*(a_old + a_new)*dt
    for n in range(len(bodies)):
        b = bodies[n]
        b.velocity[0] += 0.5 * (b.prev_acceleration[0] + b.acceleration[0]) * dt_seconds
        b.velocity[1] += 0.5 * (b.prev_acceleration[1] + b.acceleration[1]) * dt_seconds


#------ Visualization of the 3 Body system using Matplotlib ------

def calculate_energy():
    energies = [0] 
    #First value is for the total energy of the system, then each body will have its own energy value

    for i in range(len(bodies)):
        b = bodies[i]
        kinetic_energy = 0.5 * b.mass * (b.velocity[0]**2 + b.velocity[1]**2)

        # Potential energy: -G * m1 * m2 / r
        potential_energy = 0
        for j in range(len(bodies)):
            if i != j:
                dx = bodies[j].coordinates[0] - b.coordinates[0]
                dy = bodies[j].coordinates[1] - b.coordinates[1]
                r = sqrt(dx**2 + dy**2 + ε**2)
                potential_energy -= 0.5 * (G * b.mass * bodies[j].mass) / r

        # Total energy
        total_energy = kinetic_energy + potential_energy
        energies[0] += total_energy
        energies.append(total_energy)



    return energies


def Visualize(energy_reset_period=6000):
    """Visualize the n-body system with energy tracking.
    
    Args:
        energy_reset_period: Reset the energy graph every N days (default: 365)
    """
    fig, (ax, ax_energy) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [1, 1.2]})
    ax.set_aspect('equal')

    # auto limits based on initial system
    max_r = max(max(abs(x.coordinates[0]), abs(x.coordinates[1])) for x in bodies)
    limit = max(max_r * 2, 1e10)

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("x coordinates (m)")
    ax.set_ylabel("y coordinates (m)")

    colors = ["#FFBB00", "#424242F8", "#CA9B00", "#00329E", "#C92F00", "#8B7200", "#FFE054", "#00A1FF"]

    points = []
    trails = []

    trails_x = [[] for _ in bodies]
    trails_y = [[] for _ in bodies]

    for i, b in enumerate(bodies):
        p, = ax.plot([], [], 'o', color=colors[i % len(colors)], markersize=6)
        t, = ax.plot([], [], color=colors[i % len(colors)], linewidth=1)

        points.append(p)
        trails.append(t)

    info_text = ax.text(
        0.02,
        0.98,
        '',
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment='top',
        bbox=dict(facecolor='black', alpha=0.6, boxstyle='round,pad=0.4'),
        color='white'
    )

    # Energy plot setup
    initial_energies = calculate_energy()
    energy_lines = []
    energy_data = [[e] for e in initial_energies]
    time_data = [0.0]

    limit = max(abs(e) for e in initial_energies) * 1.5 if initial_energies else 1.0

    for i in range(len(initial_energies)):
        label = "Entire System" if i == 0 else bodies[i-1].name
        line, = ax_energy.plot(time_data, energy_data[i], label=label, color=colors[i-1 % len(colors)])
        energy_lines.append(line)

    ax_energy.set_xlim(0, energy_reset_period)
    ax_energy.set_ylim(-limit, limit)  # Adjust limits as needed
    ax_energy.set_xlabel("Time (days)")
    ax_energy.set_ylabel("Total Mechanical Energy (J)")
    ax_energy.legend()

    # Track loop iterations for continuous time display
    loop_count = [0]
    previous_frame = [0]
    last_reset_time = [0]
    collisions = []  # Track collision events: (frame, time, body1, body2)

    def init():
        return points + trails + energy_lines

    def update(frame):
        # Detect when animation loops back to start
        if frame < previous_frame[0]:
            loop_count[0] += 1
        previous_frame[0] = frame

        Update_()

        for i, b in enumerate(bodies):

            trails_x[i].append(b.coordinates[0])
            trails_y[i].append(b.coordinates[1])

            if len(trails_x[i]) > MAX_TRAIL:
                trails_x[i].pop(0)
                trails_y[i].pop(0)

            points[i].set_data([b.coordinates[0]], [b.coordinates[1]])
            trails[i].set_data(trails_x[i], trails_y[i])

        frame_limit_reached_n = 0
        if frame*δt >= Δt - δt:
            frame_limit_reached_n += 1

        # Update energy data
        energies = calculate_energy()
        current_time = frame + loop_count[0] * Δt

        # Reset energy graph if it exceeds the reset period
        time_since_reset = current_time - last_reset_time[0]
        if time_since_reset >= energy_reset_period:
            last_reset_time[0] = current_time
            energy_data[:] = [[] for _ in energies]
            time_data[:] = []
            time_since_reset = 0

        time_data.append(time_since_reset)
        for energy_list, line, energy in zip(energy_data, energy_lines, energies):
            energy_list.append(energy)
            line.set_data(time_data, energy_list)

        # Build info text with collision history
        info_lines = [f"Time passed: {current_time:.1f} days | {current_time / 365.25:.1f} years", f"Timestep: {δt} days"]
        
        # Show last 3 collisions
        recent_collisions = collisions[-3:] if collisions else []
        if recent_collisions:
            info_lines.append("\nRecent collisions:")
            for _, col_time, col_msg in recent_collisions:
                info_lines.append(f"  {col_msg} @ {col_time:.1f}d")
        
        info_text.set_text("\n".join(info_lines))

        ax_energy.set_xlim(0, energy_reset_period)

        # Check for collisions between all body pairs
        for i in range(len(bodies)):
            for j in range(i + 1, len(bodies)):
                dx = bodies[j].coordinates[0] - bodies[i].coordinates[0]
                dy = bodies[j].coordinates[1] - bodies[i].coordinates[1]
                r = sqrt(dx**2 + dy**2)
                
                if r < (bodies[i].radii + bodies[j].radii):
                    collision_msg = f"Collision: {bodies[i].name} & {bodies[j].name}"
                    collisions.append((frame, current_time, collision_msg))
                    print(f"Collision detected between {bodies[i].name} and {bodies[j].name} at time {frame*δt:.2f} days!")
                    # Simple collision response: merge bodies (conservation of mass and momentum)
                    total_mass = bodies[i].mass + bodies[j].mass
                    new_velocity_x = (bodies[i].velocity[0] * bodies[i].mass + bodies[j].velocity[0] * bodies[j].mass) / total_mass
                    new_velocity_y = (bodies[i].velocity[1] * bodies[i].mass + bodies[j].velocity[1] * bodies[j].mass) / total_mass
                    new_coordinates_x = (bodies[i].coordinates[0] * bodies[i].mass + bodies[j].coordinates[0] * bodies[j].mass) / total_mass
                    new_coordinates_y = (bodies[i].coordinates[1] * bodies[i].mass + bodies[j].coordinates[1] * bodies[j].mass) / total_mass
                            
                    # Update body i to be the merged body, and remove body j
                    bodies[i].mass = total_mass
                    bodies[i].velocity = [new_velocity_x, new_velocity_y]
                    bodies[i].coordinates = [new_coordinates_x, new_coordinates_y]
                    bodies.pop(j)
                    # Exit inner loop to avoid index issues after pop
                    break
        
        return points + trails + energy_lines + [info_text]
    

    ani = FuncAnimation(
        fig,
        update,
        frames=np.arange(0, Δt, δt),
        init_func=init,
        interval=2,
        blit=True,
        repeat=True
    )

    plt.show()

# ---------------- RUN PROGRAM ----------------

Visualize()
