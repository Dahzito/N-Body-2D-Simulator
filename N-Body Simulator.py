from math import*
from sys import*
from os import*
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

δt:float = 2 #days
δt_default:float = 1.5 #default timestep for normal conditions
Δt:float = 6000.0 #days

ε:float = 0 #Gravitational softening

G = 6.67430 * 10**-11

_t = 24*3600 #Conversion factor from days to seconds

MAX_TRAIL = 1000

class Body:
        def __init__(self, mass, acc, vel, coord, name, radii):
            self.mass = mass
            self.acceleration = list(acc)   # List as [acc_x, acc_y]
            self.velocity = list(vel)       # List as [vel_x, vel_y]
            self.coordinates = list(coord)  # List as [coord_x, coord_y]
            self.prev_acceleration = list(acc)  # Store previous acceleration for velocity Verlet
            self.name = name
            self.radii = radii

bodies = []

bodies.append(Body(1.99*10**30, [0,0], [0, -25080], [-1.5e11, 0], "Star 1", 6.96e8))

bodies.append(Body(5.97*10**30, [0,0], [0, 22080], [1.5e11, 0], "Star 2", 6.37e6))

bodies.append(Body(1.898*10**28, [0,0], [0, 73270], [3.0e11, 0], "Star 3", 7.15e7))
#Mass, acceleration, velocity, coordinates, name, radii;

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