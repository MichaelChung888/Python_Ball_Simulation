# %%
# from ALL_CODE import *

from functions import *

from Ball_class import Ball

from Simulation_class import *


# Testing coding

"""
This is the testing file. It will check the diffrerent functions from my
class and should reproduce the graphs handed in on my report. Also shows some
of my debugging and checking.

Run %matplotlib auto in the console for the animations and plots to show up

Run each cell in order.
"""

# %% finds postion, velocity

c1 = Ball(1, 1, 2, 4, 3, 3)
c2 = Ball(2, 1, 4, 2, 0, 1)

c3 = Ball(1, 1, 2, 4, 3, 1)
c4 = Ball(2, 1, 4, 2, 1, 5)

c5 = Ball(0, 20, 0, 0, 0, 0, True)
c6 = Ball(1, 1, 0, 0, 1, 0)

c7 = Ball(1, 1, 3, 4, 1, 0)
c8 = Ball(0, 20, 3, 20, 0, 0, True)

c1.pos()  # returns position of c1
c2.vel()  # returns velocity of c2

# %% finds time to collision and finds speed

c1 = Ball(1, 1, 2, 4, 3, 3)
c2 = Ball(2, 1, 4, 2, 0, 1)

c3 = Ball(1, 1, 2, 4, 3, 1)
c4 = Ball(2, 1, 4, 2, 1, 5)

c5 = Ball(0, 20, 0, 0, 0, 0, True)  # finds speed after -container collision
c6 = Ball(1, 1, 0, 0, 1, 0)

c7 = Ball(1, 1, 3, 4, 1, 0)
c8 = Ball(0, 20, 3, 20, 0, 0, True)

print("time to collision between c1 and c2 is:",
      c1.time_to_collision(c2), 's and speed after is', c1.collide(c2), "m/s")
print("time to collision between c3 and c4 is:",
      c3.time_to_collision(c4), 's and speed after is', c3.collide(c4), "m/s")
print("time to collision between c5 and c6 is:",
      c5.time_to_collision(c6), 's and speed after is', c5.collide(c6), "m/s")
print("time to collision between c7 and c8 is:",
      c7.time_to_collision(c8), 's and spoeed after is', c7.collide(c8), "m/s")

# %% moves the ball to a new position afte time dt

print(c1.move(2))  # new position of c1 after 2 seconds
print(c1.move(5))

# %% returns methods in Simulation class (inherited from ball class)

help(Simulation)

# %% Runs the simulaiton for N balls and returns properties

"""
Can only see the values of each properties in the system after each
collision (itartaion) so ran the simulation.
Make sure to run %matplotlib auto in console before running this cell
"""

s1 = Simulation(5)

s1.run(100, animate=True)  # 5 balls after 100 collisions

# %% - extended to a simulation of 50 balls

s2 = Simulation(50)
s2.run(100, animate=True)

# %% testing raise error

s2.run(1000, animate=True, histograms=True)

# %% histrogram plots of ball distance and interball separation

s2.run(1000, histograms=True)

# %% KE time graph

s2.KE_time_graph(100)  # 100 collisions

# %%

s2.momentum_time_graph(100)

# %% mass time graph

s2.mass_time_graph(100)

# %% Pressure temperature graph

s2.pressure_temp_graph(1000)

# %% Pressure temperature graph od different radiis


a1 = Simulation(10, r_ball=1)
a2 = Simulation(10, r_ball=2)
a3 = Simulation(10, r_ball=0.1)

a1.pressure_temp_graph2(1000, a2, a3)

# %% Maxwell - Boltzmann fitted to theoretical speeds

c1 = Simulation(150, r_container=200)
c1.Maxwell_Boltmann()
