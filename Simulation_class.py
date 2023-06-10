#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 17:30:52 2022

@author: sharonliu
"""
from Ball_class import Ball

import numpy as np

import pylab as pl

import random as rd

import matplotlib.pyplot as plt


def random_pos(r_container):
    """A function to generate a random position inside a cirlce
        using polar coordinates"""
    angle = rd.uniform(0, 2*(np.pi))  # random angle between 0 and 2pi
    xpos = rd.uniform(0,  r_container*np.cos(angle))
    ypos = rd.uniform(0,  r_container*np.sin(angle))
    return [xpos, ypos]


def distance_between_balls(all_balls, single_ball, r_ball):
    """A function to test if the balls are being initialised overlapping"""
    for i in range(len(all_balls)):
        test = np.sqrt((all_balls[i]._position[0] - single_ball[0])
                       ** 2 + (all_balls[i]._position[1] - single_ball[1])**2)
        if test <= (2*r_ball):
            return True
    return False


def pressure_array(sim, num_frames):
    """A function to return pressure and temperature as an array"""
    p = []
    temp = []
    for frame in range(num_frames):
        print("iteration: ", frame+1)
        sim.next_collision()
        p.append(sim.pressure_check())
        temp.append(sim.temperature())
        sim.change_speed()
    return p, temp


def random_speed():
    """A function to set the speed at random following a normal distribution"""
    rand1 = 0
    rand2 = 0
    mean = 0.6
    std = 0.25
    while rand1 <= 0:
        rand1 = np.random.normal(mean, std)
    while rand2 >= 0:
        rand2 = np.random.normal(-mean, std)
    x = rd.randrange(0, 2)
    if x == 0:
        return rand1
    else:
        return rand2


class Simulation(Ball):

    """This is the inheritanced class Simulation from the ParentClass Ball"""
    # run help(Simulation) to check all the inherited methods

    def __init__(self, num, r_container=12, r_ball=1):
        """
        Simulation Class formed of these attributes.

        Parameter:
            num (int): The number of balls in the collision
            r_container (int): The rad
            ius of the container
            r_ball (int): The radius of the ball
        """

        balls_array = []
        c = Ball(0, r_container, 0, 0, 0, 0, True)
        """mass of container should really be infinity by code adapted to
            account for the infinite mass)"""
        for i in range(num):
            randvelx = random_speed()
            randvely = random_speed()

            position = random_pos(r_container-2)

            if (i == 0):
                b = Ball(1, r_ball, position[0], position[1], randvelx, randvely, False)
            else:
                while(distance_between_balls(balls_array, position, r_ball)):
                    position = random_pos(r_container-(2*r_ball))
                b = Ball(1, r_ball, position[0], position[1], randvelx, randvely, False)

            balls_array.append(b)
            # positioned the balls so that they are systematically spreaded
        balls_array.append(c)

        balls_array.append(c)
        self._number = num
        self._b = balls_array
        self._c_r = r_container
        self._t_elapsed = 0
        self._momentum_change = 0

    def set_time_elapsed(self, t):
        """A function to keep record of the time elapsed so far"""
        self._t_elapsed = t

    def set_momentum_change(self, m):
        """A function to keep record of the momentum change"""
        self._momentum_change = m

    def next_collision(self):
        """A function to upates the time, speed, distance and momentum of the
        balls after one collision"""
        x = []
        for i in range(len(self._b)):
            for j in range(i+1, len(self._b)):
                t = self._b[i].time_to_collision(self._b[j])
                if (not x) or (t < x[0]):
                    x = [t, i, j]  # returns time and the balls i and j

        # find the speed
        speed = self._b[x[1]].collide(self._b[x[2]])

        # Sets the current time elapsed across current simulation frame
        self.set_time_elapsed(self._t_elapsed + x[0])

        # Checks whether ball 1/ 2 are containers and adds momentum change
        if self._b[x[1]]._container is True:
            m = self._b[x[2]]._mass
            v = self._b[x[2]]._velocity
            abs_v = np.sqrt(np.dot(v, v))
            self.set_momentum_change(self._momentum_change + 2*m*abs_v)

        if self._b[x[2]]._container is True:
            m = self._b[x[1]]._mass
            v = self._b[x[1]]._velocity
            abs_v = np.sqrt(np.dot(v, v))
            self.set_momentum_change(self._momentum_change + 2*m*abs_v)

        # updates distances of balls
        for i in range(len(self._b)):
            if self._b[i]._container is True:
                pass
            else:
                distance = self._b[i].move(x[0])
                self._b[i].set_pos(distance)
                self._b[i]._patch.center = distance

        # speeds get updated down here
        if self._b[x[1]]._container is True:
            self._b[x[2]].set_vel(speed)

        elif self._b[x[2]]._container is True:
            self._b[x[1]].set_vel(speed)

        else:
            self._b[x[1]].set_vel(speed[0])
            self._b[x[2]].set_vel(speed[1])

    def kinetic_energy_check(self):
        """A function to calculate the kinetic energy of the system"""
        Sum_e = 0
        for i in range(len(self._b)):
            ke = 0.5 * self._b[i]._mass * \
                (np.dot(self._b[i]._velocity, self._b[i]._velocity))
            Sum_e += ke
        Ke = round(Sum_e, 1)
        print('The kinetic energy of the system is', Ke, 'J')
        return Ke

    def momentum_check(self):
        """A frunction to calculate the momentum of the system"""
        Sum_p = 0

        for i in range(len(self._b)):
            p = self._b[i]._mass * \
                np.sqrt(np.dot(self._b[i]._velocity, self._b[i]._velocity))
            Sum_p += p
        P = round(Sum_p)
        print('the momentum of the balls in the system is:', P, 'Kgm/s')
        return P

    def mass_check(self):
        """The mass of the balls in system to check for mass conservation"""
        Sum_mass = 0
        for i in range(len(self._b)):
            mass = self._b[i]._mass
            Sum_mass += mass
        print('The mass of the system is', Sum_mass, 'kg')
        return Sum_mass

    def temperature(self):
        """A function to find the temperature of the system"""
        k_b = 1.38e-23  # boltzman's constant
        T = (self.kinetic_energy_check()) / (self._number * k_b)
        print('Temperature is', T)
        return T

    def pressure_check(self):
        """A function to find the pressure exerted by the ball on the wall."""
        total_p = 0
        if self._momentum_change == 0:
            return 0
        else:
            total_p += self._momentum_change
            pressure = total_p / (self._t_elapsed * 4 * np.pi * self._c_r)
            print("Pressure is:", pressure, 'Pa')
            return pressure

    def hist_dist_center(self, dist_from_center):
        """A function to finds the balls distance from container center"""
        # center = [0, 0]
        for i in range(len(self._b)):
            d = np.sqrt(np.dot(self._b[i]._position, self._b[i]._position))
            dist_from_center.append(d)

    def hist_dist_balls(self, dist_balls):
        """A function to find interball separation distance of all the balls"""
        for i in range(len(self._b)):
            for j in range(i+1, len(self._b)):
                b1 = self._b[i]._position
                b2 = self._b[j]._position
                d = np.sqrt((b1[0] - b2[0])**2 + (b1[1] - b2[1])**2)
                dist_balls.append(d)

    def run(self, num_frames, animate=False, histograms=False,
            change_velocity=False):
        """
        A function to run the collisions and show its properties.

        Parameter:
            num_frame (int): how many frames of collsiisons to run
            animate (Boolean): run %matplotlib auto in console for animation
                                to play. Also return kinetic energy, mass,
                                momentum and pressure on wall of the system
                                after each collision.
            histograms (Boolean): Return histograms of the distance the balls
                                   are from the center and one of the
                                   interball separation distance.
        """

        if animate is True and histograms is True:
            raise Exception(
                'only one of the attributes animate or pressure and histograms'
                ' can be True')

        if animate:
            pl.figure()
            ax = pl.axes(xlim=(-self._c_r, self._c_r),
                         ylim=((-self._c_r, self._c_r)))
            for i in range(len(self._b)):
                if self._b[i]._container is True:
                    ax.add_artist(self._b[i].get_patch())
                ax.add_patch(self._b[i].get_patch())
        else:
            dist_from_center = []
            dist_balls = []

        for frame in range(num_frames):
            print("iteration: ", frame+1)
            if animate:
                self.next_collision()
                #self.kinetic_energy_check()
                #self.mass_check()
                #self.momentum_check()
                #self.pressure_check()
                #self.temperature()
                if(change_velocity):
                    self.change_speed()
            else:
                self.hist_dist_center(dist_from_center)
                self.hist_dist_balls(dist_balls)
            if animate:
                pl.pause(0.05)
        if animate:
            pl.show()
        else:
            plt.subplot(1, 2, 1)  # plots the 2 histograms side by side
            plt.title(
                'distance of 100 balls from the center after 1000 collisions')
            plt.xlabel('Distance')
            plt.ylabel('Frequency')
            plt.hist(dist_from_center, bins=5)

            plt.subplot(1, 2, 2)
            plt.title(
                'inter-ball separation of 100 balls after 1000 collsions')
            plt.xlabel('distance')
            plt.ylabel('Frequency')
            plt.hist(dist_balls, bins=10)

    def change_speed(self):
        """A function to increase the speed of each ball after each collision
            so that we can investigate pressure ve temperature relation"""

        for i in range(len(self._b)):
            if(self._b[i]._container is True):
                pass
            else:
                if(self._b[i]._velocity[0] < 0):
                    x_new_vel = -0.2
                else:
                    x_new_vel = 0.2
                if(self._b[i]._velocity[1] < 0):
                    y_new_vel = -0.2
                else:
                    y_new_vel = 0.2
                speed_incr = np.array([x_new_vel, y_new_vel])
                new_speed = self._b[i]._velocity + speed_incr
                self._b[i].set_vel(new_speed)

    def KE_time_graph(self, num_frames):
        """Generates a graph of time over kinetic energy of the system
           for a given number of frames"""
        ke = []  # y-axis is k
        for frame in range(num_frames):
            print("iteration: ", frame+1)
            self.next_collision()
            ke.append(self.kinetic_energy_check())
        # print('check the ke array:', ke)

        time = []  # x axis is time
        for frame in range(num_frames):
            self.next_collision()
            time.append(self._t_elapsed)

        plt.plot(time, ke)
        plt.title('Time Against System Kinetic Energy-'
                  '100 balls after 100 collisions ')
        plt.xlabel('Time (s)')
        plt.ylabel('Kinetic energy (J)')
        plt.show()

    def mass_time_graph(self, num_frames):
        """"Generates a graph of mass over time"""
        mass = []
        for frame in range(num_frames):
            print("iteration: ", frame+1)
            self.next_collision()
            mass.append(self.mass_check())

        time = []
        for frame in range(num_frames):
            self.next_collision()
            time.append(self._t_elapsed)

        plt.plot(time, mass)
        plt.title('Time Against System Mass-'
                  '100 balls after 100 collisions ')
        plt.xlabel('Time (s)')
        plt.ylabel('Mass (kg)')
        plt.show()

    def momentum_time_graph(self, num_frames):
        """Generates a graph of time over kinetic energy of the system
           for a given number of frames"""
        momentum = []
        for frame in range(num_frames):
            print("iteration: ", frame+1)
            self.next_collision()
            momentum.append(self.momentum_check())
        # print('check the momentum array:', momentum)

        time = []
        for frame in range(num_frames):
            self.next_collision()
            time.append(self._t_elapsed)

        plt.plot(time, momentum)
        plt.title('Time Against System Momentum-'
                  '100 balls after 100 collisions ')
        plt.xlabel('Time (s)')
        plt.ylim(0, 1000)
        plt.ylabel('Momentum (kgm/s)')
        plt.show()

    def pressure_temp_graph(self, num_frames):
        """A function to plot a temperature pressure graph for a number of
        balls with its velocity increasing in each frame"""
        p, temp = pressure_array(self, num_frames)
        # print('check the pressure array:', pressure,
        # 'check the temp array:', temp)

        pressure = np.array(p)
        m, b = np.polyfit(pressure, temp, 1)  # slope of m with b intercept
        plt.plot(pressure, temp, 'x', label=(
            "gradient is ", '{0:1.2E}'.format(m)))
        plt.plot(pressure, m*pressure + b)
        plt.title('Pressure Against Temperature-'
                  '50 balls after 1000 collisions ')
        gradient = ((4/3) * np.pi * (self._c_r)**3) / (self._number * 1.38e-23)
        print('theroertical gradient is', gradient)
        plt.xlabel('Pressure (Pa)')
        plt.ylabel('Temperature (K)')
        plt.legend()
        plt.show()

    def pressure_temp_graph2(self, num_frames, sim2, sim3):
        """A function to plot a temperature pressure graph for a number of
        balls with different radiis with its velocity increasing
        in each frame"""

        p1, temp1 = pressure_array(self, num_frames)
        p2, temp2 = pressure_array(sim2, num_frames)
        p3, temp3 = pressure_array(sim3, num_frames)

        pressure1 = np.array(p1)
        pressure2 = np.array(p2)
        pressure3 = np.array(p3)
        m, b = np.polyfit(pressure1, temp1, 1)
        x, y = np.polyfit(pressure2, temp2, 1)
        c, d = np.polyfit(pressure3, temp3, 1)
        plt.plot(pressure1, temp1, 'x', label='radius 1')
        plt.plot(pressure2, temp2, 'x', label='radius 2')
        plt.plot(pressure3, temp3, 'x', label='radius 0.1')
        plt.plot(pressure1, m*pressure1 + b)
        plt.plot(pressure2, x*pressure2 + y)
        plt.plot(pressure3, c*pressure3 + d)
        plt.title('Pressure Against Temperature-'
                  '10 balls after 1000 collisions ')
        plt.xlabel('Pressure (Pa)')
        plt.ylabel('Temperature (K)')
        plt.legend()
        plt.show()

    def Maxwell_Boltmann(self, m=1):
        """Maxwell-Boltzmann speed distribution graph with experimental data"""
        x = []
        for i in range(len(self._b)):
            x.append(np.dot(self._b[i]._velocity, self._b[i]._velocity))

        k_b = 1.38e-23
        T = self.temperature()
        v = np.linspace(0, 25, 10000)
        prob = m*np.exp(-m*v**2/(2*k_b*T))/(2*np.pi*T*k_b)*2*np.pi*v
        plt.plot(v, prob, label="Maxwell–Boltzmann distribution")
        plt.legend(loc="upper right")
        plt.hist(x, bins=15, density=True)
        plt.ylabel('Probability')
        plt.xlabel('Velocity of the balls (m/s)')
        plt.xlim(0, 5)
        plt.show()
