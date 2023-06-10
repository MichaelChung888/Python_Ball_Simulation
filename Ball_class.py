#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 17:24:28 2022

"""

# Thermodynamic snooker - usings PEP 8 to increase readability and consistency

import numpy as np


import pylab as pl


class Ball():

    """Class for initialising the different parameters of a ball in 2D"""

    def __init__(self, mass, radius, posx, posy, velx, vely,
                 container=False, colour="r"):
        """
        Ball Class formed of these attributes.

        Parameter:
            mass (float): The mass of the ball
            radius (float): The radius of the ball
            posx (float): The x component of the position of the ball
            posy (float): The y component of the position of the ball
            velx (float): The x component of the velocity of the ball
            vely (float): The y component of the velocity of the ball
            Container (Boolean): Describes weather the ball in desription is
                                a particle or the contianer
            Colour (string): The colour of the ball
        """
        self._mass = mass
        self._radius = radius
        self._position = np.array([posx, posy])
        self._velocity = np.array([velx, vely])
        self._container = container
        self._colour = colour
        # created a colour attribute used for checking and debugging
        if container is True:
            self._patch = pl.Circle(
                [posx, posy], radius, fill=False, ls="solid")
        else:
            self._patch = pl.Circle([posx, posy], radius, fc=colour)

    def pos(self):
        """Returns the position of the Ball"""

        print('position is', self._position)

    def set_pos(self, new):
        """Sets the new position of the Ball"""

        self._position = new

    def vel(self):
        """Returns the velocity of the Ball"""

        print('velocity is', self._velocity)

    def set_vel(self, new):
        """Sets the new velocity of the Ball"""

        self._velocity = new

    def get_patch(self):
        """Patch to be used later in the code"""

        return self._patch

    def move(self, dt):
        """Moves the ball to the new position that it travels to in time dt"""

        r_new = self._position + (self._velocity * dt)
        return r_new

    def time_to_collision(self, other):
        """Calculates the time to collision  between two balls"""

        r = self._position - other._position
        v = self._velocity - other._velocity
        r_dot = np.dot(r, r)
        v_dot = np.dot(v, v)
        r_dot_v = np.dot(r, v)
        if self._container is True or other._container is True:
            delta_radius = self._radius - other._radius
        else:  # has to be ball on ball
            delta_radius = self._radius + other._radius

        discrim = r_dot_v**2 - v_dot*(r_dot - delta_radius**2)

        if discrim <= 0:
            return 1000  # make very large time, so always larger than other t
        elif self._container is True or other._container is True:
            return (-r_dot_v + np.sqrt(discrim))/v_dot
        else:  # has to be ball on ball
            x = (-r_dot_v - np.sqrt(discrim))/v_dot
            if x <= 0:
                return 1000
            else:
                return x
            # want the smaller of the 2 values

    def collide(self, other):
        """Collide two balls together, finding its velocity after"""
        if(self.time_to_collision(other) <= 0):
            pass
        else:
            r1 = self._position - other._position
            r2 = other._position - self._position
            v1 = self._velocity - other._velocity
            v2 = other._velocity - self._velocity
            M = self._mass + other._mass
            r1_dot = np.dot(r1, r1)
            r2_dot = np.dot(r2, r2)
            r1_dot_v1 = np.dot(r1, v1)
            r2_dot_v2 = np.dot(r2, v2)
            if(r1_dot == 0):
                r1_unit = 0
                r2_unit = r2/r2_dot  # the unit vector
            elif(r2_dot == 0):
                r1_unit = r1/r1_dot
                r2_unit = 0  # the unit vector
            else:
                r1_unit = r1/r1_dot  # the unit vector
                r2_unit = r2/r2_dot  # the unit vector

            if(self._container is True):
                return -(other._velocity)  # Speed change is -u if with wall
            elif(other._container is True):
                return -(self._velocity)
            else:
                v1 = self._velocity - (2*self._mass/M)*(r1_dot_v1)*(r1_unit)
                v2 = other._velocity - (2*other._mass/M)*(r2_dot_v2)*(r2_unit)

                return v1, v2
