#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 17:12:44 2022

"""

import numpy as np

import random as rd


"""Functions to be imported and used in the Simulation and Ball class."""


def random_pos(r_container):
    """A function to generate a random position inside a cirlce
        using polar coordinates"""
    angle = rd.uniform(0, 2*(np.pi))  # random angle between 0 and 2pi
    xpos = rd.uniform(0,  r_container*np.cos(angle))
    ypos = rd.uniform(0,  r_container*np.sin(angle))
    return [xpos, ypos]


def distance_between_balls(all_balls, single_ball, r_ball):
    """A function to test if the balls are overlapping when initialised"""
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
    x = rd.randrange(0, 2) #Choose a random value between 0 and 1
    if x == 0:
        return rand1
    else:
        return rand2
