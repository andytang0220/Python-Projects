#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 13:50:06 2020

@author: dimitrisgiovanis
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import numpy as np


def plot(pos):
    """
    A function to plot the graphene sheet or the carbon nanotube.
    Input:
       pos: A numpy array containing the coordinates of the atom included
       in the nanotube (pos.shape = (N, 3) or (N, 2)).
    """
    if pos.shape[1] == 3:

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        u = np.linspace(0, np.pi, 100)
        v = np.linspace(0, 2 * np.pi, 100)

        Xs = 0.1 * np.outer(np.cos(u), np.sin(v))
        Ys = 0.1 * np.outer(np.sin(u), np.sin(v))
        Zs = 0.1 * 2 * np.outer(np.ones(np.size(u)), np.cos(v))

        for t in range(pos.shape[0]):

            ax.plot_surface(Xs + pos[t, 0], Ys + pos[t, 1], Zs + pos[t, 2], rstride=4, cstride=4,
                            color='r');  # Plots the balls

            for t2 in range(t + 1, pos.shape[0]):
                r = sum((pos[t, :] - pos[t2, :]) ** 2);
                if r > 0.95 and r < 1.05:
                    print()
                    ax.plot([pos[t, 0], pos[t2, 0]], [pos[t, 1], pos[t2, 1]], [pos[t, 2], pos[t2, 2]], c='b')
                elif r < 0.05:
                    print('Error: Atoms overlapping!')
        plt.savefig('Nanotube.png')
    elif pos.shape[1] == 2:

        fig = plt.figure()
        # Plot Graphene
        u = np.linspace(0, np.pi, 100)
        v = np.linspace(0, 2 * np.pi, 100)
        Xs = 0.1 * np.outer(np.cos(u), np.sin(v))
        Ys = 0.1 * np.outer(np.sin(u), np.sin(v))

        for t in range(pos.shape[0]):

            plt.scatter(Xs + pos[t, 0], Ys + pos[t, 1], color='r');  # Plots the balls

            for t2 in range(t + 1, pos.shape[0]):
                r = sum((pos[t, :] - pos[t2, :]) ** 2)
                if r > 0.95 and r < 1.05:
                    print()
                    plt.plot([pos[t, 0], pos[t2, 0]], [pos[t, 1], pos[t2, 1]], c='b')
                elif r < 0.05:
                    print('Error: Atoms overlapping!')
        plt.savefig('Graphene.png')

    else:
        raise ValueError('Wrong dimensions')
