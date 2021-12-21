import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

from astropy import units as u

from ctapipe.visualization import ArrayDisplay
from ctapipe.instrument import SubarrayDescription, TelescopeDescription

# read layout file and plot the array

def plot_telescope_layout(tel_layout_file_name):
    layout_file = "/Users/alicedonini/PycharmProjects/Pointing_tool/layout-3AL4-BN15-MST.txt"
    grd_coords = pd.read_csv(layout_file, header=None, delimiter="  ", names=["x", "y", "z"], engine='python')

    # convert to m
    grd_coords = grd_coords / 100

    # create a dictionary with the tel position on the ground by tel_id and the tel description
    tel_descr = {}
    G_coords = {}
    for tel_id, coord in enumerate(grd_coords.values, 1):
        G_coords[tel_id] = coord * u.m
        tel_descr[tel_id] = TelescopeDescription.from_name(
            optics_name='MST',
            camera_name='NectarCam'
        )

    # create the subarray
    sub = SubarrayDescription(name="Baseline only MST",
                              tel_positions=G_coords,
                              tel_descriptions=tel_descr
                              )
    #sub.info()
    #sub.to_table(kind='optics')

    # display the array
    plt.figure(num=None, figsize=(7, 7), facecolor='w', edgecolor='k')
    disp = ArrayDisplay(sub, tel_scale=3.0)

    return sub


# Pointing generator

class Circle:
    """Class representing a circle"""

    def __init__(self, cx, cy, r):
        """Initialize the circle of radius r and center (cx, cy)"""
        self.cx, self.cy, self.r = cx, cy, r

    def overlap(self, cx, cy, r):
        """Does the circle overlap with another of radius r at (cx, cy)? Yes"""
        d = np.hypot(cx - self.cx, cy - self.cy)
        return d < r + self.r

def plot_circles(circles, r):
    """Plot the circle in a scatter plot"""
    fig, ax = plt.subplots(figsize=(20, 14))

    cx = [circle.cx for circle in circles]
    cy = [circle.cy for circle in circles]
    #r = [circle.r for circle in circles]

    print(cx)
    print(cy)
    for x, y in zip(cx, cy):
        circle1 = plt.Circle((x, y), r, color='r', alpha=0.15, clip_on=False)
        ax.add_patch(circle1)
        plt.scatter(x, y, c="b")
        plt.axis('scaled')
        plt.grid()
        ax.set_aspect('equal')

def place_circle(N, C_X, C_Y, R, r):
    """ Fit the smaller FoVs in a hyper FoV of radius R and center (C_X, C_Y)"""

    circles = []

    # Pick a random position, uniformly on the larger circle's interior
    while len(circles) < N:
        theta = 2 * np.pi * np.random.random()  # random angle between [0,2pi]
        cr = R * np.sqrt(np.random.random())  # random distance < R

        print("theta: ", theta)
        print("cr: ", cr)

        # eq. parametric circumf., random point on the circumference
        cx = cr * np.cos(theta)
        cy = cr * np.sin(theta)

        print("cx: ", cx)
        print("cy: ", cy)
        # The circle should fits inside the HFoV
        if cr + r <= R:
            if len(circles) == 0:
                #place the first one at the center of the HFoV for now
                #circle = Circle(cx + C_X, cy + C_Y, r)
                c = Circle(C_X, C_Y, r)
                circles.append(c)
            else:
                # The circle should overlap any other circle: if yes place it.
                if any(circle.overlap(C_X + cx, C_Y + cy, r) for circle in circles):
                    c = Circle(cx + C_X, cy + C_Y, r)
                    circles.append(c)

    print(circles)
    return circles

def main(tel_layout_file_name):

    FOV = 20
    R = FOV/2
    C_X = 180
    C_Y = 20
    fov = 8
    r = fov / 2
    N = 15
    plot_telescope_layout(tel_layout_file_name)
    circles = place_circle(N, C_X, C_Y, R, r)
    plot_circles(circles, r)
    plt.show()

if __name__ == '__main__':
    filename = "layout-3AL4-BN15-MST.txt"
    main(tel_layout_file_name=filename)
