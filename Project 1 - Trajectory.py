"""
Mac OS X must use tk version 8.6.7.  See project instructions on how to switch.
"""
from graphics import *
import time
import math

""" Write your code here """
#obtain and assign input values
radius = float(input('Enter a reasonable radius (10-20):\n'))
v0 = float(input('Enter a reasonable initial speed (30-80):\n'))
theta_raw = float(input('Enter a reasonable angle (10-90):\n'))
time_inc = float(input('Enter a time increment:\n'))

#theta input as angle, however, must convert to radians for math methods
theta = theta_raw * math.pi / 180

#defining the constant g
g = 9.8

#performing calculations for necessary constants
max_height = math.pow(v0, 2) * math.pow(math.sin(theta), 2) / (2 * g)
max_range = math.pow(v0, 2) * (math.sin(2 * theta) / g)
window_height = math.floor(max_height + (3 * radius))
window_width = math.floor(max_range + (3 * radius))
t_flight = 2 * v0 * math.sin(theta) / g
num_moves = int(t_flight // time_inc)

#create window and object based on in put values
win = GraphWin('Trajectory', window_width, window_height)
p1 = Point(0, window_height)
C = Circle(p1, radius)
C.setFill("yellow")
C.draw(win)

#initializing variables needed to display trajectory
dx = 0
dy = 0
t = time_inc
iteration = 1

#for loop: sleeps for the time increment, then calculates dx and dy
#each iteration. Moves the object according to dx and dy. Prints out the
#iteration number as well as the dx and dy that the object was moved by.
for x in range(num_moves):
    time.sleep(time_inc)
    dx = v0 * math.cos(theta) * time_inc
    dy = -((v0 * math.sin(theta)) - g * t) * time_inc
    C.move(dx, dy)
    print('Iteration', iteration, 'dx = {:.3f}, dy = {:.3f}'.format(dx, -dy))
    t += time_inc
    iteration += 1

#Close after mouse click
try:
    win.getMouse()    
    win.close()
except:
    pass
    
