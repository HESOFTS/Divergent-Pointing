#!/usr/bin/python

'''
Author: F. Di Pierro (federico.dipierro@to.infn.it)

Goal: for given telescope layout, fov-offset, direction of the center of the hyper-fov, the program produces a multiplicity map in the sky and provide individual telescope theta and phi.

version 0:
The geometry of the fov centers in the sky reproduces the position of the telescope at ground (pointing_type = a1, the only implemented one)

'''

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from pylab import *

#ottengo le opzioni: 
#--telescopes-layout: file.txt simil-corsika
#--pointing-type a1, a2, b, c, d,...
#--fov-offset
#--super-fov-theta_0
#--super-fov-phi_0

__description__ = 'Producing the azimuth and zenith angle for sim_telarray'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-tel_lay', '--telescopes_layout', type=str, required=False, default='layout.txt', help='the input .txt file (corsika.inputs style)')
PARSER.add_argument('-poi_typ', '--pointing_type', type=str, required=False, default='a1', help='Select the pointing option between a1, a2, b, ...')
PARSER.add_argument('-fov_off', '--fov_offset', type=float, required=True, default=1.0, help='offset between adjacent fovs')
PARSER.add_argument('-theta_0', '--theta_0', type=float, required=True, default=20.0, help='zenith of the general pointing center')
PARSER.add_argument('-phi_0', '--phi_0', type=float, required=True, default=180.0, help='azimuth of the general pointing center. As defined in sim_telarray: from North towards East(in corsika PHIP is from South towards East). Examples: N = 0, E = 90, S = 180, W = 270.')
PARSER.add_argument('-tel_fov', '--tel_fov', type=float, required=False, default=8.0, help='Telescopes field-of-view. Diameter, degree.')

def plot_telescope_layout(tel_layout_file_name):
	
	tel_coords = pd.read_csv(tel_layout_file_name, header=None, delimiter = "  ", names = ["x","y","z"])
	tel_coords = tel_coords/100. # cm to m !!!
	x = np.array(tel_coords.x)
	y = np.array(tel_coords.y)
		
	plt.figure(num=None, figsize=(7, 7), facecolor='w', edgecolor='k')	
	plt.scatter(-y, x, color='black', marker='s', s=30)  
	plt.xlabel('-y corsika frame (W->E) [m]')
	plt.ylabel('x corsika frame (S->N) [m]')
	plt.xlim(-600,600)
	plt.ylim(-600,600)
	for i, val in enumerate(x):
		plt.annotate(i, (-y[i],x[i]))
	
	plt.figure(num=None, figsize=(7, 7), facecolor='w', edgecolor='k')	
	plt.scatter(x, y, color='black', marker='s', s=30)  
	plt.xlabel('x corsika frame (S->N) [m]')
	plt.ylabel('y corsika frame (E->W) [m]')
	plt.xlim(-600,600)
	plt.ylim(-600,600)
	for i, val in enumerate(x):
		plt.annotate(i, (x[i],y[i]))

def plot_pointing_map(thetas,phis,plot_type):

	if plot_type =='polar':
		plt.figure(facecolor='white')	
		ax = subplot(111, polar=True)
		plt.scatter(np.radians(phis), thetas, color='red') # che pirla sto matplotlib, plotta in deg ma vuole input in rad
		ax.set_ylim(0.,90.)
		label_position=ax.get_rlabel_position()
		ax.text(np.radians(label_position+10),ax.get_rmax()/2.,'Zenith [deg]',
        rotation=label_position,ha='center',va='center')
		ax.set_theta_direction(-1) # rotate clockwise
		ax.set_theta_offset(pi/2.0)# 0, up (north)
	
	if plot_type =='cartesian':				
		plt.figure(num=None, figsize=(7, 7), facecolor='w', edgecolor='k')	
		#plt.plot(phis, thetas, 'rs')
		plt.scatter(phis, thetas, color='red')
		plt.gca().invert_yaxis()
		plt.xlabel('Azimuth [deg]') 
		plt.ylabel('Zenith [deg]')		
		for i, val in enumerate(phis):# enumerate restituisce la coppia indice, valore i-esimo
		#for i in range(0,len(phis)): # funziona anche cosi`
			plt.annotate(i, (phis[i],thetas[i]))

def my_circle_scatter(axes, x_array, y_array, radius=0.5, **kwargs):
    for x, y in zip(x_array, y_array):
        circle = plt.Circle((x,y), radius=radius, **kwargs)
        axes.add_patch(circle)
    return True

def plot_pointing_map_circles(thetas,phis,plot_type,fov):
	if plot_type =='cartesian':
		fig, ax = plt.subplots(figsize=(6.,6.), facecolor="white" )
		my_circle_scatter(ax,phis,thetas,radius=fov/2,alpha=0.2,color='r')
		plt.gca().invert_yaxis()
		plt.scatter(phis,thetas,c='b',s=4,marker='o',zorder=30)
		plt.grid(True)
		plt.axis('scaled')
		plt.xlabel('Azimuth [deg]') 
		plt.ylabel('Zenith [deg]')		
		plt.title('map of fovs and their centers')
		#plt.xlim((166, 183))
		#plt.ylim((12, 28))
		for i, val in enumerate(phis):# enumerate restituisce la coppia indice, valore i-esimo
			plt.annotate(i, (phis[i],thetas[i]))			

def scale_fovoffset(thetas,phis,offset,theta0,phi0):
	newthetas = thetas*offset + theta0
	newphis = phis*offset + phi0
	return newthetas, newphis 


#A1
def pointing_map_a1(tel_layout_file_name,phi0):
	
	tel_coords = pd.read_csv(tel_layout_file_name, header=None, delimiter = "  ", names = ["x","y","z"])
	tel_coords = tel_coords/100. # cm to m

	new_tel_coords_x = np.array((tel_coords.x - np.mean(tel_coords.x) )/ 175.)
	new_tel_coords_y = np.array((tel_coords.y - np.mean(tel_coords.y) )/ 175.)
	
	#applying roto-translation, from (x,y corsika frame) to (Az,Zd simtel frame)
	thetas_initial = new_tel_coords_x * np.cos(np.radians(phi0)) - new_tel_coords_y * np.sin(np.radians(phi0))
	phis_initial = - new_tel_coords_x * np.sin(np.radians(phi0)) - new_tel_coords_y * np.cos(np.radians(phi0))
			
	return thetas_initial, phis_initial


#A2
def pointing_map_a2():
	
	thetas_initial = np.array()
	phis_initial = np.array()
	
	return thetas_initial, phis_initial
	
#A3

#B

#leggo le posizioni dei telescopi dal file simil-input di corsika
def main(**kwargs):	
	
	if kwargs['pointing_type'] == 'a1' :
		thetas_init,phis_init = pointing_map_a1(kwargs['telescopes_layout'],kwargs['phi_0'])
	
	if kwargs['pointing_type'] == 'a2' :
		thetas_init,phis_init = pointing_map_a2()
	
	thetas,phis = scale_fovoffset(thetas_init,phis_init,kwargs['fov_offset'],kwargs['theta_0'],kwargs['phi_0'])	
		
	for i in range(0, len(thetas)):
		print('#       elif TELESCOPE == ',i+5) # 5 because the first 4 telescopes are LSTs and telescope number starts from 1 instead of 0 
		print('#         include <CTA-ULTRA6-MST-NectarCam.cfg>')
		print('         TELESCOPE_THETA='"%0.2f"%thetas[i])
		print('         TELESCOPE_PHI='"%0.2f"%phis[i])
		print(' ')

	plot_pointing_map(thetas,phis,'polar')
	
	plot_pointing_map(thetas,phis,'cartesian')
	
	plot_pointing_map_circles(thetas,phis,'cartesian',kwargs['tel_fov'])

	plot_telescope_layout(kwargs['telescopes_layout'])	
	
	plt.show()

if __name__ == '__main__':
	args = PARSER.parse_args()
	main(**args.__dict__)
