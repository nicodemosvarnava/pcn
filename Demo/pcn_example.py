#!/usr/bin/env python

# Example calculation of the partial Chern number (pcn) for a spinor slab model i.e extended in x-y plane
# and finite in the z direction. After the C(l) list is calculated
# the routine plots the pcn density and integrated pcn in a single output PDF.

from __future__ import print_function # python3 style print
from pythtb import * # import TB model class
from numpy import linalg as LA
import matplotlib.pyplot as plt
from pcn import *

# Example model for axion insulator from Phys. Rev. B 83, 245132
def cubic_model(n_layers):
   # set model parameters
   # lattice parameter implicitly set to a=1
   m   = 1.00
   a   = 0.50
   b   = 0.50
   Bx  = 0.40

   # define lattice vectors
   lat=[[ 1.0, 0.0, 0.0],[ 0.0, 1.0, 0.0],[ 0.0, 0.0, 1.0]]
   # define coordinates of orbitals
   # orb = [site a or 0, site b or 1]
   orb=[[0.0,0.0,0.0], [0.0,0.0,0.0]]
   # make 3d model
   my_model=tb_model(3,3,lat,orb,nspin=2)

   # vector with the pauli matrices
   Id = [[1, 0],[ 0, 1]]
   sx = [[0, 1],[ 1, 0]]
   sy = [[0, -1j],[1j, 0]]
   sz = [[1, 0],[0, -1]]
   sd = [[1, 0],[0, 0]]
   s_pauli = np.zeros((4,2,2),dtype=complex)
   s_pauli  = [Id,sx,sy,sz]

   # set on-site energy
   onsite1 =  (3-m)*np.array(s_pauli[0])+Bx*np.array(sz)
   onsite2 = -(3-m)*np.array(s_pauli[0])-Bx*np.array(sz)
   my_model.set_onsite([onsite1,onsite2])

   # set first neighbour hoppings
   my_model.set_hop(b*1.j*np.array(s_pauli[3]), 1, 0, [-1, 0, 0])
   my_model.set_hop(a*np.array(s_pauli[0]), 1, 0,    [ 0, -1, 0])
   my_model.set_hop(a*1.j*np.array(s_pauli[2]), 1, 0, [ 0, 0, -1])

   my_model.set_hop(-b*1.j*np.array(s_pauli[3]), 1, 0, [ 1, 0, 0])
   my_model.set_hop(-a*np.array(s_pauli[0])    , 1, 0, [ 0, 1, 0])
   my_model.set_hop(-a*1.j*np.array(s_pauli[2]), 1, 0, [ 0, 0, 1])

   #second neighbours
   my_model.set_hop(-b*np.array(s_pauli[0]), 0, 0, [ 1, 0, 0])
   my_model.set_hop(-a*np.array(s_pauli[0]), 0, 0, [ 0, 1, 0])
   my_model.set_hop(-a*np.array(s_pauli[0]), 0, 0, [ 0, 0, 1])

   my_model.set_hop(b*np.array(s_pauli[0]), 1, 1, [ 1, 0, 0])
   my_model.set_hop(a*np.array(s_pauli[0]), 1, 1, [ 0, 1, 0])
   my_model.set_hop(a*np.array(s_pauli[0]), 1, 1, [ 0, 0, 1])

   # create finite model in the z direction
   finite_model=my_model.cut_piece(n_layers,2)
   return finite_model

#Set the model and calculate PCN
################################################
#number of layers of the slab
n_layers= 10
#set the slab model
my_model = cubic_model(n_layers)
print("The model has been constructed..")
#focus on half the slab
n_layer    = int(n_layers/2) #
orbs       = 4 #number of orbital*spin per layer
n_orb001   = orbs*n_layer
n_occ      = int(n_layers*orbs/2) #assuming half filling
#calculate the partial Chern number for each orbital
print("Calling subroutine pcn.. (this may take several minutes)")
cl001 = c_l(my_model,n_occ)[:n_orb001]
#################################################

print("Plotting pcn density and integrated pcn..")
#Define rectangular functions
rect001 = np.zeros((n_layer,n_orb001),dtype=float)
for i in range(n_layer):
    for j in range(i*orbs):
        rect001[i,j] = 0.0
    for j in range(orbs):
        rect001[i,i*orbs+j] = 1.0

#Define ramp functions
ramp001 = np.zeros((n_layer,n_orb001),dtype=float)
for i in range(n_layer):
    for j in range(i*orbs):
        ramp001[i,j] = 1.0
    for j in range(orbs):
        ramp001[i,i*orbs+j] = 0.5

#pcn density
cdens    = np.zeros((n_layer),dtype=float)
for i in range(n_layer):
    cdens[i] = -np.dot(rect001[i,:],cl001[:])

#integrated pcn(partial Chern number)
cint = np.zeros((n_layer),dtype=float)
for i in range(n_layer):
    cint[i] = -np.dot(ramp001[i,:],cl001[:])

#Begin Plotting
#Define x-axis i.e. layer depth
layers = np.linspace(0.5, n_layer-0.5, num=n_layer)
labs=['(a)','(b)']
fig, ax = plt.subplots(1,2,figsize=(8.,3.))
#Plot Chern number per layer
ax[0].plot(layers,cdens,color='blue', marker='o',label = "[001]")
ax[0].axhline(0,c='k',ls='dashed')
# put title
ax[0].set_xlabel("layer")
ax[0].set_ylabel(r'$C(l) $')
ax[0].text(-0.7,0.48,labs[0],size=18.)


#Plot integrated Chern number
ax[1].plot(layers,cint,color='blue', marker='o',label = "[001]")
ax[1].set_ylim(0.25,0.55)
ax[1].axhline(0.5,c='k',ls='dashed')
# put title
ax[1].set_xlabel("layer")
ax[1].set_ylabel(r'$C_{int}(n) $')
ax[1].text(-1.0,0.54,labs[1],size=18.)

# make an PDF figure of a plot
fig.tight_layout()
plt.subplots_adjust(wspace=0.35)
fig.savefig("pcn.pdf")


print('Done.\n')
##############################
