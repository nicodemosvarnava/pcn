#!/usr/bin/env python

########################
#  pcn.py Version 1.0  #
########################

# Routine that calculates the partial Chern number C(l) for each orbital l
# of a spinor slab model (extended in two direction and finite in the third)
# provided the number of occupied bands.
#
# The routine assumes dim_k=2 while dim_r=2 (one layer) or 3 (many layers)
# and C(l) corresponds to the component of the partial Chern vector along
# the k_1 x k_2 direction.

from pythtb import * # import TB model class
import numpy as np

def c_l(my_model,n_occ):
    #length of the uniform k-mesh centered at gamma
    kmesh = 20
    n_spin=my_model._nspin
    n_orb=my_model._norb
    dim_k=my_model._dim_k
    n_bands = n_orb*n_spin
    n_unoc = n_bands - n_occ

    #calculate the number of hoppings
    n_hop = 0
    for hop in my_model._hoppings:
        for ind_1 in range(2):
            for ind_2 in range(2):
                me=hop[0][ind_1,ind_2]
                if me != 0:
                    n_hop += 1

    # generate uniform k-mesh
    kpts_uni=[]
    n_k = kmesh*kmesh
    sym_kmesh = int(kmesh/2)
    for i in range(kmesh):
        for j in range(kmesh):
            k_vec = [(-sym_kmesh+float(i))/float(kmesh),(-sym_kmesh+float(j))/float(kmesh)]
            kpts_uni.append(k_vec)

    #solve system on a uniform 2dk-grid
    print("     Calculating energy eigenvalues E_nk and eigenvectors \psi_nk..")
    (evals,evecs)=my_model.solve_all(kpts_uni,eig_vectors=True)

    #flatten evecs
    evec = np.zeros((n_bands,n_k,n_bands),dtype='complex')
    for n in range(n_bands):
        for k in range(n_k):
            evec[n,k] = [item for sublist in evecs[n,k] for item in sublist]

    #these list contain the hopping parameters needed to calculate pos_vc
    mes = np.zeros((n_k,n_hop,2),dtype='complex')
    mcc = np.zeros((n_k,n_hop,2),dtype='complex')
    pos_ij = np.zeros((n_hop,2),dtype='int')
    count = 0
    for hop in my_model._hoppings:
        for ind_1 in range(2):
            for ind_2 in range(2):
                me=hop[0][ind_1,ind_2]
                if me != 0:
                    i=hop[1]*2+ind_1
                    j=hop[2]*2+ind_2
                    # for each hopping save the orbital position i,j
                    pos_ij[count,0] = i
                    pos_ij[count,1] = j
                    # get bond vector
                    r_bra=my_model.get_orb()[hop[1],:2]
                    r_ket=my_model.get_orb()[hop[2],:2]
                    hh = hop[3][:2]
                    bond_vec=hh+(r_ket-r_bra)
                    for k in range(n_k):
                        p=np.dot(bond_vec,kpts_uni[k])
                        mee=me*np.exp(2.j*np.pi*p)
                        s=-bond_vec
                        mes[k,count,:] = mee*s
                        mcc[k,count,:] = np.conj(mee)*(-s)
                    count += 1

    #calculate pos_vc
    print("     Calculating position matrices X_vc(k) and Y_vc(k)..")
    pos_vc = np.zeros((dim_k,n_occ,n_unoc,n_k),dtype='complex')
    for v in range(n_occ):
        for cc in range(n_unoc):
            c = cc + n_occ
            for k in range(n_k):
                en_dif  =  evals[c,k] - evals[v,k]
                for l in range(dim_k):
                    for count in range(n_hop):
                        i = pos_ij[count,0]
                        j = pos_ij[count,1]
                        pos_vc[l,v,cc,k] += np.conjugate(evec[v,k,i])* mes[k,count,l] * evec[c,k,j]/en_dif
                        pos_vc[l,v,cc,k] += np.conjugate(evec[v,k,j])* mcc[k,count,l] * evec[c,k,i]/en_dif

    #covariant curvature tensor F_vv'
    print("     Calculating covariant curvature tensor F_vv(k)..")
    f_vvp = np.zeros((n_occ,n_occ,n_k),dtype='complex')
    for v in range(n_occ):
        for vv in range(n_occ):
            for c in range(n_unoc):
                for k in range(n_k):
                    f_vvp[v,vv,k] += pos_vc[0,v,c,k]*np.conjugate(pos_vc[1,vv,c,k])
    #Calculate C(l)
    print("     Calculating partial Chern number C(l)..")
    c = np.zeros((n_bands),dtype=float)
    for v in range(n_occ):
        for vv in range(n_occ):
            for k in range(n_k):
                for l in range(n_bands):
                    c[l] += -4*np.pi*np.imag(evec[v,k,l]*np.conjugate(evec[vv,k,l])*f_vvp[v,vv,k])/float(n_k)
    return c
