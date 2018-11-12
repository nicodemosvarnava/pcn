# pcn
This is an extention package of PythTB, a package used to simulate crystal structures using the tight-binding approximation http://physics.rutgers.edu/pythtb/.

pcn.py
Routine c_l(my_model,n_occ) calculates the partial Chern number C(l) for each orbital l of a spinor slab model (extended in two direction and finite in the third) provided the number of occupied bands.
The routine assumes dim_k=2 while dim_r=2(one layer) or 3(many layers) and C(l) corresponds to the component of the 
partial Chern vector along the k_1 x k_2 direction.

For the implementation of pcn see:
[1] T. Rauch, T. Olsen, D. Vanderbilt, and I. Souza, "Geometric and nongeometric contributions to the surface anomalous Hall conductivity," Phys. Rev. B 98, 115108 (2018).
[2] N. Varnava and D. Vanderbilt, "Surfaces of axion insulators," arXiv:1809.02853.


pcn_example.py
Example usage of pcn.py for the calculation of the partial Chern number of a spinor slab model of an axion isnulator. After the C(l) list is calculated, the routine plots the pcn density and integrated pcn in a single output PDF.

For the axion insulator model see:
[3] T.L. Hughes, E. Prodan, and B.A. Bernevig, "Inversion-symmetric topological insulators," Phys. Rev. B 83, 245132 (2011).


