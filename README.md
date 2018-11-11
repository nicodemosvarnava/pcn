# pcn
This is an extention package of PythTB, a package used to simulate crystal structures using the tight-binding approximation http://physics.rutgers.edu/pythtb/.

pcn.py
Routine c_l(my_model,n_occ) calculates the partial Chern number C(l) for each orbital l of a spinor slab model (extended in two direction and finite in the third) provided the number of occupied bands.
The routine assumes dim_k=2 while dim_r=2(one layer) or 3(many layers) and C(l) corresponds to the component of the 
partial Chern vector along the k_1 x k_2 direction.

For the implementation of pcn see:
[1] Phys. Rev. B 98, 115108 - "Geometric and nongeometric contributions to the surface anomalous Hall conductivity"
[2] arXiv:1809.02853 - "Surfaces of axion insulators"


pcn_example.py
Example usage of pcn.py for the calculation of the partial Chern number of a spinor slab model of an axion isnulator. After the C(l) list is calculated, the routine plots the pcn density and integrated pcn in a single output PDF.

For the axion insulator model see:
[3] Phys. Rev. B 83, 245132 - "Inversion-symmetric topological insulators" 


