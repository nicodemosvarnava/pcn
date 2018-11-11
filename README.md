# pcn
pcn.py
Routine that calculates the partial Chern number C(l) for each orbital l of a spinor slab model (extended in two direction 
and finite in the third) provided the number of occupied bands.
The routine assumes dim_k=2 while dim_r=2(one layer) or 3(many layers) and C(l) corresponds to the component of the 
partial Chern vector along the k_1 x k_2 direction.

pcn_example.py
Example calculation of the partial Chern number (pcn) for a spinor slab model. After the C(l) list is calculated
the routine plots the pcn density and integrated pcn in a single output PDF.
