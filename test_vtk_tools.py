from vtk_tools.vtk_tools import init_vtk_cell

v_cell = init_vtk_cell(28)
print(type(v_cell))

from vtk_tools.shapefunctions import sf
b = sf(72,[2,2,2])
print(b.node_order_hash)

b = sf(72,1)
print(b.node_order_hash)

b = sf(72,[2,1,2])
print(b.node_order_hash)
