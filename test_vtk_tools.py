from vtk_tools.vtk_tools import init_vtk_cell

v_cell,_ = init_vtk_cell(28)
print(type(v_cell))

from vtk_tools.shapefunctions import sf
b = sf(72,[2,2,2])
print(b.node_order_hash)

b = sf(72,1)
print(b.node_order_hash)

b = sf(72,[2,1,2])
print(b.node_order_hash)

b = sf(72,[2,2,2])
c = sf(72,[2,2,2],vtk_version='8.0.1')
print([b.node_order_hash[18],c.node_order_hash[19]])
print([b.node_order_hash[19],c.node_order_hash[18]])


b.format_shape_functions()
b.format_shape_functions(sum=False)
b.format_shape_functions(fmt='yt')
