from vtk_tools.vtk_tools import init_vtk_cell
import vtk_tools.shapefunctions as sfs
import vtk_tools.yaml_support as yams

v_cell,_ = init_vtk_cell(28)
print(type(v_cell))

b = sfs.sf(72,[2,2,2])
print(b.node_order_hash)

b = sfs.sf(72,1)
print(b.node_order_hash)

b = sfs.sf(72,[2,1,2])
print(b.node_order_hash)

b = sfs.sf(72,[2,2,2])
c = sfs.sf(72,[2,2,2],vtk_version='8.0.1')
print([b.node_order_hash[18],c.node_order_hash[19]])
print([b.node_order_hash[19],c.node_order_hash[18]])

b = sfs.sf(72,[1,1,1])
b.write_yaml('testyaml.yaml','w','test_hex8','hex8')

b = sfs.sf(72,[2,2,2])
b.write_yaml('testyaml.yaml','a','test_hex27','hex27')

yt_h = yams.yt_helper()
print("yt mesh types in autogenerator: ")
print(yt_h.available_mesh_types)
print("hex8 shapefunctions:")
print(yt_h.mesh_types['Hex8']['shape_functions'])
