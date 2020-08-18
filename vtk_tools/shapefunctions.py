import sympy as sy 
from .lagrange import LagrangPoly
from .vtk_tools import init_vtk_cell,require_vtk_min_version
import numpy as np 

require_vtk_min_version()

class sf(object):
    def __init__(self,vtk_type,element_order=1,vtk_version='9.0.1'):
        self.cell, self.cell_type = init_vtk_cell(vtk_type)        
        self.n_dims = self.cell.GetCellDimension()
        
        self.vtk_create_version=vtk_version
        self.vtk_create_major_version=int(vtk_version.split('.')[0])
        self.vtk_create_minor_version=int(vtk_version.split('.')[1])
        
        # make element order a list of length n_dims if it's not a list 
        if isinstance(element_order,int):
            element_order = [element_order] * self.n_dims
        self.element_order = element_order 
        
        if hasattr(self.cell,'SetOrder'):
            self.cell.SetOrder(*self.element_order)
        
        for nstr in ['points','edges','faces']:
            attr_name = 'n_'+nstr 
            meth = 'GetNumberOf'+nstr.capitalize()
            setattr(self,attr_name,getattr(self.cell,meth)())

        self.node_order_hash = self._build_point_hash()
        self.shape_functions = self._build_shape_funcs()
        
    def _build_point_hash(self):
        # returns a dict with keys-value pairs of vtk_node_number : ijk node number 
        pts = []
        els = np.array(self.element_order) + 1 
        node_nums = range(0,np.prod(els))

        dim0,dim1,dim2 = self.element_order
        
        for i in range(dim0+1):
            for j in range(dim1+1):
                for k in range(dim2+1):
                    pts.append(self.cell.PointIndexFromIJK(i,j,k))
            
        node_hash = dict(zip(pts,node_nums))
        if self.vtk_create_major_version < 9 and self.cell_type == 72:
            # see https://gitlab.kitware.com/vtk/vtk/-/commit/7a0b92864c96680b1f42ee84920df556fc6ebaa3
            ids2swap = [ [18,19], [30,28], [29,31]]
            for ids in  ids2swap:
                if ids[0] in pts and ids[1] in pts: 
                    node_hash[ids[0]], node_hash[ids[1]] = node_hash[ids[1]], node_hash[ids[0]]

        return node_hash
                   
    def _build_shape_funcs(self):
        
        dim0,dim1,dim2 = self.element_order
        
        shape_funcs = []
        x=sy.symbols('x')
        y=sy.symbols('y')
        z=sy.symbols('z')
        for z_i in range(dim0+1):
            for y_i in range(dim1+1):
                for x_i in range(dim2+1):
                    LP1 = LagrangPoly(x,dim0,x_i,[-1,0,1])
                    LP2 = LagrangPoly(y,dim1,y_i,[-1,0,1])
                    LP3 = LagrangPoly(z,dim2,z_i,[-1,0,1])
                    shape_funcs.append(sy.simplify(LP1 * LP2 * LP3))
        return shape_funcs 
        
