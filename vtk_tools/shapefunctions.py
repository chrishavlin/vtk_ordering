import sympy as sy 
from .lagrage import LagrangPoly
from .vtk_tools import init_vtk_cell
import numpy as np 


class sf(object):
    def __init__(self,vtk_type,element_order=1):
        self.cell = init_vtk_cell(vtk_type)
        self.n_dims = self.cell.GetCellDimension()
        
        # make element order a list of length n_dims if it's not a list 
        if isinstance(element_order,int):
            element_order = [element_order] * self.n_dims
        self.element_order = element_order 
        
        if hasattr(self.cell,'SetOrder'):
            self.cell.SetOrder(*element_order)
        
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
        
        dim0=self.element_order[0]
        dim1=self.element_order[1]
        dim2=self.element_order[2]
        
        for i in range(dim0+1):            
            for j in range(dim1+1):
                for k in range(dim2+1):
                    pts.append(self.cell.PointIndexFromIJK(i,j,k))
                    
        return dict(zip(pts,node_nums))
                   
    def _build_shape_funcs(self):
        
        shape_funcs = []
        x=sy.symbols('x')
        y=sy.symbols('y')
        z=sy.symbols('z')
        for z_i in range(3):
            for y_i in range(3):
                for x_i in range(3):
                    LP1 = LagrangPoly(x,2,x_i,[-1,0,1])
                    LP2 = LagrangPoly(y,2,y_i,[-1,0,1])
                    LP3 = LagrangPoly(z,2,z_i,[-1,0,1])
                    shape_funcs.append(sy.simplify(LP1 * LP2 * LP3))
        return shape_funcs 
        
