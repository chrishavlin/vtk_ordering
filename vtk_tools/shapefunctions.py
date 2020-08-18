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

        # dim0,dim1,dim2 = self.element_order                
        # for i in range(dim0+1):
        #     for j in range(dim1+1):
        #         for k in range(dim2+1):
        #             pts.append(self.cell.PointIndexFromIJK(i,j,k))
                        
        for ijk in self._get_ijk_permuatations():
            if hasattr(self.cell,'PointIndexFromIJK'):
                pts.append(self.cell.PointIndexFromIJK(*ijk))
            else:
                pts.append(self._get_cell_id_from_ijk(*ijk))
        
        node_hash = dict(zip(pts,node_nums))
        if self.vtk_create_major_version < 9 and self.cell_type == 72:
            # see https://gitlab.kitware.com/vtk/vtk/-/commit/7a0b92864c96680b1f42ee84920df556fc6ebaa3
            ids2swap = [ [18,19], [30,28], [29,31]]
            for ids in  ids2swap:
                if ids[0] in pts and ids[1] in pts: 
                    node_hash[ids[0]], node_hash[ids[1]] = node_hash[ids[1]], node_hash[ids[0]]

        return node_hash
                   
    def _build_shape_funcs(self):
        # builds sympy expressions for shape function evalulation 
        
        shape_funcs = []
        x=sy.symbols('x')
        y=sy.symbols('y')
        z=sy.symbols('z')                
        pos_i,pos_j,pos_k = self._get_ijk_positions(*self.element_order)
                
        for ijk in self._get_ijk_permuatations():
            LPi = LagrangPoly(x,self.element_order[0],ijk[0],pos_i)
            LPj = 1
            LPk = 1 
            
            if len(ijk) > 1: 
                LPj= LagrangPoly(y,self.element_order[1],ijk[1],pos_j)                
            if len(ijk) > 2: 
                LPk = LagrangPoly(z,self.element_order[2],ijk[2],pos_k)
                
            shape_funcs.append(sy.simplify(LPi * LPj * LPk))
            
        # dim0,dim1,dim2 = self.element_order            
        # for z_i in range(dim0+1):
        #     for y_i in range(dim1+1):
        #         for x_i in range(dim2+1):
        #             LP1 = LagrangPoly(x,dim0,x_i,[-1,0,1])
        #             LP2 = LagrangPoly(y,dim1,y_i,[-1,0,1])
        #             LP3 = LagrangPoly(z,dim2,z_i,[-1,0,1])
        #             shape_funcs.append(sy.simplify(LP1 * LP2 * LP3))
        return shape_funcs 
        
    def _get_ijk_permuatations(self):
        # returns a list of ijk values for looping over element nodes in 1d, 2d, 3d.     
        
        if self.n_dims == 1: 
            ivals = range(0,self.element_order[0]+1)
            return [[ival] for ival in ivals]
        elif self.n_dims > 1: 
            ivals = range(self.element_order[0]+1)
            jvals = range(self.element_order[1]+1)
            
            if self.n_dims == 2: 
                ig,jg = np.meshgrid(ivals,jvals,indexing='ij')
                ig = ig.ravel(order='C')
                jg = jg.ravel(order='C')
                return np.column_stack((ig,jg)).tolist()
            else:
                kvals = range(self.element_order[2]+1)
                ig,jg,kg = np.meshgrid(ivals,jvals,kvals,indexing='ij')
                kg = kg.ravel(order='C')            
                ig = ig.ravel(order='C')
                jg = jg.ravel(order='C')
                return np.column_stack((ig,jg,kg)).tolist()
    
    def _get_ijk_positions(self,order_i,order_j=None,order_k=None):
        # builds list of positions for each coordinate in parent element 
        pos_i = np.linspace(-1,1,order_i+1).tolist()
        pos_j = [] 
        pos_k = [] 

        if order_j is not None: 
            pos_j = np.linspace(-1,1,order_j+1).tolist()
        
        if order_k is not None:
            pos_k = np.linspace(-1,1,order_k+1).tolist()
            
        return pos_i,pos_j,pos_k
                
    
    def _get_cell_id_from_ijk(self,i,j=None,k=None):        
        raise NotImplementedError("Cell type does not hav an ijk mapping yet.")

    def format_shape_functions(self,outputfile = 'shapefunctions.txt',sum=True, fmt = None, progress = True):                        
                        
        def _applyformatting(shape_func,Lvalue):
            # substitution rules, convert to string 
            x, y, z = sy.symbols('x,y,z')
            vl = sy.symbols('v'+str(Lvalue))
            if fmt == 'yt':
                # x,y,z are coord[0],coord[1],coord[2] respectively
                for c in [[x,0],[y,1],[z,2]]:
                    shape_func = shape_func.replace(c[0],sy.Symbol(f'coord[{c[1]}]'))  
                    vl = sy.symbols('values['+str(Lvalue)+']')
                    
                shape_func = vl * shape_func
                shape_func = str(shape_func)
                spaces = ' '.join([' ']*4)
                if sum:
                    sum_str = ' +'
                else:
                    sum_str = ''
                shape_func = spaces + shape_func + sum_str + '\n'
                
            else:
                shape_func = str(shape_func)
                
            return shape_func        
        
        with open(outputfile,'w') as fhandle:
                
            for Lnum,sf in enumerate(self.shape_functions):
                
                sf = _applyformatting(sf,Lnum)
                if progress:
                    print(sf)                 
                fhandle.write(sf)
            
        pass 
