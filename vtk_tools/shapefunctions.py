import sympy as sy 
from .lagrange import LagrangPoly
from .vtk_tools import init_vtk_cell,require_vtk_min_version
import yaml 
from . import yaml_support 
import numpy as np 
import os 

require_vtk_min_version()

class sf(object):
    def __init__(self,vtk_type,element_order=1,vtk_version='9.0.1'):
        self.cell, self.cell_type = init_vtk_cell(vtk_type)        
        self.num_dims = self.cell.GetCellDimension()
        
        self.vtk_create_version=vtk_version
        self.vtk_create_major_version=int(vtk_version.split('.')[0])
        self.vtk_create_minor_version=int(vtk_version.split('.')[1])
        
        # make element order a list of length num_dims if it's not a list 
        if isinstance(element_order,int):
            element_order = [element_order] * self.num_dims
        self.element_order = element_order 
        
        # calculate number of vertices 
        els = np.array(self.element_order) + 1 
        self.num_vertices = int(np.prod(els))
        
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
        node_nums = range(0,self.num_vertices)
                        
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
            
        return shape_funcs 
        
    def _get_ijk_permuatations(self):
        # returns a list of ijk values for looping over element nodes in 1d, 2d, 3d.     
        
        if self.num_dims == 1: 
            ivals = range(0,self.element_order[0]+1)
            return [[ival] for ival in ivals]
        elif self.num_dims > 1: 
            ivals = range(self.element_order[0]+1)
            jvals = range(self.element_order[1]+1)
            
            if self.num_dims == 2: 
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

    def write_yaml(self,file,write_mode,mesh_name,mesh_type=None,fmt='yt'):
        """writes element info, including shape functions to a yaml file. 

        Parameters
        ----------
        file : str
            filename to write to 
        write_mode : str
            file write mode, 'w' or 'a'
        mesh_name : str
            what to name this element entry.
        mesh_type : str
            mesh type, required if using fmt='yt'.
        fmt : str
            the yaml format to use. Can be 'yt' or None (the default is 'yt').

        Examples
        --------
        
        to write two elements to a new yaml file with 'yt' formatting: 
        
        import vtk_tools.shapefunctions as sfs
        el_1 = sfs.sf(72,[1,1,1])
        el_1.write_yaml('testyaml.yaml','w','test_hex8','hex8')

        el_2 = sfs.sf(72,[2,2,2])
        el_2.write_yaml('testyaml.yaml','a','test_hex27','hex27')

        """
                
        yaml_dict = {}
        if fmt =='yt':  
            if mesh_type is None:
                raise ValueError("mesh_type must be provided for yt formatting.")      
            yaml_dict[mesh_name] = sf_yt_formatter(
                                                   mesh_type,
                                                   self.num_dims,
                                                   self.num_vertices,
                                                   self.shape_functions
                                                   ).get_dict()   
        else:                     
            formatter = sf_formatter()
            formatter.set_shape_functions(self.shape_functions)
            yaml_dict[mesh_name] = formatter.get_dict()
                        
        with open(file,write_mode) as yaml_stream:
            yaml.dump(yaml_dict,yaml_stream)    
    
class sf_formatter(dict):
    # a strict dictionary with a list of _allowed_keys
    def __init__(self,*args,**kwargs):
        super().__init__(*args, **kwargs)
        self._allowed_keys=['shape_functions']        
        
    def set_empty_keys(self):
        for keyname in self._allowed_keys:
            self.__setitem__(keyname, None)
            
    def __setitem__(self,key,value):
        if key in self._allowed_keys:
            super().__setitem__(key, value)
        else:
            raise KeyError((
                        f"keyname '{key}' is not in list of allowed keys. "
                        f"Allowed keys are: {','.join(self._allowed_keys)}"
                        ))
            
    def _format_sympy(self,shape_func):
        return str(shape_func) 
        
    def set_shape_functions(self,shape_func_list):
        """sets the shape_functions key by processing a list of sympy expressions

        Parameters
        ----------
        shape_func_list : list 
            a list of sympy expressions representing the shapefiles

        """
        
        if len(shape_func_list) != self.get('num_vertices') and len(shape_func_list)!=0:
            raise ValueError("length of shape file list must match number of vertices!")
        
        shape_functions=[self._format_sympy(sf)+',' for sf in shape_func_list]        
        self.__setitem__('shape_functions',shape_functions)
        
    def get_dict(self):
        # returns a normal dictionary 
        return {key : value for key,value in self.items()}
        
class sf_yt_formatter(sf_formatter):
    def __init__(self,mesh_type,num_dim,num_vertices,shape_func_list=[],*args,**kwargs):
        super().__init__(*args,**kwargs)        
        self._allowed_keys = ['mesh_type', 'num_dim', 'num_vertices', 
                           'num_mapped_coords', 'shape_functions']
        self.set_empty_keys()
        self.__setitem__('mesh_type',mesh_type)
        self.__setitem__('num_dim',num_dim)
        self.__setitem__('mesh_type',mesh_type)
        self.__setitem__('num_mapped_coords',num_dim)
        self.__setitem__('num_vertices',num_vertices)
        self.set_shape_functions(shape_func_list)
        
                    
    def _format_sympy(self,shape_func):
        # substitution rules, convert to string 
        x, y, z = sy.symbols('x,y,z')
    
        # x,y,z are x[0],x[1],x[2] respectively
        for c in [[x,0],[y,1],[z,2]]:
            shape_func = shape_func.replace(c[0],sy.Symbol(f'x[{c[1]}]'))  
                                
        return str(shape_func)
        
    def set_shape_functions(self,shape_func_list):
        super().set_shape_functions(shape_func_list)
        sf_list = self.get('shape_functions')
        sf_str = '['+'\n'.join(sf_list).rstrip(',') + ']'        
        self.__setitem__('shape_functions',yaml_support.literal_str(sf_str))
        
