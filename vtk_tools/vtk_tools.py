import vtk
from .vtk_cell_support import vtk_cell_hash

def listCells():
    for ID,vals in vtk_cell_hash.items():
        print(f"{vals['vtk_type']} : {ID}")
        
def init_vtk_cell(vtk_type):
    """initializes a vtk cell of a given type. 

    Parameters
    ----------
    vtk_type : int or str
        either the VTK ID or the VTK cell type name

    Returns
    -------
    vtk class : 
        a new instance of a vtk class of type vtk_type 

    To instantiate a 27 node lagrange hexahedron (VTK type 72): 
    
    >>> from vtk_tools.vtk_tools import init_vtk_cell
    >>> new_cell = init_vtk_cell('VTK_LAGRANGE_HEXAHEDRON')
    >>> type(new_cell)
    <class 'vtkmodules.vtkCommonDataModel.vtkLagrangeHexahedron'>

    Or 
        
    >> new_cell = init_vtk_cell(72)
    >>> type(new_cell)
    <class 'vtkmodules.vtkCommonDataModel.vtkLagrangeHexahedron'>
    
    
    """
    
    if isinstance(vtk_type,str):
        if hasattr(vtk,vtk_type):
            vtk_type = getattr(vtk,vtk_type)
        else:
            raise ValueError(f"vtk_type {vtk_type} does not exist.")
        
    if isinstance(vtk_type,int) and vtk_type in vtk_cell_hash.keys():
        cell_info = vtk_cell_hash[vtk_type]
        return getattr(vtk,cell_info['vtk_class'])()
    else:
        raise ValueError(f"vtk_type {vtk_type} is not a valid vtk type.")
    
# https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html        
# https://vtk.org/doc/nightly/html/vtkCellType_8h.html
