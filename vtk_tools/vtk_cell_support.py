import vtk 

# a copy of vtkCellType.h, e.g., https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
vtk_ids = """   
   // Linear cells
   VTK_EMPTY_CELL = 0,
   VTK_VERTEX = 1,
   VTK_POLY_VERTEX = 2,
   VTK_LINE = 3,
   VTK_POLY_LINE = 4,
   VTK_TRIANGLE = 5,
   VTK_TRIANGLE_STRIP = 6,
   VTK_POLYGON = 7,
   VTK_PIXEL = 8,
   VTK_QUAD = 9,
   VTK_TETRA = 10,
   VTK_VOXEL = 11,
   VTK_HEXAHEDRON = 12,
   VTK_WEDGE = 13,
   VTK_PYRAMID = 14,
   VTK_PENTAGONAL_PRISM = 15,
   VTK_HEXAGONAL_PRISM = 16,
  
   // Quadratic, isoparametric cells
   VTK_QUADRATIC_EDGE = 21,
   VTK_QUADRATIC_TRIANGLE = 22,
   VTK_QUADRATIC_QUAD = 23,
   VTK_QUADRATIC_POLYGON = 36,
   VTK_QUADRATIC_TETRA = 24,
   VTK_QUADRATIC_HEXAHEDRON = 25,
   VTK_QUADRATIC_WEDGE = 26,
   VTK_QUADRATIC_PYRAMID = 27,
   VTK_BIQUADRATIC_QUAD = 28,
   VTK_TRIQUADRATIC_HEXAHEDRON = 29,
   VTK_QUADRATIC_LINEAR_QUAD = 30,
   VTK_QUADRATIC_LINEAR_WEDGE = 31,
   VTK_BIQUADRATIC_QUADRATIC_WEDGE = 32,
   VTK_BIQUADRATIC_QUADRATIC_HEXAHEDRON = 33,
   VTK_BIQUADRATIC_TRIANGLE = 34,
  
   // Cubic, isoparametric cell
   VTK_CUBIC_LINE = 35,
  
   // Special class of cells formed by convex group of points
   VTK_CONVEX_POINT_SET = 41,
  
   // Polyhedron cell (consisting of polygonal faces)
   VTK_POLYHEDRON = 42,
  
   // Higher order cells in parametric form
   VTK_PARAMETRIC_CURVE = 51,
   VTK_PARAMETRIC_SURFACE = 52,
   VTK_PARAMETRIC_TRI_SURFACE = 53,
   VTK_PARAMETRIC_QUAD_SURFACE = 54,
   VTK_PARAMETRIC_TETRA_REGION = 55,
   VTK_PARAMETRIC_HEX_REGION = 56,
  
   // Higher order cells
   VTK_HIGHER_ORDER_EDGE = 60,
   VTK_HIGHER_ORDER_TRIANGLE = 61,
   VTK_HIGHER_ORDER_QUAD = 62,
   VTK_HIGHER_ORDER_POLYGON = 63,
   VTK_HIGHER_ORDER_TETRAHEDRON = 64,
   VTK_HIGHER_ORDER_WEDGE = 65,
   VTK_HIGHER_ORDER_PYRAMID = 66,
   VTK_HIGHER_ORDER_HEXAHEDRON = 67,
  
   // Arbitrary order Lagrange elements (formulated separated from generic higher order cells)
   VTK_LAGRANGE_CURVE = 68,
   VTK_LAGRANGE_TRIANGLE = 69,
   VTK_LAGRANGE_QUADRILATERAL = 70,
   VTK_LAGRANGE_TETRAHEDRON = 71,
   VTK_LAGRANGE_HEXAHEDRON = 72,
   VTK_LAGRANGE_WEDGE = 73,
   VTK_LAGRANGE_PYRAMID = 74,
  
   // Arbitrary order Bezier elements (formulated separated from generic higher order cells)
   VTK_BEZIER_CURVE = 75,
   VTK_BEZIER_TRIANGLE = 76,
   VTK_BEZIER_QUADRILATERAL = 77,
   VTK_BEZIER_TETRAHEDRON = 78,
   VTK_BEZIER_HEXAHEDRON = 79,
   VTK_BEZIER_WEDGE = 80,
   VTK_BEZIER_PYRAMID = 81,
"""
vtk_ids = vtk_ids.splitlines() 

vtk_id_hash = {}
for ln in vtk_ids: 
    if '=' in ln : 
        vtk_name,vtk_id = ln.rstrip(',').replace(' ','').split('=')
        vtk_id_hash[int(vtk_id)]=vtk_name
        
        
vtk_att_hash={} 
for vtk_attr in dir(vtk):
    vtk_att_hash[vtk_attr.lower()]=vtk_attr


def _build_vtk_cell_hash():
    """builds a dictionary from vtk cell type number to vtk name and class

    Returns
    -------
    dict
        keys are integers corresponding to vtk type numbers. Each dict entry is a dict 
        with 'vtk_type' and 'vtk_class' strings corresponding to the vtk name and the 
        python-vtk api function. 
        
    $ VTKhash = build_vtk_cell_hash()
    $ VTKhash[72]

    """
    unsupported = []
    CellTypeDict={}
    for ID,att in vtk_id_hash.items():
        attclass = att.lower().split('_')
        attclass = ''.join(attclass)
        
        if attclass in vtk_att_hash.keys():
            attclass = vtk_att_hash[attclass]
            CellTypeDict[ID]={'vtk_type':att,'vtk_class':attclass}
            
        elif 'tetrahedron' in attclass: 
            attclass = attclass.replace('tetrahedron','tetra')
            if attclass in vtk_att_hash.keys():
                attclass = vtk_att_hash[attclass]             
                CellTypeDict[ID]={'vtk_type':att,'vtk_class':attclass}
            else:
                unsupported.append(ID)
        else:
            unsupported.append(ID)

    return CellTypeDict, unsupported

vtk_cell_hash, vtk_unsupported_ids = _build_vtk_cell_hash()
