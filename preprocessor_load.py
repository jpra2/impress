
def init_mesh(mesh_name):

    import time
    import os
    from .preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
    # from preprocessor.meshHandle.dualCoarseMesh import DualCoarseMesh as dual
    from .preprocessor.geoUtil import geoTools as gtool
    #import sys
    #import imp

    #print(sys.path)
    #sys.path.append('/mesh')
    #foobar = imp.load_source('20.h5m', '/mesh')

    start = time.time()
    M = msh(mesh_name, dim = 3)
    # M = msh('mesh/malha03.msh', dim = 2)


    end = time.time()
    print("The preprocessing step lasted {0}s".format(end-start))

    return M
