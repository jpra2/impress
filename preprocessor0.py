
# Run preprocessor

# docker run -t -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"

# docker run -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"
# %load_ext autoreload
# %autoreload 2
import time
from .preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
# from preprocessor.meshHandle.dualCoarseMesh import DualCoarseMesh as dual
from .preprocessor.geoUtil import geoTools as gtool
#import sys
#import imp

#print(sys.path)
#sys.path.append('/mesh')
#foobar = imp.load_source('20.h5m', '/mesh')


name_mesh = 'mesh/27x27x27.msh'

start = time.time()
M = msh(name_mesh, dim = 3)
# M = msh('mesh/malha03.msh', dim = 2)


end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
#l = dual(M)
