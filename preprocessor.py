
# Run preprocessor

# docker run -t -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"

# docker run -it -v  /home/arturcastiel/projetos:/pytest desenvolvimento:latest bash -c "cd /pytest; bash"
# %load_ext autoreload
# %autoreload 2
import time
from preprocessor.meshHandle.multiscaleMesh import FineScaleMeshMS as msh
from preprocessor.meshHandle.dualCoarseMesh import DualCoarseMesh as dual
import preprocessor.geoUtil.geoTools as gtool
#import sys
#import imp

# adaptar para 2d
# faces a direita e a esquerda corase
# numerar internos

#print(sys.path)
#sys.path.append('/mesh')
#foobar = imp.load_source('20.h5m', '/mesh')

start = time.time()
M = msh('mesh/semi2.msh', dim = 3)
end = time.time()
print("The preprocessing step lasted {0}s".format(end-start))
#l = dual(M)
