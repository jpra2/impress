"""
Schemes of partitioning
"""

#import pdb
import numpy as np
from numba import jit

#para definir funções de particionamento deve-se atentar as saidas
#é necessário retornar um vetor tag que associa uma partição a cada celula fina
#também é necessario retornar o uma matrix
# n x 3 com o centro de cada volume coarse, n -> numero de elementos do volume coarse

#a função tagAdjust remove renumera os elementos da malha coarse de maneira que eles permanecam
# continuos contando de 0 ao N de volumes coarse, a função tb remover os centros indesejados

#para adicionar novos esquemas basta criar um esquema com numeracao sequencial
#ex. scheme2  , mais shceme3
#a seção campo de leitura do arquivo msCoarse.ini deverá ter nome correspondente
#ex [Coarsening_2_Input] e [Coarsening_3_Input]
#todos atributos serao passados para funcao correspondente na ordem definida
#por default a leitura dos elementos é float, caso necessário. converta para int
#ex int(nx)

def scheme1(centerCoord, num_of_vol, rx,ry,rz ,nx = 3, ny = 3, nz =3, *argv, **kwargs):

    #input : centerCoord - > array with the center of elements
    #        num_of_vol = number of volumes

    #        rx,ry,rz - (min, max) values of x,y,z of the phyisical domain
    #        nx, ny, nz
    # msh -> objeto da clase meshUtil
    #centerCoord = msh.readData("CENTER")

    nx = int(nx)
    ny = int(ny)
    nz = int(nz)
    if (rz[1] == 0)  &  (rz[0] == 0):
        nz = 1
        rz = (-1,1)

    box = np.array([0, (rx[1] - rx[0])/nx, 0,(ry[1] - ry[0]) /ny, 0,(rz[1] - rz[0])/(nz+0)]).reshape(3,2)
    cent_coord_El1 = box.sum(axis =1)/2
    tag = np.zeros(num_of_vol).astype("int")
    coarseCenters = np.zeros((nx*ny*nz,3))
    index = 0
    init_coords = np.array([rx[0],ry[0],rz[0]])
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                inc = np.multiply(box[:,1], np.array([x,y,z]))

                #cent = cent_coord_El1 + inc
                coarseCenters[index] = cent_coord_El1 + inc
                # pdb.set_trace()

                #inc = np.array([(nx) * x, (ny) * y, (nz) * z])
                boxMin = box[:,0] + inc + init_coords
                boxMax = box[:,1] + inc + init_coords
                point = checkinBox(centerCoord,x=(boxMin[0], boxMax[0]), y=(boxMin[1], boxMax[1]) , z=(boxMin[2], boxMax[2]))
                tag[point] = index
                index += 1
    return tagAdjust(tag,coarseCenters)

def scheme2(centerCoord, num_of_vol, rx,ry,rz ,nx = 3, ny = 3, nz =3 ):
    #input : centerCoord - > array with the center of elements
    #        num_of_vol = number of volumes

    #        rx,ry,rz - (min, max) values of x,y,z of the phyisical domain
    #        nx, ny, nz
    # msh -> objeto da clase meshUtil
    #centerCoord = msh.readData("CENTER")
    print("ESQUEMA TIPO 2 PARTICIONAMENTO")
    nx = int(nx)
    ny = int(ny)
    nz = int(nz)
    box = np.array([0, (rx[1] - rx[0])/nx, 0,(ry[1] - ry[0]) /ny, 0,(rz[1] - rz[0])/(nz+0)]).reshape(3,2)
    cent_coord_El1 = box.sum(axis =1)/2
    tag = np.zeros(num_of_vol).astype("int")
    coarseCenters = np.zeros((nx*ny*nz,3))
    index = 0

    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                inc = np.multiply(box[:,1], np.array([x,y,z]))
                #cent = cent_coord_El1 + inc
                coarseCenters[index] = cent_coord_El1 + inc
                #inc = np.array([(nx) * x, (ny) * y, (nz) * z])
                boxMin = box[:,0] + inc
                boxMax = box[:,1] + inc
                point = checkinBox(centerCoord,x=(boxMin[0], boxMax[0]), y=(boxMin[1], boxMax[1]) , z=(boxMin[2], boxMax[2]))
                tag[point] = index
                index += 1
    #return tagAdjust(tag,coarseCenters)

def scheme3(centerCoord, num_of_vol, rx,ry,rz ,nx = 3, ny = 3, nz =3 ):
    #input : centerCoord - > array with the center of elements
    #        num_of_vol = number of volumes

    #        rx,ry,rz - (min, max) values of x,y,z of the phyisical domain
    #        nx, ny, nz
    # msh -> objeto da clase meshUtil
    #centerCoord = msh.readData("CENTER")
    print("ESQUEMA 3")
    nx = int(nx)
    ny = int(ny)
    nz = int(nz)
    box = np.array([0, (rx[1] - rx[0])/nx, 0,(ry[1] - ry[0]) /ny, 0,(rz[1] - rz[0])/(nz+0)]).reshape(3,2)
    cent_coord_El1 = box.sum(axis =1)/2
    tag = np.zeros(num_of_vol).astype("int")
    coarseCenters = np.zeros((nx*ny*nz,3))
    index = 0
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                inc = np.multiply(box[:,1], np.array([x,y,z]))
                #cent = cent_coord_El1 + inc
                coarseCenters[index] = cent_coord_El1 + inc
                #inc = np.array([(nx) * x, (ny) * y, (nz) * z])
                boxMin = box[:,0] + inc
                boxMax = box[:,1] + inc
                point = checkinBox(centerCoord,x=(boxMin[0], boxMax[0]), y=(boxMin[1], boxMax[1]) , z=(boxMin[2], boxMax[2]))
                tag[point] = index
                index += 1
    return tagAdjust(tag,coarseCenters)

############################
## adicionado por jp
def scheme4(centerCoord, num_of_vol, rx,ry,rz ,nx = 3, ny = 3, nz =3, *argv, **kwargs):
    name_coord_nodes = 'coord_nodes'
    coord_nodes = kwargs.get(name_coord_nodes, None)
    if coord_nodes is None:
        raise ValueError(f'\n a chave {name_coord_nodes} deve conter as coordenadas dos nós \n')

    centerCoordTuple = (tuple(x) for x in centerCoord)
    map_center_coord = dict(zip(centerCoordTuple, range(len(centerCoord))))

    boxes_points, cent_boxes = get_boxes(nx, ny, nz, rx, ry, rz)
    # boxes_points, cent_boxes = get_boxes_points(coord_nodes, boxes_points)
    fine_elements_in_boxes = get_fine_elements_in_boxes(boxes_points, centerCoord)
    # fine_elements_in_boxes, cent_boxes = test_fine_elements(fine_elements_in_boxes, centerCoord, boxes_points, cent_boxes)
    fineTag = get_fineTag(fine_elements_in_boxes, map_center_coord, centerCoord, cent_boxes)

    return fineTag, cent_boxes
###########################


#checa se um ponto esta dentro de um cubo

@jit(parallel=True)
def checkinBox(coords, x , y, z):
    tag1 = (coords[:,0] > x[0])   &  (coords[:,0] < x[1])
    tag2 = (coords[:,1] > y[0])   &  (coords[:,1] < y[1])
    tag3 = (coords[:,2] > z[0])   &  (coords[:,2] < z[1])
    return tag1 & tag2 & tag3



#função para corrigir os tags em caso que as malhas geradoras possuam volumes grossos sem celulas dentro
#e remover as respectivas coordenadas do centro dos volumes das malhas primais
#  @jit(parallel=True, cache=True)
def tagAdjust(tag, coarseCenter):
    # msh -> objeto da clase meshUtil
    fineTag =  tag
    elementsOriginal = [*set(tag)]
    elementsNovo = [*set(range(len(elementsOriginal)))]
    elementsMissing = set(range(len(coarseCenter))) - set(elementsOriginal)
    for elo, eln in zip(elementsOriginal,elementsNovo):
        if elo != eln:
            pointer = (tag == elo)
            fineTag[pointer] = eln

    return fineTag.astype(int) , np.delete(coarseCenter, [*elementsMissing], axis = 0)

####################
## adicionado por jp
def get_boxes(nx, ny, nz, rx, ry, rz):

    rs = [rx[1], ry[1], rz[1]]
    ns = [nx, ny, nz]
    rests = []
    divs = []
    points = []

    for i in range(3):
        rests.append(int(rs[i] % ns[i]))
        divs.append(int(rs[i] // ns[i]))

    for i in range(3):
        d = divs[i]
        r = rests[i]
        n = ns[i]
        ri = rs[i]
        if r == 0:
            pp = np.linspace(0.0, ri, n+1)
        else:
            pp = [0.0]
            for j in range(n-1):
                pp.append(pp[j] + d)
            pp.append(ri)
            pp = np.array(pp)

        points.append(pp)

    points = np.array(points)

    boxes = []
    # cent_boxes = []

    # n2 = len(points[2])-1
    # n1 = len(points[1])-1
    # n0 = len(points[0])-1
    # ps0 = points[0]
    # ps1 = points[1]
    # ps2 = points[2]

    # boxes = np.array([np.array([np.array([ps0[i], ps1[j], ps2[k]]), np.array([ps0[i+1], ps1[j+1], ps2[k+1]])]) for k in range(n2) for j in range(n1) for i in range(n0)])
    # boxes3 = [np.array([ps0[i+1], ps1[j+1], ps2[k+1]]) for k in range(n2) for j in range(n1) for i in range(n0)]


    for k in range(len(points[2]) - 1):
        for j in range(len(points[1]) - 1):
            for i in range(len(points[0]) - 1):
                p1 = np.array([points[0][i], points[1][j], points[2][k]])
                p2 = np.array([points[0][i+1], points[1][j+1], points[2][k+1]])
                boxes.append(np.array([p1, p2]))

    boxes = np.array(boxes)
    cent_boxes = np.mean(boxes, axis=1)
    return boxes, cent_boxes

def get_boxes_points(coord_nodes, boxes):

    coord2 = coord_nodes
    ids = np.arange(len(coord2))

    boxes_points = np.zeros(boxes.shape)
    cent_boxes = []

    for i, box in enumerate(boxes):
        p0 = box[0]
        p1 = box[1]

        pp0 = get_prox(coord2, p0, ids)
        pp1 = get_prox(coord2, p1, ids)

        boxes_points[i][0] = pp0
        boxes_points[i][1] = pp1

    cent_boxes = np.mean(boxes_points, axis=1)
    return boxes_points, cent_boxes

def get_prox(coord2, p, ids):
    tt = coord2 - p
    n = np.linalg.norm(tt, axis=1)
    minn = n.min()
    ind = ids[n == minn][0]
    pp = coord2[ind]
    return pp

def get_fine_elements_in_boxes(boxes_points, centerCoord):

    fine_elements_in_boxes = []
    ids = np.arange(len(centerCoord))

    for box in boxes_points:
        indices = get_box(centerCoord, box, ids)
        fine_elements_in_boxes.append(centerCoord[indices])

    return np.array(fine_elements_in_boxes)

def get_box_dep0(all_centroids, limites):

    '''
    all_centroids->coordenadas dos centroides do conjunto
    limites-> diagonal que define os volumes objetivo (numpy array com duas coordenadas)
    Retorna os indices cujo centroide está dentro de limites
    '''

    inds0 = np.where(all_centroids[:,0] > limites[0,0])[0]
    inds1 = np.where(all_centroids[:,1] > limites[0,1])[0]
    inds2 = np.where(all_centroids[:,2] > limites[0,2])[0]
    c1 = set(inds0) & set(inds1) & set(inds2)
    inds0 = np.where(all_centroids[:,0] < limites[1,0])[0]
    inds1 = np.where(all_centroids[:,1] < limites[1,1])[0]
    inds2 = np.where(all_centroids[:,2] < limites[1,2])[0]
    c2 = set(inds0) & set(inds1) & set(inds2)
    inds_vols = list(c1 & c2)
    return inds_vols

def get_box(all_centroids, limites, ids):

    '''
    all_centroids->coordenadas dos centroides do conjunto
    limites-> diagonal que define os volumes objetivo (numpy array com duas coordenadas)
    Retorna os indices cujo centroide está dentro de limites
    '''

    inds_vols = ids[
            (all_centroids[:,0] > limites[0,0]) & (all_centroids[:,1] > limites[0,1]) & \
            (all_centroids[:,2] > limites[0,2]) & (all_centroids[:,0] < limites[1,0]) & \
            (all_centroids[:,1] < limites[1,1]) & (all_centroids[:,2] < limites[1,2])]
    return inds_vols

def test_fine_elements(fine_elements_in_boxes, centerCoord, boxes, center_boxes):

    centers2 = centerCoord

    gg = fine_elements_in_boxes.flatten().reshape(len(centerCoord), 3)
    gg2 = np.setdiff1d(centers2, gg)

    if len(gg2) > 0:
        idsc = np.arange(len(center_boxes))
        for p in gg2:
            diff = center_boxes - p
            minn = np.linlalg.norm(diff, axis=1)
            indice = idsc[diff == minn][0]
            dd = list(fine_elements_in_boxes[indice])
            dd.append(p)
            fine_elements_in_boxes[indice] = np.array(dd)
            center = np.mean(fine_elements_in_boxes[indice], axis=0)
            center_boxes[indice] = center

    return fine_elements_in_boxes, center_boxes

def get_fineTag(fine_elements_in_boxes, map_center_coord, centerCoord, center_boxes):

    fineTag = np.zeros(len(centerCoord), dtype=np.int32)
    for i, elements in enumerate(fine_elements_in_boxes):
        els = (tuple(x) for x in elements)
        indices = np.array([map_center_coord[k] for k in els])
        fineTag[indices] = np.repeat(i, len(indices))

    return fineTag
####################
