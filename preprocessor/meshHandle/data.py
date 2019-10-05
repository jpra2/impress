## criado por joao paulo
'''
classe Data criada para armazenar um dict linkando o nome da variavel
a variavel do impress e tambem as informacoes dos dados.
Essa classe será utilizada para automatizar alguns processos.
'''
import numpy as np
from .. import directories as direc
import pickle

class Data:
    '''
    armazena as variaveis do impress e dicts que linkam
    as entidades do moab aos seus respectivos ids globais
    '''

    def __init__(self, n_nodes, n_edges, n_faces, n_volumes, fine_scale_mesh_obj):
        '''
        n_nodes: numero de nos
        n_edges: numero de arestas
        n_faces: numero de faces
        n_volumes: numero de volumes
        fine_scale_mesh_obj: objeto multiscalemeshMS
        '''

        self.variables = dict()
        self.info_data = dict()
        self.len_entities = dict()

        self.len_entities[direc.entities_lv0[0]] = n_nodes
        self.len_entities[direc.entities_lv0[1]] = n_faces
        self.len_entities[direc.entities_lv0[2]] = n_edges
        self.len_entities[direc.entities_lv0[3]] = n_volumes
        self.mesh = fine_scale_mesh_obj
        self.name_variables = direc.path_local_variables
        self.name_info_data = direc.path_local_info_data

    def get_info_data(self, name, data_size, data_format, entity, level=0):
        '''
        name: nome da variavel
        impress_variable: variavel do impress
        data_size: tamanho da variavel
        data_format: formato da variavel
        entity: entidade onde foi setada
        level: nivel

        dict que armazena as informacoes das variaveis do impress
        '''
        # self.variables[name] = impress_variable
        info = dict()

        info[direc.names_datas[0]] = int(data_size)
        info[direc.names_datas[1]] = data_format
        info[direc.names_datas[2]] = entity
        info[direc.names_datas[3]] = level

        self.info_data[name] = info

    def init_datas(self, names=None):
        '''
        zera todas as variaveis do impress e armazena em self.variables
        names deve ser obrigatoriamente uma lista de strings
        '''

        variables = dict()

        for name, infos in self.info_data.items():
            n = infos[direc.names_datas[0]]
            format = infos[direc.names_datas[1]]
            entity = infos[direc.names_datas[2]]
            n_entity = self.len_entities[entity]

            if format == direc.data_formats[0]:
                data = np.zeros(n_entity)
            elif format == direc.data_formats[1]:
                data = np.zeros(n_entity, dtype=np.int32)
            if n > 1:
                data = np.repeat(data, n).reshape([n_entity, n])

            variables[name] = data

        self.variables = variables

    def init_dicts(self):
        '''
        dict_elements = dict com os dicts que linkam as entidades do pymoab aos
        seus ids globais
        mesma coisa para os centroides
        '''


        self.dict_elements = dict()

        try:
            dict_elem = dict(zip(self.mesh.core.all_volumes, self.mesh.volumes.all))
            self.dict_elements[direc.entities_lv0[3]] = dict_elem
        except:
            pass

        dict_elem = dict(zip(self.mesh.core.all_faces, self.mesh.faces.all))
        self.dict_elements[direc.entities_lv0[2]] = dict_elem
        dict_elem = dict(zip(self.mesh.core.all_edges, self.mesh.edges.all))
        self.dict_elements[direc.entities_lv0[1]] = dict_elem
        dict_elem = dict(zip(self.mesh.core.all_nodes, self.mesh.nodes.all))
        self.dict_elements[direc.entities_lv0[0]] = dict_elem

        self.elements_lv0 = dict()

        internal_faces = self.mesh.faces.internal
        self.elements_lv0[direc.entities_lv0_0[0]] = internal_faces
        vols_viz_faces = self.mesh.faces.bridge_adjacencies(self.mesh.faces.all, 2, 3)
        self.elements_lv0[direc.entities_lv0_0[1]] = vols_viz_faces

        vols_viz_internal_faces = np.zeros((len(internal_faces), 2), dtype=np.int32)
        for i, f in enumerate(internal_faces):
            vols_viz_internal_faces[i] = vols_viz_faces[f]
        self.elements_lv0[direc.entities_lv0_0[2]] = vols_viz_internal_faces

        boundary_faces = np.setdiff1d(self.mesh.faces.all, internal_faces)

        vols_viz_boundary_faces = np.zeros((len(boundary_faces)), dtype=np.int32)
        for i, f in enumerate(boundary_faces):
            vols_viz_boundary_faces[i] = vols_viz_faces[f]

        self.elements_lv0[direc.entities_lv0_0[3]] = vols_viz_boundary_faces

        self.centroids = dict()
        '''
        self.centroids = dicts para linkar o nome das entidades aos seus centroides
        '''

        try:
            self.centroids[direc.entities_lv0[3]] = self.mesh.volumes.center(self.mesh.volumes.all)
        except:
            pass

        self.centroids[direc.entities_lv0[2]] = self.mesh.faces.center(self.mesh.faces.all)
        self.centroids[direc.entities_lv0[1]] = self.mesh.edges.center(self.mesh.edges.all)
        self.centroids[direc.entities_lv0[0]] = self.mesh.nodes.center(self.mesh.nodes.all)

    def update_to_mesh(self, names=None):
        if names:
            for name in names:
                command = 'self.mesh.' + name + '[:] = ' + 'self.variables["' + name + '"]'
                exec(command)

        else:
            for name in self.variables.keys():
                command = 'self.mesh.' + name + '[:] = ' + 'self.variables["' + name + '"]'
                exec(command)

    def load_from_mesh(self, names=None):
        if names:
            for name in names:
                command = 'self.variables["' + name + '"] = self.mesh.' + name + '[:]'
                exec(command)

        else:
            for name in self.variables.keys():
                command = 'self.variables["' + name + '"] = self.mesh.' + name + '[:]'
                exec(command)

    def export_to_npz(self):

        name_variables = self.name_variables
        name_info_data = self.name_info_data

        np.savez(name_variables, **self.variables)

    def load_from_npz(self):
        name_variables = self.name_variables
        name_info_data = self.name_info_data

        arq = np.load(name_variables)

        for name, variable in arq.items():
            self.variables[name] = variable
