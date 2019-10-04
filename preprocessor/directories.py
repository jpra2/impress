import os

impress_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

entities_lv0 = ['nodes', 'edges', 'faces', 'volumes']
entities_lv0_0 = ['internal_faces', 'vols_viz_faces']
names_datas = ['data_size', 'data_format', 'entity', 'level']
data_formats = ['float', 'int', 'bool']
name_variables = 'variables.npz'

flying_impress = 'flying'
path_local_variables = flying_impress + '/' + name_variables
