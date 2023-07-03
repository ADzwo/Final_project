import re
import json
from tuples import *

class Genome:
    path: list # Path: vertex index: int, orientation: int, used: bool
    
    def __init__(self, path):
        self.path = path


class Vertex:
    length: int # length of the represented sequence
    occurrences: list # occurrence: genome index, index on its path

    def __init__(self, length, occurrences=None):
        self.v_length = length
        if not occurrences:
            self.occurrences = []
        else:
            self.occurrences = occurrences

    def add_occurrence(self, g_idx, nr_on_path):
        self.occurrences.append(occurrence(g_idx, nr_on_path))


class Graph:
    genomes: list # list of Genome objects
    vertices: list # list of Vertex objects
    name: str

    def __init__(self, graph_file_path):
        self.vertices = []
        self.genomes = []
        vertex_name_to_idx = {} # dict to save {vertex id from gfa file: index in Graph.vertices}
        genome_name_to_idx = {}
        sequences = []
        self.name = re.split(r'[\/]', graph_file_path)[-1]
        self.name = re.split(r'\.', self.name)[0]

        # find S lines and fill vertex_name_to_idx dict 
        with open(graph_file_path,'r') as graph_file:
            for line in graph_file.readlines():
                if line.startswith('S'):
                    v_name, v_sequence = self.add_vertex(line)
                    if v_name is not None:
                        vertex_name_to_idx[v_name] = len(self.vertices)-1
                        sequences.append(v_sequence)
        with open(f'vertex_sequences/{self.name}.txt', 'w') as f:
            f.writelines(sequences)
            del sequences
        # find P lines and fill genome_name_to_idx dict 
        with open(graph_file_path,'r') as graph_file:
            for line in graph_file.readlines():
                if line.startswith('P'): # or 'W'
                    g = line.strip().split() # P, name, vertices' names, overlaps
                    path = [] # list[Path] - vertex, orientation, used
                    p_length = 0
                    for v_pos, vertex in enumerate(g[2].split(',')):
                        v_name = vertex[:-1]
                        if v_name not in vertex_name_to_idx:
                            raise ValueError(f'File {graph_file_path} contains undefined sequence {v_name} in a genome.')
                        v_idx = vertex_name_to_idx[v_name]
                        v_orientation = 1 if vertex[-1]=='+' else -1
                        self.vertices[v_idx].add_occurrence(len(self.genomes), v_pos) # genome, nr_on_path
                        p_length += self.vertices[v_idx].v_length
                        path.append(Path(v_idx, v_orientation, False, p_length))
                    genome_name_to_idx[g[1]] = len(self.genomes)
                    self.genomes.append(Genome(path))
        # save dictionaries to .json files
        with open(f'vertex_name_to_idx/{self.name}.json', 'w') as f:
            json.dump(vertex_name_to_idx, f)
        with open(f'genome_name_to_idx/{self.name}.json', 'w') as f:
            json.dump(genome_name_to_idx, f)
          
    def add_vertex(self, line):
        v = line.strip().split()
        if len(v)==3:
            v_name = v[1]
            self.vertices.append(Vertex(len(v[2])))
            v_sequence = v[2].upper()+'\n'
            return v_name, v_sequence
        return None, None

    def find_seeds(self, v_idx, PARAM_a):
        '''
        Function finds all occurrences of vertex v_idx, which have not been used before.
        If such occurrences exist, function returns a tuple containing:
            - seeds of collinear walks - list of CollinearWalk objects, describing the occurrence, 
            along with its genome and its orientation relative to the carrying path,
            - carrying path seed orientation - integer (1 or -1), describing the (BAZWZGLĘDNA) orientation of  <--- TO FIX
            the first seed in the list (carrying path seed).
        '''
        v = self.vertices[v_idx]
        collinear_seeds = [] # list to store occurrences of v not used before
        carrying_seed_orientation = None
        if len(v.occurrences) <= PARAM_a: # pruning, if necessary
            for g_idx, i in v.occurrences: # occurrences of the vertex
                genome = self.genomes[g_idx]
                g_path_pos = genome.path[i] # vertex index: int, orientation: int, used: bool
                assert g_path_pos.vertex==v_idx, f'v_idx should be saved in path. Got {v_idx=}, path_v_idx={g_path_pos.vertex}.'
                if g_path_pos.used==False:
                    if not collinear_seeds:
                        orient = 1
                        carrying_seed_orientation = g_path_pos.orientation
                    else:
                        orient = g_path_pos.orientation*carrying_seed_orientation
                    collinear_seeds.append(CollinearWalk(g_idx, i, i, orient))
                    # walk_start_end_check(collinear_seeds[-1], len(genome.path))
            if carrying_seed_orientation is None:
                assert not collinear_seeds, f'If carrying_seed_orientation is None, collinear_seeds must be empty. Got {collinear_seeds=}.'

        return collinear_seeds, carrying_seed_orientation