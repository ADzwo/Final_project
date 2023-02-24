import json
import random
from collections import namedtuple
src = '/root/agesia/magisterka/'

PARAM_m = 50
PARAM_b = 200

CollinearPath = namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
Path = namedtuple('Path', ['vertex', 'orientation', 'used'])
ExtensionVertex = namedtuple('ExtensionsVertex', ['coverage', 'distance', 'path_nr', 'start', 'end'])


class Genome:
    path: list[tuple] # vertex index: int, orientation: int, used: bool
    
    def __init__(self, path):
        self.path = path


class Vertex:
    length: int # legth of the represented sequence
    occurences: list[tuple] # genome index, position on its path

    def __init__(self, length, occurences=None):
        self.length = length
        if not occurences:
            self.occurences = []
        else:
            self.occurences = occurences
    def add_occurence(self, genome_idx, position):
        self.occurences.append((genome_idx, position))


class Graph:
    genomes: list # list of Genome objects
    vertices: list # list of Vertex objects

    def __init__(self, graph_path):
        self.vertices = []
        self.genomes = []
        vertex_name_to_idx = {} # dict to save {vertex id from gfa file: index in Graph.vertices}
        genome_name_to_idx = {}
        with open(graph_path,'r') as graph_file:
            for line in graph_file.readlines():
                if line.startswith('S'):
                    v = line.strip().split('\t')
                    vertex_name_to_idx[v[1]] = len(self.vertices)
                    self.vertices.append(Vertex(len(v[2])))
                elif line.startswith('P'):
                    g = line.strip().split('\t') # P, name, vertices' names, overlaps
                    path = [] # list(tuple(int, int, bool)) # vertex, orientation, used
                    for v_pos, vertex in enumerate(g[2].split(',')):
                        v_name = vertex[:-1]
                        v_idx = vertex_name_to_idx[v_name]
                        v_orientation = vertex[-1]
                        self.vertices[v_idx].add_occurence(len(self.genomes), v_pos) # genome, position
                        path.append(Path(v_idx, v_orientation, False))
                    genome_name_to_idx[g[1]] = len(self.genomes)
                    self.genomes.append(Genome(path))
        with open('vertex_name_to_idx.json', 'w') as f:
            json.dump(vertex_name_to_idx, f)
        with open('genome_name_to_idx.json', 'w') as f:
            json.dump(genome_name_to_idx, f)
    
    def find_collinear_blocks(self):
        collinear_blocks = []
        vertex_indices_shuffled = random.shuffle(list(range(len(self.vertices))))
        for v_idx in vertex_indices_shuffled: # select a vertex
            v = self.vertices[vertex_indices_shuffled]
            v_occurences = []
            for g_nr, pos in v.occurences: # occurences of the vertex
                genome = self.genomes[g_nr]
                path_v_idx, orientation, used = genome.path[pos] # vertex index: int, orientation: int, used: bool
                assert path_v_idx==v_idx, f'v_idx should be saved in path. Got {v_idx=}, {path_v_idx=}.'
                if not used:
                    if v_occurences:
                        orient = 1 if orientation==v_occurences[0][-1] else -1 # chyyyyba sprawdzam, czy orientacja 
                                                                            # jest taka jak dla carrying path
                    else:
                        orient = orientation
                    v_occurences.append(CollinearPath(g_nr, path_v_idx, path_v_idx+1, orient))
            if v_occurences:
                new_block = CollinearBlock(v_idx, v_occurences)

            while new_block.scoring_function(self)>=0:
                w0 = new_block.carrying_path[-1]
                Q = BlockExtensions(new_block.collinear_paths, w0)
                t, r = Q.most_visited_vertex() # t - index of the most frequently visited vertex; r - path from w0 to t (genome, start, end, orient)
                for i in range(len(r)):
                    new_block.carrying_path.append(r[i].vertex)
                    new_block.carrying_path_orientations.append(r[i].orient)
                    Q.update_collinear_walks(r[i])
                
        pass
                

class CollinearBlock:
    carrying_path: list[tuple] # vertex index: int, orientation: int
    collinear_paths: list[namedtuple] # genome index: int,
                                # start index from genome's path: int,
                                # end index from genome's path: int,
                                # orientaion with respect to the path: int
    def __init__(self, v_initial_idx, v_occurences):
        # self.carrying_path = [(v_initial_idx, v_occurences[0][-1])] # vertex index, orientation
        self.carrying_path = [v_initial_idx] # vertex index
        self.carrying_path_orientations = [v_occurences[0][-1]] # orientation
        self.collinear_paths = v_occurences # one of the collinear paths starts with the same vertex as the carrying path

    def scoring_function(self, graph):
        score = 0
        for path in self.collinear_paths:
            if path.end - path.start >= PARAM_m:
                genome = graph.genomes[path.genome]
                i = path.start
                q1 = 0
                bubble = 0
                chain = 0
                while genome.path[i].vertex not in self.carrying_path:
                    q1 += 1 # calculate length of the hanging end at the beginning of the path
                    if i==path.end:
                        break
                    if q1 > PARAM_b:
                        return -1
                    i += 1
                else:
                    chain += 1
                
                while i<=path.end:
                    if genome.path[i].vertex in self.carrying_path:
                        chain += bubble + 1
                        bubble = 0
                    else:
                        bubble += 1
                        if bubble > PARAM_b:
                            return -1 # instead of -infty (from SibeliaZ algorithm)
                score += chain - (q1 + bubble)**2

        return score

                    
class BlockExtensions:
    extensions: list[namedtuple] # b-extensions
    vertices: dict # {vertex index: (coverage, distance, *shortest_walk)}
                    # coverage - number of occurrences in extensions
                    # distance - length of the shortest walk from w_0
                    # shortest_walk (path number, range)
    def __init__(self, collinear_paths, graph, w0):
        self.extensions = []
        self.vertices = {}
        for path_nr, path in enumerate(collinear_paths):
            genome_path_len = len(path.genome.path)
            if path.end >= genome_path_len-1: # if we've reached the end of the genome
                self.extensions.append(CollinearPath(-1, -1, -1, -1)) # <- to improve
            else:
                end_idx = min(path.end+PARAM_b, len(genome_path_len)-1)
                self.extensions.append(CollinearPath(path.genome, path.end+1, end_idx, path.orient))
                for i in range(path.end+1, end_idx+1):
                    v_idx = path.genome
                    if v_idx not in self.vertices:
                        self.vertices[v_idx] = ExtensionVertex(1, i-path.end, path_nr, path.end, i)
                    else:
                        self.vertices[v_idx].coverage += 1
                        if self.vertices[v_idx].distance > i-path.end and graph[path.genome].path[path.end-1].vertex==w0:
                            self.vertices[v_idx].distance = i-path.end
                            self.vertices[v_idx].path_nr = path_nr
                            self.vertices[v_idx].start = path.end
                            self.vertices[v_idx].end = i
    
    def most_visited_vertex(self, graph):
        '''
        Function returns 
         - vertex t reachable from w0 via a genomic walk, visited by the most extensions
         - the shortest walk r from w0 to t.
        '''
        t = None
        highest_coverage = 0
        for v in self.vertices:
            if self.vertices[v].distance > 0 and (self.vertices[v].coverage > highest_coverage or t is None):
                # check if v is reachable from w0 and if its coverage is greater than the highest one by now
                highest_coverage = self.vertices[v].coverage
                t = v
        r = [] # list of consecutive vertices and their orientations
        path = self.extensions[self.vertices[t].path_nr] # namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
        genome = graph.genomes[path.genome]
        for i in range(path.start, path.end+1):
            genome_path_element = genome.path[i] # namedtuple('Path', ['vertex', 'orientation', 'used'])
            r.append((genome_path_element.vertex, genome_path_element.orientation*path.orient)) # <- not sure which orientation to use
        return t, r


    def update_collinear_walks(vertex_info:tuple, graph): # vertex_info - tuple(vertex index, orientation)
        # for occurence in graph.vertices(vertex_info.vertex).occurences:
        pass

Graph(src+'data/chr13_KI270842v1_alt_10_5_2.gfa')
