import json
import random
from collections import namedtuple
src = '/root/agesia/magisterka/'

PARAM_m = 50
PARAM_b = 200

CollinearPath = namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
Path = namedtuple('Path', ['vertex', 'orientation', 'used'])
ExtensionVertex = namedtuple('ExtensionVertex', ['coverage', 'distance', 'path_nr', 'start', 'end'])
ShortestPath = namedtuple('ShortestPath', ['vertex', 'orientation'])

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

    def is_vertex_on_path(self, path, w0):
        idx = path.end
        genome_path = self.genomes[path.genome].path
        while idx>=0:
            if genome_path[idx].vertex==w0:
                return True
            idx -= 1
        return False
    
    def mark_vertices_as_used(self, block):
        for collinear_path in block.collinear_paths:
            genome_path = self.genomes[collinear_path.genome].path
            for i in range(collinear_path.start, collinear_path.end+1):
                v_idx = genome_path[i].vertex
                self.vertices[v_idx].used = True

    def find_collinear_blocks(self):
        collinear_blocks = []
        vertex_indices_shuffled = random.shuffle(list(range(len(self.vertices))))
        
        for v_idx in vertex_indices_shuffled: # select a vertex --- seed of a new CollinearBlock
            best_score = -1
            best_block = None
            v = self.vertices[v_idx]
            v_occurences = [] # list to store occurences of v not used before
            for g_nr, pos in v.occurences: # occurences of the vertex
                genome = self.genomes[g_nr]
                path_v_idx, orientation, used = genome.path[pos] # vertex index: int, orientation: int, used: bool
                assert path_v_idx==v_idx, f'v_idx should be saved in path. Got {v_idx=}, {path_v_idx=}.'
                if not used:
                    if v_occurences:
                        orient = 1 if orientation==v_occurences[0][-1] else -1 # chyyyyba sprawdzam, czy orientacja 
                                                                            # jest taka jak dla carrying path
                    else:
                        orient = 1
                        carrying_path_seed_orientation = orientation

                    v_occurences.append(CollinearPath(g_nr, path_v_idx, path_v_idx+1, orient))
            if v_occurences:
                new_block = CollinearBlock(v_idx, v_occurences, carrying_path_seed_orientation)

            while new_block.scoring_function(self)>=0:
                w0 = new_block.carrying_path[-1]
                Q = BlockExtensions(new_block.collinear_paths, w0)
                t, r = Q.most_visited_vertex() # t - index of the most frequently visited vertex; r - path from w0 to t (vertex, orientation)
                for i in range(len(r)):
                    new_block.carrying_path.append(r[i].vertex)
                    new_block.carrying_path_orientations.append(r[i].orientation)
                    Q.update_collinear_walks(r[i])
                    new_score = new_block.scoring_function(self)
                    if new_score>best_score:
                        best_block = new_block
                        best_score = new_score
            
            if best_score>0:
                collinear_blocks.append(best_block)
                self.mark_vertices_as_used(best_block)
                
        return collinear_blocks
                

class CollinearBlock:
    carrying_path: list[int] # vertex index: int
    carrying_path_orientations: list[int] # orientation: int
    collinear_paths: list[namedtuple] # genome index: int,
                                    # start index from genome's path: int,
                                    # end index from genome's path: int,
                                    # orientaion with respect to the path: int
    def __init__(self, v_initial_idx, v_occurences, carrying_path_seed_orientation):
        # self.carrying_path = [(v_initial_idx, v_occurences[0][-1])] # vertex index, orientation
        self.carrying_path = [v_initial_idx] # vertex index
        self.carrying_path_orientations = [carrying_path_seed_orientation] # orientation
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
                    # shortest_walk (path number, start, end)
    def __init__(self, collinear_paths, graph, w0):
        self.extensions = []
        self.vertices = {}
        for path_nr, collinear_path in enumerate(collinear_paths):
            genome_id = collinear_path.genome
            genome = graph.genomes[genome_id]
            genome_path_len = len(genome.path)
            if collinear_path.end >= genome_path_len-1: # if we've reached the end of the genome
                self.extensions.append(CollinearPath(-1, -1, -1, -1)) # <- to improve
            else:
                end = min(collinear_path.end+PARAM_b, genome_path_len-1)
                start = collinear_path.end+1
                # self.extensions.append(CollinearPath(genome_id, start, end, collinear_path.orient))
                for i in range(start, end):
                    if genome.path[i].used:
                        end = i
                        break
                    v_idx = genome.path[i].vertex
                    if v_idx not in self.vertices:
                        if graph.is_vertex_on_path(collinear_path, w0):
                            self.vertices[v_idx] = ExtensionVertex(1, i-start+1, path_nr, start, i) # namedtuple('ExtensionVertex', ['coverage', 'distance', 'path_nr', 'start', 'end'])
                        else:
                            self.vertices[v_idx] = ExtensionVertex(1, -1, -1, -1, -1)
                    else:
                        self.vertices[v_idx].coverage += 1
                        if self.vertices[v_idx].distance > i-start+1:
                            if graph.is_vertex_on_path(collinear_path, w0):
                                self.vertices[v_idx].distance = i-start+1
                                self.vertices[v_idx].path_nr = path_nr
                                self.vertices[v_idx].start = start
                                self.vertices[v_idx].end = i
                self.extensions.append(CollinearPath(genome_id, start, end, collinear_path.orient))
    
    def most_visited_vertex(self, graph):
        '''
        Function returns 
         - index t of a vertex reachable from w0 via a genomic walk, visited by the most extensions
         - the shortest walk r from w0 to t.
        '''
        t = None
        highest_coverage = -1
        for v_idx in self.vertices:
            if self.vertices[v_idx].distance > 0 and self.vertices[v_idx].coverage > highest_coverage:
                # check if v is reachable from w0 and if its coverage is greater than the highest one by now
                highest_coverage = self.vertices[v_idx].coverage
                t = v_idx
        if t is None:
            print('t is None!')
            return None, []
        #####
        r = [] # list of consecutive vertices and their orientations
        extension = self.extensions[self.vertices[t].path_nr] # namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
        genome = graph.genomes[extension.genome]
        for i in range(extension.start, extension.end+1):
            genome_path_element = genome.path[i] # namedtuple('Path', ['vertex', 'orientation', 'used'])
            r.append(ShortestPath(genome_path_element.vertex, genome_path_element.orientation*extension.orient)) # <- not sure which orientation to use
        return t, r


    def update_collinear_walks(vertex_info:tuple, graph): # vertex_info - tuple(vertex index, orientation)
        # for occurence in graph.vertices(vertex_info.vertex).occurences:
        pass

Graph(src+'data/chr13_KI270842v1_alt_10_5_2.gfa')
