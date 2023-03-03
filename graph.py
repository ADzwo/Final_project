import os
import re
import json
import random
from collections import namedtuple
import pandas as pd
from datetime import date
SRC = '/root/agesia/magisterka/'

PARAM_m = 50
PARAM_b = 200
SORT_SEEDS = ['no', 'nr_occurences', 'length'][2]

CollinearPath = namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
Path = namedtuple('Path', ['vertex', 'orientation', 'used'])
ExtensionVertex = namedtuple('ExtensionVertex', ['coverage', 'distance', 'path_nr', 'start', 'end'])
ShortestPath = namedtuple('ShortestPath', ['vertex', 'orientation'])
Occurence = namedtuple('Occurence', ['genome', 'position'])

class Genome:
    path: list[tuple] # vertex index: int, orientation: int, used: bool
    
    def __init__(self, path):
        self.path = path


class Vertex:
    length: int # legth of the represented sequence
    occurences: list[namedtuple] # genome index, position on its path

    def __init__(self, length, occurences=None):
        self.length = length
        if not occurences:
            self.occurences = []
        else:
            self.occurences = occurences

    def add_occurence(self, genome_idx, position):
        self.occurences.append(Occurence(genome_idx, position))


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
                elif line.startswith('P'): # or 'W'
                    g = line.strip().split('\t') # P, name, vertices' names, overlaps
                    path = [] # list(tuple(int, int, bool)) # vertex, orientation, used
                    for v_pos, vertex in enumerate(g[2].split(',')):
                        v_name = vertex[:-1]
                        v_idx = vertex_name_to_idx[v_name]
                        v_orientation = 1 if vertex[-1]=='+' else -1
                        self.vertices[v_idx].add_occurence(len(self.genomes), v_pos) # genome, position
                        path.append(Path(v_idx, v_orientation, False))
                    genome_name_to_idx[g[1]] = len(self.genomes)
                    self.genomes.append(Genome(path))
        graph_name = re.split(r'[\.\/]', graph_path)[-2]
        with open(f'{SRC}vertex_name_to_idx/{graph_name}.json', 'w') as f:
            json.dump(vertex_name_to_idx, f)
        with open(f'{SRC}genome_name_to_idx/{graph_name}.json', 'w') as f:
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
            genome_id = collinear_path.genome
            genome_path = self.genomes[genome_id].path
            for i in range(collinear_path.start, collinear_path.end+1):
                path = genome_path[i] # Path = namedtuple('Path', ['vertex', 'orientation', 'used'])
                genome_path[i] = Path(path.vertex, path.orientation, True)
        print(f'Marked as used: {len(block.collinear_paths)*(collinear_path.end+1-collinear_path.start)}.')

    def find_collinear_blocks(self):
        collinear_blocks = []

        if SORT_SEEDS=='nr_occurences':
            vertex_indices_shuffled = sorted(list(range(len(self.vertices))), key=lambda x: len(self.vertices[x].occurences))
        elif SORT_SEEDS=='length':
            vertex_indices_shuffled = sorted(list(range(len(self.vertices))), key=lambda x: len(self.vertices[x].occurences))
        else:
            assert SORT_SEEDS=='no', f'SORT_SEEDS must be one of "no", "nr_occurences", "length", got {SORT_SEEDS}'
            vertex_indices_shuffled = list(range(len(self.vertices)))
            random.shuffle(vertex_indices_shuffled)
        
        for v_idx in vertex_indices_shuffled: # select a vertex --- seed of a new CollinearBlock
            best_score = -1
            best_block = None
            v = self.vertices[v_idx]
            collinear_seeds = [] # list to store occurences of v not used before
            for g_idx, pos in v.occurences: # occurences of the vertex
                genome = self.genomes[g_idx]
                _, orientation, used = genome.path[pos] # vertex index: int, orientation: int, used: bool
                # assert _==v_idx, f'v_idx should be saved in path. Got {v_idx=}, path_v_idx={_}.'
                if used==False:
                    if not collinear_seeds:
                        orient = 1
                        carrying_seed_orientation = orientation
                    else:
                        orient = orientation*carrying_seed_orientation 
                    collinear_seeds.append(CollinearPath(g_idx, pos, pos, orient)) # namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
                    assert collinear_seeds[-1].start>=0 and collinear_seeds[-1].end<len(genome.path), f'start={collinear_seeds[-1].start}, end={collinear_seeds[-1].end}, genome len={len(genome.path)}'
            if not collinear_seeds:
                break
            print('NEW SEED')
            new_block = CollinearBlock(v_idx, collinear_seeds, carrying_seed_orientation)
            new_score = 0
            while new_score>=0:
                w0 = new_block.carrying_path[-1]
                Q = BlockExtensions(new_block.collinear_paths, self, w0) # collinear_paths, graph, w0
                r = Q.extend_carrying_path(self) # t - index of the most frequently visited vertex; r - path from w0 to t (vertex, orientation)
                if not r:
                    break
                for i in range(len(r)):
                    new_block.carrying_path.append(r[i].vertex)
                    new_block.carrying_path_orientations.append(r[i].orientation)
                    new_block.update_collinear_walks(r[i], Q, self)
                    new_score = new_block.scoring_function(self)
                    if new_score>best_score:
                        best_block = new_block
                        # print(f'{len(best_block.collinear_paths)=}')
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
            if abs(path.end - path.start) >= PARAM_m:
                genome = graph.genomes[path.genome]
                i = path.start
                q1 = 0
                while i<=path.end and genome.path[i].vertex not in self.carrying_path:
                    q1 += 1 # calculate length of the hanging end at the beginning of the path
                    if i==path.end:
                        chain += 1
                        break
                    if q1 > PARAM_b:
                        return -1
                    i += 1
                
                bubble = 0
                chain = 0
                while i<=path.end:
                    if genome.path[i].vertex in self.carrying_path:
                        chain += bubble + 1
                        bubble = 0
                    else:
                        bubble += 1
                        if bubble > PARAM_b:
                            return -1 # instead of -infty (from SibeliaZ algorithm)
                    i += 1
                score += chain - (q1 + bubble)**2
        return score

    def update_collinear_walks(self, w_info, block_extensions, graph): # vertex_info - tuple(vertex index, orientation)
        
        # w bloku szukam ścieżek, których przedłużenia zawierają to wystąpienie w
            # extension: namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
            # occurence: namedtuple('Occurence', ['genome', 'position'])
        # wybieram ścieżkę, która kończy się najpóźniej
            # jeśli znalazłam taką ścieżkę, to przedłużam ją o jej przedłużenie obcięte do w
            # jeśli nie znalazłam, to dodaję to wystąpienie do bloku jako nową ścieżkę

        w = graph.vertices[w_info.vertex]
        w_orient_on_carrying_path = None
        for occurence in w.occurences:
            path_to_extend = None
            for e_idx, extension in enumerate(block_extensions.extensions):
                if extension.genome==occurence.genome and extension.start<=occurence.position<=extension.end:
                    if not path_to_extend:
                        path_to_extend = e_idx
                    else: # there can be multiple paths containing this occurence of w
                        if extension.orient==1:
                            if extension.end>block_extensions.extensions[path_to_extend].end:
                                path_to_extend = e_idx
                        else:
                            if extension.start<block_extensions.extensions[path_to_extend].start:
                                path_to_extend = e_idx
            if path_to_extend is not None:
                path = self.collinear_paths[path_to_extend]
                if path.orient==1:
                    self.collinear_paths[path_to_extend] = CollinearPath(path.genome, path.start, block_extensions.extensions[path_to_extend].end, path.orient)
                else:
                    self.collinear_paths[path_to_extend] = CollinearPath(path.genome, block_extensions.extensions[path_to_extend].start, path.end, path.orient)
            elif graph.genomes[occurence.genome].path[occurence.position].used==False:
                if w_orient_on_carrying_path is None:
                    for v_idx, v in enumerate(self.carrying_path):
                        if v==w_info.vertex:
                            w_orient_on_carrying_path = self.carrying_path_orientations[v_idx]
                            break
                occurence_orient = graph.genomes[occurence.genome].path[occurence.position].orientation
                orient = w_orient_on_carrying_path* occurence_orient
                self.collinear_paths.append(CollinearPath(occurence.genome, occurence.position, occurence.position, orient))
                
        
        # # ZAŚLEPKA <- dodaję całe przedłużenie naraz
        # for i in range(len(self.collinear_paths)):
        #     path = self.collinear_paths[i] # namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
        #     if path.orient==1:
        #         self.collinear_paths[i] = CollinearPath(path.genome, path.start, block_extensions.extensions[i].end, path.orient)
        #     else:
        #         self.collinear_paths[i] = CollinearPath(path.genome, block_extensions.extensions[i].start, path.end, path.orient)


                    
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
            g_idx = collinear_path.genome
            genome = graph.genomes[g_idx]
            genome_path_len = len(genome.path)
            if collinear_path.end>genome_path_len-1 or collinear_path.start<0:
                    raise ValueError(f'collinear_path.end > genome_path_len-1, {collinear_path.start=}, {collinear_path.end=}, {genome_path_len-1=}. {collinear_path.orient=}')
            if collinear_path.orient==1 and collinear_path.end==genome_path_len-1: # if we've reached the end of the genome
                self.extensions.append(CollinearPath(-1, -1, -1, -1)) # <- to improve
            elif collinear_path.orient==-1 and collinear_path.start==0:
                self.extensions.append(CollinearPath(-1, -1, -1, -1))
            else:
                if collinear_path.orient==1:
                    distal = min(collinear_path.end+PARAM_b, genome_path_len-1)
                    proximal = collinear_path.end+1
                else:
                    proximal = collinear_path.start-1 # start is closer to the collinear path than end
                    distal = max(collinear_path.start-PARAM_b, 0)
                for i in range(proximal, distal+1):
                    if genome.path[i].used:
                        distal = i
                        break
                    v_idx = genome.path[i].vertex
                    distance = abs(i-proximal+1)
                    if v_idx not in self.vertices:
                        if graph.is_vertex_on_path(collinear_path, w0):
                            if collinear_path.orient==1:
                                self.vertices[v_idx] = ExtensionVertex(1, distance, path_nr, proximal, i) # namedtuple('ExtensionVertex', ['coverage', 'distance', 'path_nr', 'start', 'end'])
                            else:
                                self.vertices[v_idx] = ExtensionVertex(1, distance, path_nr, i, proximal)
                        else:
                            self.vertices[v_idx] = ExtensionVertex(1, -1, -1, -1, -1)
                    else:
                        coverage = self.vertices[v_idx].coverage + 1
                        if self.vertices[v_idx].distance > distance and graph.is_vertex_on_path(collinear_path, w0):
                            if collinear_path.orient==1:
                                self.vertices[v_idx] = ExtensionVertex(coverage, distance, path_nr, proximal, i)
                            else:
                                self.vertices[v_idx] = ExtensionVertex(coverage, distance, path_nr, i, proximal)
                        else:
                            self.vertices[v_idx] = ExtensionVertex(coverage, self.vertices[v_idx].distance, self.vertices[v_idx].path_nr,
                                                                   self.vertices[v_idx].start, self.vertices[v_idx].end)
                if collinear_path.orient==1:
                    assert proximal<=distal, f'{proximal=}, {distal=}'
                    self.extensions.append(CollinearPath(g_idx, proximal, distal, collinear_path.orient))
                else:
                    assert proximal>=distal, f'{proximal=}, {distal=}'
                    self.extensions.append(CollinearPath(g_idx, distal, proximal, collinear_path.orient))
                assert self.extensions[-1].start>=0
                assert self.extensions[-1].end<genome_path_len
    
    def extend_carrying_path(self, graph):
        '''
        Function finds 
         - index t of a vertex reachable from w0 via a genomic walk, visited by the most extensions.
        Function returns 
         - the shortest walk r from w0 to t.
        '''
        t = None
        highest_coverage = -1
        for v_idx in self.vertices:
            if self.vertices[v_idx].distance>0 and self.vertices[v_idx].coverage>highest_coverage:
                # check if v is reachable from w0 and if its coverage is greater than the highest one by now
                highest_coverage = self.vertices[v_idx].coverage
                t = v_idx
        if t is None:
            return []
        
        r = [] # list of consecutive vertices and their orientations
        extension = self.extensions[self.vertices[t].path_nr] # namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
        genome = graph.genomes[extension.genome]
        for i in range(extension.start, extension.end+1):
            genome_path_element = genome.path[i] # namedtuple('Path', ['vertex', 'orientation', 'used'])
            r.append(ShortestPath(genome_path_element.vertex, genome_path_element.orientation*extension.orient)) # <- not sure which orientation to use
        return r


def save_blocks(blocks:list[CollinearPath], graph_name):
    df_all = pd.DataFrame()
    for b_nr, block in enumerate(blocks):
        df = pd.DataFrame.from_records(block.collinear_paths, columns=CollinearPath._fields)
        df['block_nr'] = b_nr
        df_all = pd.concat([df_all, df])
    today = str(date.today()).replace('-', '_')
    
    if SORT_SEEDS=='nr_occurences':
        name = f'{graph_name}_{today}_sort_by_nr_occurences.csv'
    elif SORT_SEEDS=='length':
        name = f'{graph_name}_{today}_sort_by_length.csv'
    else:
        name = f'{graph_name}_{today}.csv'
    df_all.to_csv(f'blocks/{name}', index=False)

for graph_path in os.listdir(SRC+'data'):
    for SORT_SEEDS in ['no', 'nr_occurences', 'length']:
        print(f'{graph_path.upper()}, {SORT_SEEDS=}')
        g = Graph(SRC+'data/'+graph_path)
        blocks = g.find_collinear_blocks()
        save_blocks(blocks, graph_path.split('.')[0])




# graph_path = 'chr13_KI270842v1_alt_10_5_2.gfa'
# g = Graph(src+'data/'+graph_path)
# blocks = g.find_collinear_blocks()
# save_blocks(blocks)

# IMPORTANT NOTES
# start, end and orient once again
    # always: start<end
    # if orient==1, we extend end; otherwise - we extend start.