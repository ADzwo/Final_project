import os
import re
import json
import random
from collections import namedtuple
import pandas as pd
import numpy as np
from datetime import date
SRC = '/root/agesia/magisterka/'

PARAM_m = 50
PARAM_b = 200
SORT_SEEDS = ['no', 'nr_occurences', 'length'][2]
assert SORT_SEEDS in {'no', 'nr_occurences', 'length'}, f'SORT_SEEDS must be one of "no", "nr_occurences", "length", got {SORT_SEEDS}'

CollinearPath = namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
Path = namedtuple('Path', ['vertex', 'orientation', 'used'])
ExtensionVertex = namedtuple('ExtensionVertex', ['coverage', 'distance', 'path_nr', 'start', 'end'])
ShortestPath = namedtuple('ShortestPath', ['vertex', 'orientation'])
Occurence = namedtuple('Occurence', ['genome', 'position'])
# uwaga: gdzieś w genomie można trzymać info o pozycji końcu wierzchołka (łatwo porównywalne z b!) <--- FIX THIS
# Po zakończeniiu dużego kroku: zapamiętać rozszerzenia i z nich korzystać.
# Pamiętać też, które wierzchołki się końćzą w starym t (nowym w0).
class Genome:
    path: list[namedtuple] # Path: vertex index: int, orientation: int, used: bool
    
    def __init__(self, path):
        self.path = path


class Vertex:
    length: int # length of the represented sequence
    occurences: list[namedtuple] # Occurence: genome index, position on its path

    def __init__(self, length, occurences=None):
        self.length = length
        if not occurences:
            self.occurences = []
        else:
            self.occurences = occurences

    def add_occurence(self, g_idx, position):
        self.occurences.append(Occurence(g_idx, position))


class Graph:
    genomes: list # list of Genome objects
    vertices: list # list of Vertex objects

    def __init__(self, graph_path):
        self.vertices = []
        self.genomes = []
        vertex_name_to_idx = {} # dict to save {vertex id from gfa file: index in Graph.vertices}
        genome_name_to_idx = {}
        # find S lines and fill vertex_name_to_idx dict 
        with open(graph_path,'r') as graph_file:
            for line in graph_file.readlines():
                if line.startswith('S'):
                    v = line.strip().split()
                    vertex_name_to_idx[v[1]] = len(self.vertices) # {name: id}
                    self.vertices.append(Vertex(len(v[2])))
        # find P lines and fill genome_name_to_idx dict 
        with open(graph_path,'r') as graph_file:
            for line in graph_file.readlines():
                if line.startswith('P'): # or 'W'
                    g = line.strip().split() # P, name, vertices' names, overlaps
                    path = [] # list[Path] - vertex, orientation, used
                    for v_pos, vertex in enumerate(g[2].split(',')):
                        v_name = vertex[:-1]
                        v_idx = vertex_name_to_idx[v_name]
                        v_orientation = 1 if vertex[-1]=='+' else -1
                        self.vertices[v_idx].add_occurence(len(self.genomes), v_pos) # genome, position
                        path.append(Path(v_idx, v_orientation, False))
                    genome_name_to_idx[g[1]] = len(self.genomes)
                    self.genomes.append(Genome(path))
        # save dictionaries to .json files
        graph_name = re.split(r'[\.\/]', graph_path)[-2]
        with open(f'{SRC}vertex_name_to_idx/{graph_name}.json', 'w') as f:
            json.dump(vertex_name_to_idx, f)
        with open(f'{SRC}genome_name_to_idx/{graph_name}.json', 'w') as f:
            json.dump(genome_name_to_idx, f)

    def find_collinear_blocks(self):
        collinear_blocks = []

        if SORT_SEEDS=='nr_occurences':
            vertex_indices_shuffled = sorted(list(range(len(self.vertices))), key=lambda x: len(self.vertices[x].occurences), reverse=True)
        elif SORT_SEEDS=='length':
            vertex_indices_shuffled = sorted(list(range(len(self.vertices))), key=lambda x: self.vertices[x].length, reverse=True)
        else:
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
                assert _==v_idx, f'v_idx should be saved in path. Got {v_idx=}, path_v_idx={_}.'
                if used==False:
                    if not collinear_seeds:
                        orient = 1
                        carrying_seed_orientation = orientation
                    else:
                        orient = orientation*carrying_seed_orientation
                    collinear_seeds.append(CollinearPath(g_idx, pos, pos, orient)) # namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
                    path_start_end_check(collinear_seeds[-1], len(genome.path))
            if not collinear_seeds:
                continue
            # print('NEW SEED')
            new_block = CollinearBlock(v_idx, collinear_seeds, carrying_seed_orientation)
            new_score = 0
            while new_score>=0:
                w0_idx = new_block.carrying_path[-1]
                Q = BlockExtensions(new_block.collinear_paths, self, w0_idx) # collinear_paths, graph, w0
                r = Q.get_carrying_path_extension(self) # r - path from w0 to t (vertex, orientation), where t - index of the most frequently visited vertex;
                if not r:
                    break
                for wi in r:
                    new_block.carrying_path.append(wi.vertex)
                    new_block.carrying_path_orientations.append(wi.orientation)
                    new_block.update_collinear_walks(wi, Q, self)
                    new_score = new_block.scoring_function(self)
                    if new_score>best_score:
                        best_block = new_block
                        best_score = new_score
            
            if best_score>0:
                collinear_blocks.append(best_block)
                mark_vertices_as_used(self, best_block)
                
        return collinear_blocks
                

class CollinearBlock:
    carrying_path: list[int] # vertex index: int
    carrying_path_orientations: list[int] # orientation: int
    collinear_paths: list[namedtuple] # genome index: int,
                                    # start index from genome's path: int,
                                    # end index from genome's path: int,
                                    # orientaion with respect to the path: int
    def __init__(self, v_initial_idx, v_occurences, carrying_path_seed_orientation):
        self.carrying_path = [v_initial_idx] # vertex index
        self.carrying_path_orientations = [carrying_path_seed_orientation] # orientation
        self.collinear_paths = v_occurences # one of the collinear paths starts with the same vertex as the carrying path

    def scoring_function(self, graph):
        score = 0
        carrying_len = len(self.carrying_path)
        for path in self.collinear_paths:
            if path_length(path, graph)>=PARAM_m:
                i = 0 # position on carrying path
                p_idx = 0 # position on path
                q1 = 0
                bubble = 0
                chain = 0
                # 1) calculate length of the hanging end at the beginning of carrying path
                while i<carrying_len and self.carrying_path[i] not in path:
                    q1 += graph.vertices[self.carrying_path[i]].length
                    if q1 > PARAM_b:
                        return -1
                    i += 1
                if i<carrying_len:
                    chain += graph.vertices[self.carrying_path[i]].length
                    p_idx = i
                    i += 1
                
                # 2) calculate length of the chain and the other hanging end
                while i<carrying_len:
                    v_idx = self.carrying_path[i]
                    idx_on_path, bubble_collinear = find_vertex_on_path_after_idx(graph, path, v_idx, p_idx)
                    if idx_on_path is not None:
                        if bubble_collinear>PARAM_b:
                            return -1
                        chain += bubble + graph.vertices[v_idx].length
                        bubble = 0
                    else:
                        bubble += graph.vertices[self.carrying_path[i]].length
                        if bubble>PARAM_b:
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

        # To FIX: Nie patrzeć na extensions, tylko na paths + b.
        # (Inaczej może być tak, że ścieżka się skończy, jej stare Q nie obejmuje wspólnego wierzchołka z carrying path, 
        # ale ten wierzchołek jest osiągalny w obrębie bubble.)
        w = graph.vertices[w_info.vertex]
        w_orient_on_carrying_path = None
        for occurence in w.occurences:
            # 1) search for paths, whose extensions contain the occurence of w
            path_to_extend = None
            occurence_position = None
            for e_idx, extension in enumerate(block_extensions.extensions):
                if extension.genome==occurence.genome and extension.start<=occurence.position<=extension.end:
                    if not path_to_extend:
                        path_to_extend = e_idx
                        occurence_position = occurence.position
                    # 2) if there are multiple paths containing this occurence of w, 
                    # select the one ending further
                    else:
                        if extension.orient==1:
                            if self.collinear_paths[e_idx].end>self.collinear_paths[path_to_extend].end:
                                path_to_extend = e_idx
                                occurence_position = occurence.position
                        else:
                            if self.collinear_paths[e_idx].start<self.collinear_paths[path_to_extend].start:
                                path_to_extend = e_idx
                                occurence_position = occurence.position
            # 3a) if the path is found, extend it till w
            if path_to_extend is not None:
                path = self.collinear_paths[path_to_extend]
                if path.orient==1:
                    self.collinear_paths[path_to_extend] = CollinearPath(path.genome, path.start, occurence_position, path.orient)
                else:
                    self.collinear_paths[path_to_extend] = CollinearPath(path.genome, occurence_position, path.end, path.orient)
                path_start_end_check(self.collinear_paths[path_to_extend], len(graph.genomes[path.genome].path))
            # 3b) if the path is not found, occurence becomes a new collinear path (provided it is not used)
            elif graph.genomes[occurence.genome].path[occurence.position].used==False:
                if w_orient_on_carrying_path is None:
                    for v_idx, v in enumerate(self.carrying_path):
                        if v==w_info.vertex:
                            w_orient_on_carrying_path = self.carrying_path_orientations[v_idx]
                            break
                occurence_orientation = graph.genomes[occurence.genome].path[occurence.position].orientation
                orient = w_orient_on_carrying_path * occurence_orientation
                self.collinear_paths.append(CollinearPath(occurence.genome, occurence.position, occurence.position, orient))
                
                    
class BlockExtensions:
    extensions: list[namedtuple] # b-extensions
    vertices: dict # {vertex index: (coverage, distance, *shortest_walk)}
                    # coverage - number of occurrences in extensions
                    # distance - length of the shortest walk from w_0
                    # shortest_walk (path number, start, end)
    def __init__(self, collinear_paths, graph, w0_idx):
        self.extensions = []
        self.vertices = {}
        
        for path_nr, collinear_path in enumerate(collinear_paths):
            g_idx = collinear_path.genome
            genome = graph.genomes[g_idx]
            genome_length = len(genome.path)
            path_start_end_check(collinear_path, genome_length)
            w0_positions = find_vertex_on_path(graph, collinear_path, w0_idx) # from end or start - depending on orient --- TO FIX; and check if dist<=b
                # including the last position on the path!
                # Można zapamiętać w poprzednim kroku, czy one się przedłużyły do obecnego w0, czyli poprzedniego t
            if collinear_path.orient==1 and collinear_path.end==genome_length-1: # if we've reached the end of the genome
                self.extensions.append(CollinearPath(-1, -1, -1, -1)) # <- how to improve this?
                continue
            elif collinear_path.orient==-1 and collinear_path.start==0:
                self.extensions.append(CollinearPath(-1, -1, -1, -1))
                continue
            else:
                if collinear_path.orient==1: # proximal - the end of extension proximal to the collinear path
                    proximal = collinear_path.end+1
                else:
                    proximal = collinear_path.start-1
                assert 0<=proximal<genome_length, f'0<=proximal<genome_length is not true! {proximal=}, {genome_length=}'
                p_length = 0
                i = proximal
                while genome.path[i].used==False:
                    v_idx = genome.path[i].vertex
                    p_length += graph.vertices[v_idx].length
                    if p_length>=PARAM_b:
                        # this way, we get all extensions of length x+y, 
                        # where x<b and y is the length of the last vertex of the extension
                        break
                    w0_positions_further = [w0_pos for w0_pos in w0_positions if (i-w0_pos)*collinear_path.orient>=0]
                    if w0_positions_further:
                        distances = [path_length(CollinearPath(g_idx, w0_pos, i, 1), graph) for w0_pos in w0_positions_further] # orient and the order of proximal and i doesn't matter
                        w0_position = np.argmin(distances)
                        distance = distances[w0_position]
                        if v_idx not in self.vertices:
                            self.vertices[v_idx] = ExtensionVertex(1, distance, path_nr, min(proximal, i), max(proximal, i))
                            path_start_end_check(self.vertices[v_idx], genome_length)
                        else:
                            coverage = self.vertices[v_idx].coverage + 1
                            if self.vertices[v_idx].distance>distance or self.vertices[v_idx].distance<0:
                                self.vertices[v_idx] = ExtensionVertex(coverage, distance, path_nr, min(proximal, i), max(proximal, i))
                                path_start_end_check(self.vertices[v_idx], genome_length)
                            else:
                                self.vertices[v_idx] = ExtensionVertex(coverage, *(self.vertices[v_idx][1:]))
                    elif v_idx not in self.vertices:
                        self.vertices[v_idx] = ExtensionVertex(1, -1, -1, -1, -1)
                    else:
                        coverage = self.vertices[v_idx].coverage + 1
                        self.vertices[v_idx] = ExtensionVertex(coverage, *(self.vertices[v_idx][1:]))
                    
                    if i in {0, genome_length-1}:
                        break
                    i += collinear_path.orient # going forwards or backwards on the genome
                
                if collinear_path.orient==1:
                    self.extensions.append(CollinearPath(g_idx, proximal, i, collinear_path.orient))
                else:
                    self.extensions.append(CollinearPath(g_idx, i, proximal, collinear_path.orient))
                path_start_end_check(self.extensions[-1], genome_length)

    
    def get_carrying_path_extension(self, graph):
        '''
        Function finds 
         - index t of a vertex with field distance > 0 (reachable from w0) via a genomic walk, visited by the most extensions.
        Function returns 
         - the shortest walk r from w0 to t.
        '''
        highest_coverage = -1
        w0_to_t = None
        for v_idx in self.vertices:
            if self.vertices[v_idx].distance>0 and self.vertices[v_idx].coverage>highest_coverage:
                # check if v is reachable from w0 and if its coverage is greater than the highest one by now
                w0_to_t = self.vertices[v_idx]
                highest_coverage = w0_to_t.coverage
        if w0_to_t is None:
            return []
        
        r = [] # list of consecutive vertices and their orientations
        extension = self.extensions[w0_to_t.path_nr] # namedtuple('CollinearPath', ['genome', 'start', 'end', 'orient'])
        genome = graph.genomes[extension.genome]
        for i in range(w0_to_t.start, w0_to_t.end+1):
            g_path_pos = genome.path[i] # namedtuple('Path', ['vertex', 'orientation', 'used'])
            r.append(ShortestPath(g_path_pos.vertex, g_path_pos.orientation*extension.orient)) # <- not sure if orientation is ok
        return r
    
def is_vertex_on_path(graph, path, w0_idx):
    i = path.end
    genome_path = graph.genomes[path.genome].path
    while i>=0:
        if genome_path[i].vertex==w0_idx:
            return True
        i -= 1
    return False

def find_vertex_on_path(graph, path, w0_idx):
    positions = []
    i = path.end
    genome_path = graph.genomes[path.genome].path
    while i>=0:
        if genome_path[i].vertex==w0_idx:
            positions.append(i)
        i -= 1
    return positions

def find_vertex_on_path_after_idx(graph, path, w0_idx, p_idx):
    i = p_idx
    length = 0
    genome_path = graph.genomes[path.genome].path
    while i<=path.end:
        length += graph.vertices[genome_path[i].vertex].length
        if genome_path[i].vertex==w0_idx:
            return i, length
        i += 1
    return None, None

def mark_vertices_as_used(graph, block):
    for collinear_path in block.collinear_paths:
        g_idx = collinear_path.genome
        genome_path = graph.genomes[g_idx].path
        for i in range(collinear_path.start, collinear_path.end+1):
            g_path_pos = genome_path[i] # Path = namedtuple('Path', ['vertex', 'orientation', 'used'])
            genome_path[i] = Path(*(g_path_pos[:-1]), True)
    print(f'Marked as used: {len(block.collinear_paths)*(collinear_path.end+1-collinear_path.start)}.')

def path_length(path, graph):
        length = 0
        g_idx = path.genome
        genome_path = graph.genomes[g_idx].path
        for i in range(path.start, path.end+1):
            v_idx = genome_path[i].vertex
            vertex = graph.vertices[v_idx]
            length += vertex.length
        return length

def path_start_end_check(path, genome_length):
    assert path.end<genome_length, f'end >= genome_length; {path=}, {genome_length=}'
    assert path.start>=0, f'start < 0; {path=}, {genome_length-1=}'
    assert path.start<=path.end, f'end < start; {path=}'

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

nr_blocks = []
for graph_path in os.listdir(SRC+'data'):
    for SORT_SEEDS in ['no', 'nr_occurences', 'length']:
        print(f'{graph_path.upper()}, {SORT_SEEDS=}')
        g = Graph(SRC+'data/'+graph_path)
        blocks = g.find_collinear_blocks()
        nr_blocks.append(len(blocks))
        save_blocks(blocks, graph_path.split('.')[0])
print(f'Number of blocks for consecutive options:\n{nr_blocks}')

# SORT_SEEDS = 'no'
# for graph_path in os.listdir(SRC+'data'):
#     for i in range(10):
#         print(f'{graph_path.upper()}, {SORT_SEEDS=}')
#         g = Graph(SRC+'data/'+graph_path)
#         blocks = g.find_collinear_blocks()
#         nr_blocks.append(len(blocks))
#         save_blocks(blocks, graph_path.split('.')[0])

# IMPORTANT NOTES
# start, end and orient once again
    # always: start<end
    # if orient==1, we extend end; otherwise - we extend start.
# Do all S lines need to be before P lines in .gfa?
    # If not - modify reading part. <--- DONE

# Scoring function: bubble <= b only for the carrying path. >m only for path.

# Pomysł (nie jak w artykule): update'ować Q dla danej ścieżki po małym kroku, jeśli znalazłam wspólny wierzchołek z carrying path.