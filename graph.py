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

CollinearWalk = namedtuple('CollinearWalk', ['genome', 'start', 'end', 'orient'])
Path = namedtuple('Path', ['vertex', 'orientation', 'used', 'length']) # length --- length of the suffix of the genome, ending at the last character of this vertex
PathFromW0 = namedtuple('ExtensionVertex', ['distance', 'walk_nr', 'start', 'end'])
ShortestWalk = namedtuple('ShortestWalk', ['vertex', 'orientation'])
Occurence = namedtuple('Occurence', ['genome', 'nr_on_path'])
# Po zakończeniiu dużego kroku: zapamiętać rozszerzenia i z nich korzystać. <- FIX THIS
# Pamiętać też, które wierzchołki się końćzą w starym t (nowym w0).
class Genome:
    path: list[namedtuple] # Path: vertex index: int, orientation: int, used: bool
    
    def __init__(self, path):
        self.path = path


class Vertex:
    length: int # length of the represented sequence
    occurences: list[namedtuple] # Occurence: genome index, index on its path

    def __init__(self, length, occurences=None):
        self.length = length
        if not occurences:
            self.occurences = []
        else:
            self.occurences = occurences

    def add_occurence(self, g_idx, nr_on_path):
        self.occurences.append(Occurence(g_idx, nr_on_path))


class Graph:
    genomes: list # list of Genome objects
    vertices: list # list of Vertex objects

    def __init__(self, graph_file_path):
        self.vertices = []
        self.genomes = []
        vertex_name_to_idx = {} # dict to save {vertex id from gfa file: index in Graph.vertices}
        genome_name_to_idx = {}
        # find S lines and fill vertex_name_to_idx dict 
        with open(graph_file_path,'r') as graph_file:
            for line in graph_file.readlines():
                if line.startswith('S'):
                    v = line.strip().split()
                    vertex_name_to_idx[v[1]] = len(self.vertices) # {name: id}
                    self.vertices.append(Vertex(len(v[2])))
        # find P lines and fill genome_name_to_idx dict 
        with open(graph_file_path,'r') as graph_file:
            for line in graph_file.readlines():
                if line.startswith('P'): # or 'W'
                    g = line.strip().split() # P, name, vertices' names, overlaps
                    path = [] # list[Path] - vertex, orientation, used
                    length = 0
                    for v_pos, vertex in enumerate(g[2].split(',')):
                        v_name = vertex[:-1]
                        v_idx = vertex_name_to_idx[v_name]
                        v_orientation = 1 if vertex[-1]=='+' else -1
                        self.vertices[v_idx].add_occurence(len(self.genomes), v_pos) # genome, nr_on_path
                        length += self.vertices[v_idx].length
                        path.append(Path(v_idx, v_orientation, False, length))
                    genome_name_to_idx[g[1]] = len(self.genomes)
                    self.genomes.append(Genome(path))
        # save dictionaries to .json files
        graph_name = re.split(r'[\.\/]', graph_file_path)[-2]
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
            for g_idx, i in v.occurences: # occurences of the vertex
                genome = self.genomes[g_idx]
                g_path_pos = genome.path[i] # vertex index: int, orientation: int, used: bool
                assert g_path_pos.vertex==v_idx, f'v_idx should be saved in path. Got {v_idx=}, path_v_idx={g_path_pos.vertex}.'
                if g_path_pos.used==False:
                    if not collinear_seeds:
                        orient = 1
                        carrying_seed_orientation = g_path_pos.orientation
                    else:
                        orient = g_path_pos.orientation*carrying_seed_orientation
                    collinear_seeds.append(CollinearWalk(g_idx, i, i, orient)) # namedtuple('CollinearWalk', ['genome', 'start', 'end', 'orient'])
                    walk_start_end_check(collinear_seeds[-1], len(genome.path))
            if not collinear_seeds:
                continue
            print('NEW SEED')
            new_block = CollinearBlock(v_idx, collinear_seeds, carrying_seed_orientation)
            new_score = 0
            while new_score>=0:
                w0_idx = new_block.carrying_path[-1]
                Q = BlockExtensions(new_block.collinear_walks, self, w0_idx) # collinear_walks, graph, w0
                r = Q.get_carrying_path_extension(self, new_block.collinear_walks) # r - shortest walk from w0 to t (vertex, orientation), where t - index of the most frequently visited vertex;
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
    collinear_walks: list[namedtuple] # genome index: int,
                                    # start index from genome's path: int,
                                    # end index from genome's path: int,
                                    # orientaion with respect to the path: int
    def __init__(self, v_initial_idx, v_occurences, carrying_path_seed_orientation):
        self.carrying_path = [v_initial_idx] # vertex index
        self.carrying_path_orientations = [carrying_path_seed_orientation] # orientation
        self.collinear_walks = v_occurences # one of the collinear walks starts with the same vertex as the carrying path

    def scoring_function(self, graph):
        score = 0
        carrying_len = len(self.carrying_path)
        for walk in self.collinear_walks:
            if walk_length(walk, graph)>=PARAM_m:
                i = 0 # nr_on_path on carrying path
                p_idx = 0 # nr_on_path on path
                q1 = 0
                bubble = 0
                chain = 0
                # 1) calculate length of the hanging end at the beginning of carrying path
                while i<carrying_len and self.carrying_path[i] not in walk:
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
                    idx_on_walk = find_vertex_on_path(graph, walk, v_idx, proximal=p_idx)
                    if idx_on_walk is not None:
                        bubble_collinear = walk_length(walk, graph, start=p_idx, end=idx_on_walk)
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
        w = graph.vertices[w_info.vertex]
        w_orient_on_carrying_path = None
        for occurence in w.occurences:
            # 1) search for walks, whose extensions contain the occurence of w
            walk_to_extend = None
            occurence_nr_on_path = None
            for e_idx, extension in enumerate(block_extensions.extensions):
                if extension.genome==occurence.genome and occurence.nr_on_path in range(extension.start, extension.end+1):
                    if walk_to_extend is None:
                        walk_to_extend = e_idx
                        occurence_nr_on_path = occurence.nr_on_path
                    # 2) if there are multiple walks containing this occurence of w, 
                    # select the one ending further
                    else:
                        if extension.orient==1:
                            if self.collinear_walks[e_idx].end>self.collinear_walks[walk_to_extend].end:
                                walk_to_extend = e_idx
                                occurence_nr_on_path = occurence.nr_on_path
                        else:
                            if self.collinear_walks[e_idx].start<self.collinear_walks[walk_to_extend].start:
                                walk_to_extend = e_idx
                                occurence_nr_on_path = occurence.nr_on_path
            # 3a) if the path is found, extend it till w AND update the extension!
            if walk_to_extend is not None:
                walk = self.collinear_walks[walk_to_extend]
                ######################
                if walk.orient==1:
                    self.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, walk.start, occurence_nr_on_path, walk.orient)
                else:
                    self.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, occurence_nr_on_path, walk.end, walk.orient)
                g_len = len(graph.genomes[walk.genome].path)
                walk_start_end_check(self.collinear_walks[walk_to_extend], g_len)
                distal = occurence_nr_on_path + walk.orient
                length = 0
                if distal>=g_len or distal<0:
                    block_extensions.extensions[walk_to_extend] = CollinearWalk(-1, -1, -1, -1)
                else:
                    while distal<g_len-1 and distal>0:
                        length += graph.vertices[distal].length
                        if length>PARAM_b:
                            break # distal is at most in distance PARAM_b from the collinear walk
                        distal += walk.orient
                    if walk.orient==1:
                        block_extensions.extensions[walk_to_extend] = CollinearWalk(walk.genome, occurence_nr_on_path+1, distal, walk.orient)
                    else:
                        block_extensions.extensions[walk_to_extend] = CollinearWalk(walk.genome, distal, walk.end-1, walk.orient)
                    walk_start_end_check(block_extensions.extensions[walk_to_extend], g_len)
                # ######################
                # if walk.orient==1:
                #     self.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, walk.start, occurence_nr_on_path, walk.orient)
                #     distal = occurence_nr_on_path + walk.orient
                #     length = 0
                #     g_len = len(graph.genomes[walk.genome].path)
                #     if distal>=g_len:
                #         block_extensions.extensions[walk_to_extend] = CollinearWalk(-1, -1, -1, -1)
                #     else:
                #         while length<=PARAM_b and distal<g_len-1:
                #             length += graph.vertices[distal].length
                #             distal += walk.orient
                #         block_extensions.extensions[walk_to_extend] = CollinearWalk(walk.genome, occurence_nr_on_path+1, distal, walk.orient)
                # else:
                #     self.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, occurence_nr_on_path, walk.end, walk.orient)
                #     distal = occurence_nr_on_path + walk.orient
                #     length = 0
                #     if distal<0:
                #         block_extensions.extensions[walk_to_extend] = CollinearWalk(-1, -1, -1, -1)
                #     else:
                #         while length<=PARAM_b and distal>=0:
                #             length += graph.vertices[distal].length
                #             distal += walk.orient
                #         block_extensions.extensions[walk_to_extend] = CollinearWalk(walk.genome, distal, walk.end-1, walk.orient)
                # walk_start_end_check(self.collinear_walks[walk_to_extend], len(graph.genomes[walk.genome].path))
                # ######################
                
            # 3b) if such walk is not found, occurence becomes a new collinear path (provided it is not used)
            elif graph.genomes[occurence.genome].path[occurence.nr_on_path].used==False:
                if w_orient_on_carrying_path is None:
                    for v_idx, v in enumerate(self.carrying_path):
                        if v==w_info.vertex:
                            w_orient_on_carrying_path = self.carrying_path_orientations[v_idx]
                            break
                occurence_orientation = graph.genomes[occurence.genome].path[occurence.nr_on_path].orientation
                orient = w_orient_on_carrying_path * occurence_orientation
                self.collinear_walks.append(CollinearWalk(occurence.genome, occurence.nr_on_path, occurence.nr_on_path, orient))
                
                    
class BlockExtensions:
    extensions: list[namedtuple] # list of extensions of type CollinearWalk
    coverage: dict # {vertex index: coverage}, where coverage - number of occurrences in extensions
    shortest_walk: dict  # {vertex index: (distance, *shortest_walk)}
                    # distance - length of the shortest walk from w_0
                    # shortest_walk (walk number, start, end)
    def __init__(self, collinear_walks, graph, w0_idx):
        self.extensions = []
        self.coverage = {}
        self.shortest_walk = {}
        
        for walk_nr, collinear_walk in enumerate(collinear_walks):
            g_idx = collinear_walk.genome
            genome = graph.genomes[g_idx]
            g_len = len(genome.path)
            walk_start_end_check(collinear_walk, g_len)
                # including the last position on the path!
                # Można zapamiętać w poprzednim kroku, czy one się przedłużyły do obecnego w0, czyli poprzedniego t
                # Uwaga: wybrać jedno wystąpienie
            if collinear_walk.orient==1 and collinear_walk.end==g_len-1: # if we've reached the end of the genome
                self.extensions.append(CollinearWalk(-1, -1, -1, -1))
                continue
            elif collinear_walk.orient==-1 and collinear_walk.start==0:
                self.extensions.append(CollinearWalk(-1, -1, -1, -1))
                continue
            else:
                if collinear_walk.orient==1:
                    proximal = collinear_walk.end+1 # proximal - the end of extension proximal to the collinear walk
                    to_search_from = proximal-collinear_walk.orient if proximal>0 else proximal
                else:
                    proximal = collinear_walk.start-1
                    to_search_from = proximal-collinear_walk.orient if proximal<g_len-1 else proximal
                assert 0<=proximal<g_len, f'0<=proximal<genome_length is not true! {proximal=}, {g_len=}'
                assert 0<=to_search_from<g_len, f'0<=to_search_from<genome_length is not true! {to_search_from=}, {g_len=}'
                
                w0_nr_on_path = find_vertex_on_path_till_b(graph, collinear_walk, w0_idx, proximal=to_search_from) # taking the closest position to the collinear walk
                p_length = 0
                i = proximal
                if w0_nr_on_path is None:
                    while genome.path[i].used==False:
                        v_idx = genome.path[i].vertex
                        p_length += graph.vertices[v_idx].length
                        if p_length>PARAM_b: # is it ok? (TO FIX?)
                            i -= collinear_walk.orient
                            break
                        self.coverage[v_idx] = self.coverage.get(v_idx, 0) + 1
                        if i in {0, g_len-1}:
                            break
                        i += collinear_walk.orient
                else:
                    while genome.path[i].used==False:
                        v_idx = genome.path[i].vertex
                        p_length += graph.vertices[v_idx].length
                        if p_length>PARAM_b: # is it ok? (TO FIX?)
                            i -= collinear_walk.orient
                            break
                        self.coverage[v_idx] = self.coverage.get(v_idx, 0) + 1
                        is_w0_before_v = (i-w0_nr_on_path)*collinear_walk.orient>=0
                        if is_w0_before_v:
                            distance = walk_length(CollinearWalk(g_idx, w0_nr_on_path, i, 1), graph) # orient and the order of proximal and i doesn't matter
                            if v_idx not in self.shortest_walk or self.shortest_walk[v_idx].distance>distance:
                                self.shortest_walk[v_idx] = PathFromW0(distance, walk_nr, min(proximal, i), max(proximal, i))
                                walk_start_end_check(self.shortest_walk[v_idx], g_len)
                        if i in {0, g_len-1}:
                            break
                        i += collinear_walk.orient # going forwards or backwards on the genome
            self.extensions.append(CollinearWalk(g_idx, min(proximal, i), max(proximal, i), collinear_walk.orient))

    
    def get_carrying_path_extension(self, graph, collinear_walks):
        '''
        Function finds 
         - index t of a vertex with field distance > 0 (reachable from w0) via a genomic walk, visited by the most extensions.
        Function returns 
         - the shortest walk r from w0 to t.
        '''
        highest_coverage = -1
        w0_to_t = None
        for v_idx in self.shortest_walk: # for all vertices reachable from w0 within distance PARAM_b
            if self.coverage[v_idx]>highest_coverage: # if coverage is greater than the highest one by now
                w0_to_t = self.shortest_walk[v_idx]
                highest_coverage = self.coverage[v_idx]
        if w0_to_t is None:
            return []
        
        r = [] # list of consecutive vertices and their orientations
        walk = collinear_walks[w0_to_t.walk_nr]
        genome = graph.genomes[walk.genome]
        for i in range(w0_to_t.start, w0_to_t.end+1):
            g_path_pos = genome.path[i]
            r.append(ShortestWalk(g_path_pos.vertex, g_path_pos.orientation*walk.orient))
        return r

def find_vertex_on_path(graph:Graph, walk:CollinearWalk, v_to_find:int, proximal=None, distal=None):
    '''
    Function searches for an occurance of vertex with index v_to_find on walk walk.
    The search begins at index proximal and ends at index distal of walk's genome's path.
    Function returns the first index found or None, if no index has been found.
    '''
    if proximal is None:
        proximal = walk.start if walk.orient==1 else walk.end
    if distal is None:
        distal = walk.end if walk.orient==1 else walk.start
    
    genome_path = graph.genomes[walk.genome].path
    for i in range(proximal, distal+1):
        if genome_path[i].vertex==v_to_find:
            return i
    return None

def find_vertex_on_path_till_b(graph:Graph, walk:CollinearWalk, v_to_find:int, proximal=None, distal=None):
    '''
    Function searches for an occurance of vertex with index v_to_find on walk walk.
    The search begins at index proximal and ends at index distal of walk's genome's path.
    Function returns the first index found or None, if no index has been found.
    '''
    if proximal is None:
        proximal = walk.start if walk.orient==1 else walk.end
    if distal is None:
        distal = walk.end if walk.orient==1 else walk.start
    
    genome_path = graph.genomes[walk.genome].path
    length = 0
    for i in range(proximal, distal+1):
        v_idx = genome_path[i].vertex
        if v_idx==v_to_find:
            return i
        length += graph.vertices[v_idx].length
        if length>PARAM_b:
            break
    return None

def mark_vertices_as_used(graph, block):
    '''
    Functions marks vertex occurences of CollinearBlock block as used, 
    i.e. it changes each vertex occurence's parameter used to True.
    '''
    nr_used = 0
    for collinear_walk in block.collinear_walks:
        g_idx = collinear_walk.genome
        genome_path = graph.genomes[g_idx].path
        for i in range(collinear_walk.start, collinear_walk.end+1):
            g_path_pos = genome_path[i]
            genome_path[i] = Path(*(g_path_pos[:-2]), True, g_path_pos[-1])
        nr_used += collinear_walk.end+1-collinear_walk.start
    print(f'Marked as used: {nr_used}.')

def walk_length(walk, graph, start=None, end=None):
        g_idx = walk.genome
        genome_path = graph.genomes[g_idx].path
        if start is None:
            start = walk.start
        if end is None:
            end = walk.end
        if start==0:
            return graph.vertices[genome_path[end].vertex].length
        else:
            return graph.vertices[genome_path[end].vertex].length - graph.vertices[genome_path[start-1].vertex].length

def walk_start_end_check(walk, genome_length):
    assert walk.end<genome_length, f'end >= genome_length; {walk=}, {genome_length=}'
    assert walk.start>=0, f'start < 0; {walk=}, {genome_length-1=}'
    assert walk.start<=walk.end, f'end < start; {walk=}'

def save_blocks(blocks:list[CollinearWalk], graph_name):
    df_all = pd.DataFrame()
    for b_nr, block in enumerate(blocks):
        df = pd.DataFrame.from_records(block.collinear_walks, columns=CollinearWalk._fields)
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
for graph_file_path in os.listdir(SRC+'data'):
    for SORT_SEEDS in ['nr_occurences', 'length', 'no']:
        print(f'{graph_file_path.upper()}, {SORT_SEEDS=}')
        g = Graph(SRC+'data/'+graph_file_path)
        blocks = g.find_collinear_blocks()
        nr_blocks.append(len(blocks))
        save_blocks(blocks, graph_file_path.split('.')[0])
print(f'Number of blocks for consecutive options:\n{nr_blocks}')

# SORT_SEEDS = 'no'
# for graph_file_path in os.listdir(SRC+'data'):
#     for i in range(10):
#         print(f'{graph_file_path.upper()}, {SORT_SEEDS=}')
#         g = Graph(SRC+'data/'+graph_file_path)
#         blocks = g.find_collinear_blocks()
#         nr_blocks.append(len(blocks))
#         save_blocks(blocks, graph_file_path.split('.')[0])

# IMPORTANT NOTES
# start, end and orient once again
    # always: start<end
    # if orient==1, we extend end; otherwise - we extend start.

# Do all S lines need to be before P lines in .gfa?
    # If not - modify reading part. <--- DONE

# Scoring function: bubble <= b only for the carrying path. >m only for path.

# Pomysł (nie jak w artykule): 
    # update'ować Q dla danej ścieżki po małym kroku, 
    # jeśli znalazłam wspólny wierzchołek z carrying path.