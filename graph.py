import os
import re
import math
import json
import random
from collections import namedtuple
import pandas as pd
import numpy as np
from datetime import date

PARAM_m = 50
PARAM_b = 200
SORT_SEEDS = ['no', 'nr_occurrences', 'length'][2]
assert SORT_SEEDS in {'no', 'nr_occurrences', 'length'}, f'SORT_SEEDS must be one of "no", "nr_occurrences", "length", got {SORT_SEEDS}'

CollinearWalk = namedtuple('CollinearWalk', ['genome', 'start', 'end', 'orient'])
Path = namedtuple('Path', ['vertex', 'orientation', 'used', 'length']) # length --- length of the suffix of the genome, ending at the last character of this vertex
PathFromW0 = namedtuple('ExtensionVertex', ['distance', 'walk_nr', 'w0_nr_on_path', 't_nr_on_path'])
CarryingPathExtension = namedtuple('CarryingPathExtension', ['vertex', 'orientation'])
occurrence = namedtuple('occurrence', ['genome', 'nr_on_path'])
Score = namedtuple('Score', ['q1', 'q3', 'idx_g_path', 'idx_c_path']) # last indices used to calculate score_so_far, for genome and carrying path

class Genome:
    path: list[namedtuple] # Path: vertex index: int, orientation: int, used: bool
    
    def __init__(self, path):
        self.path = path


class Vertex:
    length: int # length of the represented sequence
    occurrences: list[namedtuple] # occurrence: genome index, index on its path

    def __init__(self, length, occurrences=None):
        self.length = length
        if not occurrences:
            self.occurrences = []
        else:
            self.occurrences = occurrences

    def add_occurrence(self, g_idx, nr_on_path):
        self.occurrences.append(occurrence(g_idx, nr_on_path))


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
                        self.vertices[v_idx].add_occurrence(len(self.genomes), v_pos) # genome, nr_on_path
                        length += self.vertices[v_idx].length
                        path.append(Path(v_idx, v_orientation, False, length))
                    genome_name_to_idx[g[1]] = len(self.genomes)
                    self.genomes.append(Genome(path))
        # save dictionaries to .json files
        graph_name = re.split(r'[\.\/]', graph_file_path)[-2]
        with open(f'vertex_name_to_idx/{graph_name}.json', 'w') as f:
            json.dump(vertex_name_to_idx, f)
        with open(f'genome_name_to_idx/{graph_name}.json', 'w') as f:
            json.dump(genome_name_to_idx, f)

    def find_collinear_blocks(self):
        collinear_blocks = []

        if SORT_SEEDS=='nr_occurrences':
            vertex_indices_shuffled = sorted(list(range(len(self.vertices))), key=lambda x: len(self.vertices[x].occurrences), reverse=True)
        elif SORT_SEEDS=='length':
            vertex_indices_shuffled = sorted(list(range(len(self.vertices))), key=lambda x: self.vertices[x].length, reverse=True)
        else:
            vertex_indices_shuffled = list(range(len(self.vertices)))
            random.shuffle(vertex_indices_shuffled)
        
        for v_idx in vertex_indices_shuffled: # select a vertex --- seed of a new CollinearBlock
            v = self.vertices[v_idx]
            collinear_seeds = [] # list to store occurrences of v not used before
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
                    collinear_seeds.append(CollinearWalk(g_idx, i, i, orient)) # namedtuple('CollinearWalk', ['genome', 'start', 'end', 'orient'])
                    walk_start_end_check(collinear_seeds[-1], len(genome.path))
            if not collinear_seeds:
                continue
            best_score = -1
            best_block = None
            new_block = CollinearBlock(v_idx, collinear_seeds, carrying_seed_orientation)
            new_score = 0
            while new_score>=0:
                w0_idx = new_block.carrying_path[-1]
                Q = BlockExtensions(new_block.collinear_walks, self, w0_idx) # collinear_walks, graph, w0
                r = Q.get_carrying_path_extension(self, new_block.collinear_walks) # r - shortest walk from w0 to t (vertex, orientation), where t - index of the most frequently visited vertex;
                if not r:
                    break
                if r[-1].vertex==w0_idx:
                    print('\tr[-1].vertex==w0_idx')
                    break
                if r[0].vertex==w0_idx:
                    print('\tr[0].vertex==w0_idx')
                    break
                for wi in r:
                    new_block.carrying_path.append(wi.vertex)
                    new_block.carrying_path_orientations.append(wi.orientation)
                    new_block.update_collinear_walks(wi, Q, self)
                    new_score = new_block.scoring_function(self)
                    new_block.carrying_path_length_so_far += self.vertices[wi.vertex].length
                    if math.isinf(new_score):
                        print('-INFINITY!!!!!!!!!!!!!')
                        break
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
    scores: list
    def __init__(self, seed_idx, seed_occurrences, carrying_seed_orientation):
        self.carrying_path = [seed_idx] # vertex index
        self.carrying_path_orientations = [carrying_seed_orientation] # orientation
        self.collinear_walks = seed_occurrences # one of the collinear walks starts with the same vertex as the carrying path
        self.scores = [Score(0, 0, w.start, 1) for w in self.collinear_walks] # ['q1', 'q3', 'idx_g_path', 'idx_c_path']
        self.carrying_path_length_so_far = 0

    def scoring_function(self, graph):
        score = 0
        for w_idx, walk in enumerate(self.collinear_walks):
            p = walk_length(walk, graph)
            if p>=PARAM_m:
                s = self.scores[w_idx]
                if s.q1>PARAM_b or s.q3>PARAM_b:
                    print('Return -inf!!!!!!!!!!!')
                    return -math.inf
                score += p - (s.q1 + s.q3)**2
        return score

    def update_collinear_walks(self, wi_info, block_extensions, graph): # vertex_info - tuple(vertex index, orientation)
        wi = graph.vertices[wi_info.vertex]
        wi_orient_on_carrying_path = None # calculated only if needed, once
        walks_updated_score = set()
        for occurrence in wi.occurrences:
            # 1) search for walks, whose extensions contain the occurrence of w
            g_idx = occurrence.genome
            g_path = graph.genomes[g_idx].path
            g_len = len(g_path)
            o_nr_on_path = occurrence.nr_on_path
            if g_path[o_nr_on_path].used==True:
                continue
            walk_to_extend = None
            for e_idx, extension in block_extensions.extensions.items():
                if extension.genome==g_idx and extension.start<=o_nr_on_path<=extension.end:
                    if walk_to_extend is None:
                        walk_to_extend = e_idx
                    # 2) if there are multiple walks containing this occurrence of w, 
                    # select the one ending further
                    else:
                        if extension.orient==1:
                            if self.collinear_walks[e_idx].end>self.collinear_walks[walk_to_extend].end:
                                walk_to_extend = e_idx
                        else:
                            if self.collinear_walks[e_idx].start<self.collinear_walks[walk_to_extend].start:
                                walk_to_extend = e_idx
            # 3a) if the path is found, extend it till w AND update the extension!
            if walk_to_extend is not None:
                extension = block_extensions.extensions[walk_to_extend]
                walk = self.collinear_walks[walk_to_extend]
                if walk.orient==1:
                    assert walk.end < o_nr_on_path # walk.start < o_nr_on_path
                    self.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, walk.start, o_nr_on_path, walk.orient)
                else:
                    assert o_nr_on_path < walk.start # o_nr_on_path < walk.end
                    self.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, o_nr_on_path, walk.end, walk.orient)
                walk_start_end_check(self.collinear_walks[walk_to_extend], g_len)
                block_extensions.update_extension(self.collinear_walks[walk_to_extend], graph, collinear_walk_nr=walk_to_extend)
                old_score = self.scores[walk_to_extend]
                if old_score.q3>PARAM_b and walk_length(self.collinear_walks[walk_to_extend], graph)>=PARAM_m:
                    self.scores = [Score(old_score.q3, 0, 0, 0)]
                    print('Too high q3 ---> score<0')
                    return
                self.scores[walk_to_extend] = Score(old_score.q1, 0, o_nr_on_path, len(self.carrying_path))
                walks_updated_score.add(walk_to_extend)
                
            # 3b) if such walk is not found, occurrence becomes a new collinear path (provided it is not used)
            else:
                if wi_orient_on_carrying_path is None:
                    for i, v_idx in enumerate(self.carrying_path):
                        if v_idx==wi_info.vertex:
                            wi_orient_on_carrying_path = self.carrying_path_orientations[i]
                            break
                occurrence_orientation = g_path[o_nr_on_path].orientation
                orient = wi_orient_on_carrying_path * occurrence_orientation
                self.collinear_walks.append(CollinearWalk(g_idx, o_nr_on_path, o_nr_on_path, orient))
                walk_start_end_check(self.collinear_walks[-1], g_len)
                # Create extension for the new collinear walk
                block_extensions.update_extension(self.collinear_walks[-1], graph, collinear_walk_nr=len(self.collinear_walks)-1)
                if wi.length>=PARAM_m and self.carrying_path_length_so_far>PARAM_b:
                    self.scores = [Score(self.carrying_path_length_so_far, 0, 0, 0)]
                    print('Too high q1 ---> score<0')
                    return # We don't check any other occurrences and set block's score to -1
                self.scores.append(Score(self.carrying_path_length_so_far, 0, o_nr_on_path, len(self.carrying_path))) # ['q1', 'q3', 'idx_g_path', 'idx_c_path']
                
        # Update scores
        for w_idx in range(len(self.collinear_walks)):
            if w_idx not in walks_updated_score:
                old_score = self.scores[w_idx]
                if old_score.q3+wi.length>PARAM_b and walk_length(self.collinear_walks[w_idx], graph)>=PARAM_m:
                    self.scores = [Score(old_score.q3+wi.length, 0, 0, 0)]
                    print(f'Too high q1 ---> score<0 (walk introduced too far.)')
                    return
                self.scores[w_idx] = Score(old_score.q1, old_score.q3+wi.length, old_score.idx_g_path, len(self.carrying_path))
                  
class BlockExtensions:
    extensions: dict # list of extensions of type CollinearWalk
    coverage: dict # {vertex index: coverage}, where coverage - number of occurrences in extensions
    shortest_walk: dict  # {vertex index: (distance, *shortest_walk)}
                    # distance - length of the shortest walk from w_0
                    # shortest_walk (walk number, start, end)
    def __init__(self, collinear_walks, graph, w0_idx):
        self.extensions = {}
        self.coverage = {}
        self.shortest_walk = {}
        
        for walk_nr, collinear_walk in enumerate(collinear_walks):
            g_idx = collinear_walk.genome
            g_path = graph.genomes[g_idx].path
            g_len = len(g_path)
            if collinear_walk.orient==1 and collinear_walk.end==g_len-1: # if we've reached the end of the genome
                continue
            elif collinear_walk.orient==-1 and collinear_walk.start==0:
                continue
            else:
                if collinear_walk.orient==1:
                    to_search_from = collinear_walk.end
                    to_end_search = g_len - 1
                else:
                    to_search_from = collinear_walk.start
                    to_end_search = 0
                proximal = to_search_from + collinear_walk.orient # proximal - the end of extension proximal to the collinear walk
                assert 0<=proximal<g_len, f'0<=proximal<genome_length is not true! {proximal=}, {g_len=}'
                assert 0<=to_search_from<g_len, f'0<=to_search_from<genome_length is not true! {to_search_from=}, {g_len=}'
                
                w0_nr_on_path = find_vertex_on_path_till_b(graph, collinear_walk, w0_idx, 
                                                           proximal=to_search_from,
                                                           distal=to_end_search) # taking the closest position to the collinear walk
                # Shouldn't we take all relevant occurrences of w0? <--- TO FIX?
                p_length = 0
                i = proximal
                if w0_nr_on_path is None: # only increase coverage
                    while True:
                        if g_path[i].used==True:
                            if i!=proximal:
                                i -= collinear_walk.orient
                            break
                        v_idx = g_path[i].vertex
                        p_length += graph.vertices[v_idx].length
                        if p_length>PARAM_b: # is it ok? (TO FIX?)
                            if i!=proximal:
                                i -= collinear_walk.orient
                            break
                        self.coverage[v_idx] = self.coverage.get(v_idx, 0) + 1
                        if i in {0, g_len-1}:
                            break
                        i += collinear_walk.orient
                else:
                    while True:
                        if g_path[i].used==True:
                            if i!=proximal:
                                i -= collinear_walk.orient
                            break
                        v_idx = g_path[i].vertex
                        p_length += graph.vertices[v_idx].length
                        if p_length>PARAM_b: # is it ok? (TO FIX?)
                            if i!=proximal:
                                i -= collinear_walk.orient
                            break
                        self.coverage[v_idx] = self.coverage.get(v_idx, 0) + 1
                        is_w0_before_v = (i-w0_nr_on_path)*collinear_walk.orient>=0
                        if is_w0_before_v:
                            distance = walk_length(CollinearWalk(g_idx, w0_nr_on_path, i, 1), graph)
                            if v_idx not in self.shortest_walk or self.shortest_walk[v_idx].distance>distance:
                                self.shortest_walk[v_idx] = PathFromW0(distance, walk_nr, w0_nr_on_path+collinear_walk.orient, i)
                                assert self.shortest_walk[v_idx].w0_nr_on_path<g_len, f'{self.shortest_walk[v_idx].w0_nr_on_path=}, {g_len=}'
                                assert self.shortest_walk[v_idx].t_nr_on_path<g_len, f'{self.shortest_walk[v_idx].t_nr_on_path=}, {g_len=}'
                        if i in {0, g_len-1}:
                            break
                        i += collinear_walk.orient # going forwards or backwards on the genome
            self.extensions[walk_nr] = CollinearWalk(g_idx, min(proximal, i), max(proximal, i), collinear_walk.orient)
            if collinear_walk.orient==1:
                assert self.extensions[walk_nr].start==collinear_walk.end+1, f'{self.extensions[walk_nr]=}, {collinear_walk=}, {g_len=}'
            else:
                assert self.extensions[walk_nr].end==collinear_walk.start-1, f'{self.extensions[walk_nr]=}, {collinear_walk=}, {g_len=}'
            walk_start_end_check(self.extensions[walk_nr], g_len)

    
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
        
        r = []
        walk = collinear_walks[w0_to_t.walk_nr] 
        genome = graph.genomes[walk.genome]
        to_start = w0_to_t.w0_nr_on_path + walk.orient
        to_end = w0_to_t.t_nr_on_path + walk.orient
        for i in range(to_start, to_end):
            g_path_pos = genome.path[i]
            r.append(CarryingPathExtension(g_path_pos.vertex, g_path_pos.orientation*walk.orient))
        return r

    def update_extension(self, walk, graph, collinear_walk_nr):
        g_idx = walk.genome
        genome = graph.genomes[g_idx]
        g_len = len(genome.path)
        proximal = walk.end+1 if walk.orient==1 else walk.start-1
        if proximal>=g_len or proximal<0 or genome.path[proximal].used==True or graph.vertices[genome.path[proximal].vertex].length>PARAM_b:
            if collinear_walk_nr in self.extensions:
                del self.extensions[collinear_walk_nr]
        else:
            distal = proximal
            while True:
                g_path_pos = genome.path[distal]
                extension = CollinearWalk(g_idx, min(proximal, distal), max(proximal, distal), walk.orient)
                if g_path_pos.used==True or walk_length(extension, graph)>PARAM_b:
                    if distal!=proximal:
                        distal -= walk.orient
                    break
                distal += walk.orient
                if distal<0 or distal>g_len-1:
                    distal -= walk.orient
                    break
            self.extensions[collinear_walk_nr] = CollinearWalk(g_idx, min(proximal,distal), max(proximal,distal), walk.orient)
            if walk.orient==1:
                assert self.extensions[collinear_walk_nr].start==walk.end+1, f'{self.extensions[collinear_walk_nr]=}, {walk=}, {g_len=}'
            else:
                assert self.extensions[collinear_walk_nr].end==walk.start-1, f'{self.extensions[collinear_walk_nr]=}, {walk=}, {g_len=}'
            walk_start_end_check(self.extensions[collinear_walk_nr], g_len)

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
    for i in range(proximal, distal+walk.orient):
        if genome_path[i].vertex==v_to_find:
            return i
    return None

def find_vertex_on_path_till_b(graph:Graph, walk:CollinearWalk, v_to_find:int, proximal, distal):
    '''
    Function searches for an occurance of vertex with index v_to_find on walk walk.
    The search begins at index proximal and ends at index distal of walk's genome's path.
    Function returns the first index found or None, if no index has been found.
    '''
    genome_path = graph.genomes[walk.genome].path
    length = 0
    for i in range(proximal, distal+walk.orient):
        v_idx = genome_path[i].vertex
        if genome_path[i].used==True:
            break
        length += graph.vertices[v_idx].length
        if length>PARAM_b: # maybe change it?
            break
        if v_idx==v_to_find:
            return i
    return None

def mark_vertices_as_used(graph, block):
    '''
    Functions marks vertex occurrences of CollinearBlock block as used, 
    i.e. it changes each vertex occurrence's parameter used to True.
    '''
    nr_used = 0
    for collinear_walk in block.collinear_walks:
        g_idx = collinear_walk.genome
        genome_path = graph.genomes[g_idx].path
        for i in range(collinear_walk.start, collinear_walk.end+1):
            g_path_pos = genome_path[i]
            genome_path[i] = Path(*(g_path_pos[:-2]), True, g_path_pos[-1])
        nr_used += walk_length(collinear_walk, graph)
    print(f'Marked as used: {nr_used}.')

def walk_length(walk, graph, start=None, end=None):
        g_idx = walk.genome
        genome_path = graph.genomes[g_idx].path
        if start is None:
            start = walk.start
        if end is None:
            end = walk.end
        assert start <= end
        if start==0:
            return genome_path[end].length
        else:
            return genome_path[end].length - genome_path[start-1].length

def walk_start_end_check(walk, genome_length):
    assert walk.end<genome_length, f'end >= genome_length; {walk=}, {genome_length=}'
    assert walk.start>=0, f'start < 0; {walk=}, {genome_length-1=}'
    assert walk.start<=walk.end, f'end < start; {walk=}'

def save_blocks(blocks:list[CollinearWalk], graph_name, graph):
    df_all = pd.DataFrame()
    walk_lengths = []
    for b_nr, block in enumerate(blocks):
        for walk in block.collinear_walks:
            walk_lengths.append(walk_length(walk, graph))
        df = pd.DataFrame.from_records(block.collinear_walks, columns=CollinearWalk._fields)
        df['block_nr'] = b_nr
        df_all = pd.concat([df_all, df])
        df_all['walk_length'] = walk_lengths
    today = str(date.today()).replace('-', '_')
    
    if SORT_SEEDS=='nr_occurrences':
        name = f'{graph_name}_{today}_sort_by_nr_occurrences.csv'
    elif SORT_SEEDS=='length':
        name = f'{graph_name}_{today}_sort_by_length.csv'
    else:
        name = f'{graph_name}_{today}.csv'
    df_all.to_csv(f'blocks/{name}', index=False)


nr_blocks = []
for graph_file_path in os.listdir('data'):
    for SORT_SEEDS in ['nr_occurrences', 'length', 'no'][:2]:
        print(f'{graph_file_path.upper()}, {SORT_SEEDS=}')
        g = Graph('data/'+graph_file_path)
        blocks = g.find_collinear_blocks()
        nr_blocks.append(len(blocks))
        save_blocks(blocks, graph_file_path.split('.')[0], g)
print(f'Number of blocks for consecutive options:\n{nr_blocks}')

# additional check
SORT_SEEDS = 'no'
for graph_file_path in os.listdir('data'):
    for i in range(10):
        print(f'{graph_file_path.upper()}, {SORT_SEEDS=}')
        g = Graph('data/'+graph_file_path)
        blocks = g.find_collinear_blocks()
        nr_blocks.append(len(blocks))
        save_blocks(blocks, graph_file_path.split('.')[0], g)
print(f'Number of blocks for consecutive options:\n{nr_blocks}')


'''
IMPORTANT NOTES
start, end and orient once again
    - always: start<end
    - if orient==1, we extend end; otherwise - we extend start.

Pomysły (nie jak w artykule): 
    - Update'ować Q dla danej ścieżki po małym kroku, jeśli znalazłam wspólny wierzchołek z carrying path. <--- DONE
    - Iść po ścieżce i szukać czegokolwiek, co jest w r. Wtedy sprawdzać, na czym najlepiej się zatrzymać.

TO FIX:
    - Po zakończeniu dużego kroku: zapamiętać rozszerzenia i z nich korzystać.
    - Po dużym kroku zapamiętać scory poszczególnych ścieżek do końca łańcucha i później doliczać tylko odtąd.
    - Budując shortest_walk w przedłużeniach: chyyyba trzeba wziąć pod uwagę 
    wszystkie wystąpienia w0 na przedłużeniu. Wtedy nie przegapimy najkrótszej drogi.
    - Dlaczego są zapętlenia? Może trzeba sprawdzać rzeczywiście przedłużenia szukając r?

Not necessary:
    Pamiętać po dużym kroku, które wierzchołki się końćzą w starym t (nowym w0).
'''