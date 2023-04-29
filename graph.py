import os
import sys
import re
import math
import json
import random
from collections import namedtuple, defaultdict
import pandas as pd
import numpy as np
from datetime import date

PARAM_m = 50
PARAM_b = 200
PARAM_a = 30 # abundance pruning parameter; 150 used for mice dataset in SibeliaZ paper

complement_dict = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}

SORT_SEEDS = ['no', 'nr_occurrences', 'length'][2]
assert SORT_SEEDS in {'no', 'nr_occurrences', 'length'}, f'SORT_SEEDS must be one of "no", "nr_occurrences", "length", got {SORT_SEEDS}'

CollinearWalk = namedtuple('CollinearWalk', ['genome', 'start', 'end', 'orient'])
Path = namedtuple('Path', ['vertex', 'orientation', 'used', 'p_length']) # length --- length of the suffix of the genome, ending at the last character of this vertex
PathFromW0 = namedtuple('ExtensionVertex', ['distance', 'walk_nr', 'w0_nr_on_path', 't_nr_on_path'])
CarryingPathExtension = namedtuple('CarryingPathExtension', ['vertex', 'orientation'])
occurrence = namedtuple('occurrence', ['genome', 'nr_on_path'])
Score = namedtuple('Score', ['q1', 'q3'])

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

    def __init__(self, graph_file_path):
        self.vertices = []
        self.genomes = []
        vertex_name_to_idx = {} # dict to save {vertex id from gfa file: index in Graph.vertices}
        genome_name_to_idx = {}
        sequences = []
        graph_name = re.split(r'[\.\/]', graph_file_path)[-2]

        # find S lines and fill vertex_name_to_idx dict 
        with open(graph_file_path,'r') as graph_file:
            for line in graph_file.readlines():
                if line.startswith('S'):
                    v = line.strip().split()
                    vertex_name_to_idx[v[1]] = len(self.vertices) # {name: id}
                    self.vertices.append(Vertex(len(v[2])))
                    sequences.append(v[2]+'\n')
        with open(f'vertex_sequences/{graph_name}.txt', 'w') as f:
            f.writelines(sequences)
            del sequences
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
                        length += self.vertices[v_idx].v_length
                        path.append(Path(v_idx, v_orientation, False, length))
                    genome_name_to_idx[g[1]] = len(self.genomes)
                    self.genomes.append(Genome(path))
        # save dictionaries to .json files
        with open(f'vertex_name_to_idx/{graph_name}.json', 'w') as f:
            json.dump(vertex_name_to_idx, f)
        with open(f'genome_name_to_idx/{graph_name}.json', 'w') as f:
            json.dump(genome_name_to_idx, f)

    def find_collinear_blocks(self):
        collinear_blocks = []

        if SORT_SEEDS=='nr_occurrences':
            vertex_indices_order = sorted(list(range(len(self.vertices))), key=lambda x: len(self.vertices[x].occurrences), reverse=True)
        elif SORT_SEEDS=='length':
            vertex_indices_order = sorted(list(range(len(self.vertices))), key=lambda x: self.vertices[x].v_length, reverse=True)
        else:
            vertex_indices_order = list(range(len(self.vertices)))
            random.shuffle(vertex_indices_order)
        
        for v_idx in vertex_indices_order: # select a vertex --- seed of a new CollinearBlock
            v = self.vertices[v_idx]
            # Check if junction with occ len > a
            if len(v.occurrences) > PARAM_a:
                continue
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
                    collinear_seeds.append(CollinearWalk(g_idx, i, i, orient))
                    walk_start_end_check(collinear_seeds[-1], len(genome.path))
            if not collinear_seeds:
                continue
            
            blocks_forward_backward = {}
            for forward_backward, seeds in zip([1, -1], [collinear_seeds, [CollinearWalk(*s[:-1], -s.orient) for s in collinear_seeds]]):
                best_score = -1
                best_block = None
                new_block = CollinearBlock(v_idx, seeds.copy(), carrying_seed_orientation)
                for walk in new_block.collinear_walks:
                    assert self.genomes[walk.genome].path[walk.start].vertex in new_block.carrying_path, f'{self.genomes[walk.genome].path[walk.start].vertex=},\n{new_block.carrying_path=}'
                    assert self.genomes[walk.genome].path[walk.end].vertex in new_block.carrying_path, f'{self.genomes[walk.genome].path[walk.start].vertex=},\n{self.genomes[walk.genome].path[walk.end].vertex=},\n{new_block.carrying_path=}'
                new_score = 0
                while new_score>=0:
                    assert len(new_block.match_carrying) == len(new_block.collinear_walks)
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
                        new_block.carrying_path_length_so_far += self.vertices[wi.vertex].v_length
                        if math.isinf(new_score):
                            print('-INF!!!!!!!!!')
                            break
                        if new_score>best_score:
                            best_block = new_block
                            best_score = new_score 
                    for s_nr, seed in enumerate(collinear_seeds):
                        assert seed.genome==new_block.collinear_walks[s_nr].genome
                        
                if best_score>0:
                    blocks_forward_backward[forward_backward] = best_block

            if blocks_forward_backward:
                final_block = merge_forward_backward_blocks(blocks_forward_backward, len(collinear_seeds))
                mark_vertices_as_used(self, final_block)
                collinear_blocks.append(final_block)
                
        return collinear_blocks
                

class CollinearBlock:
    carrying_path: list # vertex index: int
    carrying_path_orientations: list # orientation: int
    collinear_walks: list # genome index: int,
                                    # start index from genome's path: int,
                                    # end index from genome's path: int,
                                    # orientaion with respect to the path: int
    scores: list
    carrying_path_length_so_far: int
    match_carrying: list # a list of the same length as block.collinear_walks, 
                        # containing lists of tuples 
                        # (position on carrying path - from 0 to len(carrying_path)-1, 
                        #  position on walk.genome - between walk.start and walk.end)
    def __init__(self, seed_idx, seed_occurrences, carrying_seed_orientation):
        self.carrying_path = [seed_idx] # vertex index
        self.carrying_path_orientations = [carrying_seed_orientation] # orientation
        self.collinear_walks = seed_occurrences # one of the collinear walks starts with the same vertex as the carrying path
        self.scores = [Score(0, 0) for w in self.collinear_walks] # ['q1', 'q3']
        self.carrying_path_length_so_far = 0
        self.match_carrying = [[(0, w.start)] for w in self.collinear_walks]

    def scoring_function(self, graph):
        score = 0
        for w_idx, walk in enumerate(self.collinear_walks):
            p = walk_length(walk, graph)
            if p>=PARAM_m:
                s = self.scores[w_idx]
                if s.q1>PARAM_b or s.q3>PARAM_b:
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
                assert g_path[walk.end].vertex in self.carrying_path
                assert g_path[walk.start].vertex in self.carrying_path
                assert g_path[o_nr_on_path].vertex in self.carrying_path
                if walk.orient==1:
                    self.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, walk.start, o_nr_on_path, walk.orient)
                else:
                    self.collinear_walks[walk_to_extend] = CollinearWalk(walk.genome, o_nr_on_path, walk.end, walk.orient)
                walk_start_end_check(self.collinear_walks[walk_to_extend], g_len)
                block_extensions.update_extension(self.collinear_walks[walk_to_extend], graph, walk_idx=walk_to_extend)
                old_score = self.scores[walk_to_extend]
                self.scores[walk_to_extend] = Score(old_score.q1, 0)
                walks_updated_score.add(walk_to_extend)
                self.match_carrying[walk_to_extend].append((len(self.carrying_path)-1, o_nr_on_path))
                
            # 3b) if such a walk is not found, occurrence becomes a new collinear path (provided it is not used)
            else:
                if wi_orient_on_carrying_path is None:
                    for i, v_idx in enumerate(self.carrying_path):
                        if v_idx==wi_info.vertex:
                            wi_orient_on_carrying_path = self.carrying_path_orientations[i]
                            break
                occurrence_orientation = g_path[o_nr_on_path].orientation
                orient = wi_orient_on_carrying_path * occurrence_orientation
                assert g_path[o_nr_on_path].vertex in self.carrying_path
                self.collinear_walks.append(CollinearWalk(g_idx, o_nr_on_path, o_nr_on_path, orient))
                walk_start_end_check(self.collinear_walks[-1], g_len)
                # Create extension for the new collinear walk
                block_extensions.update_extension(self.collinear_walks[-1], graph, walk_idx=len(self.collinear_walks)-1)
                self.scores.append(Score(self.carrying_path_length_so_far, 0))
                walks_updated_score.add(len(self.collinear_walks)-1)
                self.match_carrying.append([(len(self.carrying_path)-1, o_nr_on_path)])
                
        # Update scores
        for w_idx in range(len(self.collinear_walks)):
            if w_idx not in walks_updated_score:
                old_score = self.scores[w_idx]
                self.scores[w_idx] = Score(old_score.q1, old_score.q3+wi.v_length)
                  
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
                    to_search_from = collinear_walk.end # we start searching at the most distal position of the walk
                    to_end_search = g_len - 1
                else:
                    to_search_from = collinear_walk.start
                    to_end_search = 0
                proximal = to_search_from + collinear_walk.orient # proximal - the end of extension proximal to the collinear walk
                assert 0<=proximal<g_len, f'0<=proximal<genome_length is not true! {proximal=}, {g_len=}'
                assert 0<=to_search_from<g_len, f'0<=to_search_from<genome_length is not true! {to_search_from=}, {g_len=}'
                
                if g_path[proximal].used==True:
                    break
                w0_nr_on_path = find_vertex_on_path_till_b(graph, collinear_walk, w0_idx, 
                                                           proximal=to_search_from,
                                                           distal=to_end_search) # taking the closest position to the collinear walk
                                # Shouldn't we take all relevant occurrences of w0? <--- TO FIX?
                p_length = 0
                i = proximal
                if w0_nr_on_path is None: # only increase coverage
                    while True:
                        if g_path[i].used==True:
                            i -= collinear_walk.orient
                            break
                        v_idx = g_path[i].vertex
                        self.coverage[v_idx] = self.coverage.get(v_idx, 0) + 1
                        p_length += graph.vertices[v_idx].v_length
                        if p_length>=PARAM_b:
                            break
                        if i in {0, g_len-1}:
                            break
                        i += collinear_walk.orient
                else:
                    while True:
                        if g_path[i].used==True:
                            i -= collinear_walk.orient
                            break
                        v_idx = g_path[i].vertex
                        self.coverage[v_idx] = self.coverage.get(v_idx, 0) + 1
                        p_length += graph.vertices[v_idx].v_length
                        if p_length>=PARAM_b:
                            break
                        is_w0_before_v = ((i-w0_nr_on_path)*collinear_walk.orient) > 0
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
        for i in range(w0_to_t.w0_nr_on_path+walk.orient, w0_to_t.t_nr_on_path+walk.orient):
            g_path_pos = genome.path[i]
            r.append(CarryingPathExtension(g_path_pos.vertex, g_path_pos.orientation*walk.orient))
        return r

    def update_extension(self, walk, graph, walk_idx):
        g_idx = walk.genome
        genome = graph.genomes[g_idx]
        g_len = len(genome.path)
        proximal = walk.end+1 if walk.orient==1 else walk.start-1
        if proximal>=g_len or proximal<0 or genome.path[proximal].used==True:
            if walk_idx in self.extensions:
                del self.extensions[walk_idx]
        else:
            distal = proximal
            ext_length = 0
            while True:
                g_path_pos = genome.path[distal]
                if g_path_pos.used==True:
                    distal -= walk.orient
                    break
                ext_length += graph.vertices[g_path_pos.vertex].v_length
                if ext_length>=PARAM_b:
                    break
                distal += walk.orient
                if distal<0 or distal>g_len-1:
                    distal -= walk.orient
                    break
            self.extensions[walk_idx] = CollinearWalk(g_idx, min(proximal,distal), max(proximal,distal), walk.orient)
            if walk.orient==1:
                assert self.extensions[walk_idx].start==walk.end+1, f'{self.extensions[walk_idx]=}, {walk=}, {g_len=}'
            else:
                assert self.extensions[walk_idx].end==walk.start-1, f'{self.extensions[walk_idx]=}, {walk=}, {g_len=}'
            walk_start_end_check(self.extensions[walk_idx], g_len)

def find_vertex_on_path_till_b(graph:Graph, walk:CollinearWalk, v_to_find:int, proximal, distal):
    '''
    Function searches for an occurance of vertex with index v_to_find on walk walk.
    The search begins at index proximal and ends at index distal of walk's genome's path.
    Function returns the first index found or None, if no index has been found.
    '''
    genome_path = graph.genomes[walk.genome].path
    length = 0
    # nr_on_path = None
    assert (distal-proximal)*walk.orient>=0
    for i in range(proximal, distal+walk.orient):
        v_idx = genome_path[i].vertex
        if genome_path[i].used==True:
            break
        if v_idx==v_to_find:
            return i
            # nr_on_path = i
        length += graph.vertices[v_idx].v_length
        if length>=PARAM_b: # maybe change it?
            break
    return None # nr_on_path

def mark_vertices_as_used(graph, block):
    '''
    Functions marks vertex occurrences of CollinearBlock block as used, 
    i.e. it changes each vertex occurrence's parameter used to True.
    '''
    for walk in block.collinear_walks:
        g_idx = walk.genome
        genome_path = graph.genomes[g_idx].path
        for i in range(walk.start, walk.end+1):
            g_path_pos = genome_path[i]
            genome_path[i] = Path(*(g_path_pos[:-2]), True, g_path_pos[-1])

def walk_length(walk, graph, start=None, end=None):
        if start is None:
            start = walk.start
        if end is None:
            end = walk.end
        assert start <= end, f'Function walk_length got start > end. {start=}, {end=}'
        g_idx = walk.genome
        genome_path = graph.genomes[g_idx].path
        if start==0:
            return genome_path[end].p_length
        else:
            return genome_path[end].p_length - genome_path[start-1].p_length

def walk_start_end_check(walk, genome_length):
    assert walk.end<genome_length, f'end >= genome_length; {walk=}, {genome_length=}'
    assert walk.start>=0, f'start < 0; {walk=}, {genome_length-1=}'
    assert walk.start<=walk.end, f'end < start; {walk=}'

def merge_forward_backward_blocks(block_dict, nr_seeds):
    assert len(block_dict)<3, 'block list should contain at most two blocks (forward and backward)'
    if len(block_dict)==1:
        return list(block_dict.values())[0]
    if len(block_dict)==2:
        forward_block = block_dict[1]
        backward_block = block_dict[-1]
        for w_idx in range(nr_seeds):
            walk_f = forward_block.collinear_walks[w_idx]
            walk_b = backward_block.collinear_walks[w_idx]
            assert walk_f.genome==walk_b.genome, f'genome should be the same for the seed walks of forward and backward block,\n {walk_f=},\n {walk_b=}'
            forward_block.collinear_walks[w_idx] = CollinearWalk(walk_f.genome, min(walk_f.start,walk_b.start), 
                                                                 max(walk_f.end,walk_b.end), walk_f.orient)
        for walk in backward_block.collinear_walks[nr_seeds:]:
            forward_block.collinear_walks.append(CollinearWalk(*walk[:-1], -walk.orient))
        forward_block.carrying_path = list(reversed(backward_block.carrying_path)) + forward_block.carrying_path
        forward_block.carrying_path_orientations = list(reversed(backward_block.carrying_path_orientations)) + forward_block.carrying_path_orientations
    return forward_block


def reverse_complement(seq):
    seq_complement = ''
    for s in seq[::-1]:
        seq_complement += complement_dict[s]
    return seq_complement

class PoaVertex:
    sequence: str # one of 'A', 'C', 'G', 'T'
    v_idx: int
    v_pos: int
    next_poa: set # set of poa vertex indices
    walks: dict # {w_idx: index on walk.genome}

    def __init__(self, sequence, v_idx, v_pos):
        assert len(sequence)==1, f'Each vertex of POA graph should contain exactly ine character. Got {sequence}.'
        self.sequence = sequence
        self.v_idx = v_idx
        self.v_pos = v_pos
        self.next_poa = set()
        self.walks = {}
        self.score_column = list

class PoaGraph:
    poa_vertices: list # list of poa vertices
    carrying_to_poa: defaultdict # {v_idx_nr_on_carrying_path: [poa_v_idx_0, poa_v_idx_1]}
    score_columns: dict
    score_sources: dict
    align_end: list
    match: int
    gap: int
    mismatch: int

    def __init__(self, block, graph, sequences, match, gap, mismatch):
        self.poa_vertices = [] 
        self.carrying_to_poa = defaultdict(list)
        self.score_columns = {}
        self.score_sources = {}
        self.align_end = []

        self.gap = gap
        self.match = match
        self.mismatch = mismatch
        
        # Add carrying_path
        poa_v_old = None
        for v_idx, v_orientation in zip(block.carrying_path, block.carrying_path_orientations):
            if v_orientation==1:
                seq = sequences[v_idx].strip()
            else:
                seq = reverse_complement(sequences[v_idx].strip())
            for i, s in enumerate(seq):            
                poa_v = PoaVertex(s, v_idx, i)
                if poa_v_old is not None:
                    poa_v_old.next_poa.add(len(self.poa_vertices))
                self.carrying_to_poa[i].append(len(self.poa_vertices))
                self.poa_vertices.append(poa_v)
                poa_v_old = poa_v

        for w_idx in range(len(block.collinear_walks)):
            self.add_collinear_walk(block, graph, w_idx)

    def add_collinear_walk(self, block, graph, w_idx):        
        # Add a collinear walk
        walk = block.collinear_walks[w_idx]
        matches = block.match_carrying[w_idx]
        
        to_start = walk.start+1 if walk.orient==1 else walk.end-1 # end and beginning of each walk are already in the graph
        to_end = walk.end+1 if walk.orient==1 else walk.start-1
        g_path = graph.genomes[walk.genome].path

        poa_idx_old = self.carrying_to_poa[0][-1]
        i_on_g_path = to_start
        i_on_match = 1
        w_len_aligned = 0
        seq = ''
        while i_on_g_path!=to_end: # to_end does not belong to walk
            v_idx = g_path[i_on_g_path].vertex
            if matches[i_on_match][1]!=i_on_g_path: # there is no match on this position of walk
                if walk.orient*g_path[i_on_g_path].orientation==1:
                    seq += sequences[v_idx].strip()
                else:
                    seq += reverse_complement(sequences[v_idx].strip())
            elif seq=='': # there is a match and we do not need to align anything
                self.poa_vertices[poa_idx_old].add(self.carrying_to_poa[matches[i_on_match][0]][0])
                for poa_idx in self.carrying_to_poa[matches[i_on_match][0]]: # poa_idx representing this position on carrying path
                    poa_v = self.poa_vertices[poa_idx]
                    poa_v.walks[w_idx] = w_len_aligned
                    w_len_aligned += 1
                poa_idx_old = poa_idx
                i_on_match += 1
            else:
                match_from = matches[i_on_match-1]
                match_to = matches[i_on_match]
                collinear_from = match_from[1]
                collinear_to = match_to[1]
                poa_idx_from = self.carrying_to_poa[match_from[0]][-1] # last match
                poa_idx_to = self.carrying_to_poa[match_to[0]][0] # next match
                w_len_aligned += len(seq)
                print(f'{poa_idx_from=}, {poa_idx_to=}')
                poa_idx_old = self.align_to_graph(seq, w_idx, poa_idx_from, poa_idx_to, w_len_aligned, collinear_from, collinear_to, graph)
                seq = ''
                # w_len_aligned += len(self.carrying_to_poa[match_to[0]])
                # poa_idx_old = poa_idx_to
                # i_on_match += 1
            i_on_g_path += walk.orient

    def align_to_graph(self, seq, w_idx, poa_idx_from, poa_idx_to, w_len_aligned, collinear_from, collinear_to, graph):
        w_align_to = set(self.poa_vertices[poa_idx_from].walks.keys()).intersection(set(self.poa_vertices[poa_idx_to].walks.keys()))
        
        self.score_columns[poa_idx_from] = np.arange(len(seq)+1) * self.gap
        self.forward_steps(poa_idx_from, poa_idx_to, w_align_to, seq)
        self.backtracking(w_idx, seq, poa_idx_to, w_len_aligned, collinear_from, collinear_to, graph)

        self.align_end = []
        self.score_columns = {}
        self.score_sources = {}
    
    def update_alignment_score(self, poa_idx, poa_idx_next, seq):
        scores = self.score_columns[poa_idx]
        score_column = np.zeros(len(seq)+1)
        score_column[0] = scores[0] + self.gap
        score_sources = np.zeros((len(seq)+1, 3))
        score_sources[0,:] = (-1, 0, poa_idx)
        poa_next = self.poa_vertices[poa_idx_next]
        for i in range(len(seq)):
            is_match = seq[i]==poa_next.sequence
            new_score = scores[i]+self.match if is_match else scores[i]+self.mismatch
            score_source = (-1, -1, poa_idx)
            if scores[i+1]+self.gap>new_score:
                new_score = scores[i+1]+self.gap
                score_source = (-1, 0, poa_idx)
            if score_column[-1]+self.gap>new_score:
                new_score = score_column[-1]+self.gap
                score_source = (0, -1, poa_idx_next)
            score_column[i+1] = new_score
            score_sources[i+1] = score_source

        if poa_idx_next in self.score_columns:
            self.score_columns[poa_idx_next] = np.maximum(self.score_columns[poa_idx_next], score_column)
            self.score_sources[poa_idx_next] = np.where(self.score_columns[poa_idx_next]>score_column, 
                                                        self.score_sources[poa_idx_next], score_sources)
        else:
            self.score_columns[poa_idx_next] = score_column
            self.score_sources[poa_idx_next] = score_sources
        
        assert len(self.score_columns[poa_idx_next]) == len(self.score_columns[poa_idx])
          

    def forward_steps(self, poa_idx_from, poa_idx_to, w_align_to, seq):
        poa_idx = poa_idx_from
        for poa_idx_next in self.poa_vertices[poa_idx].next_poa:
            poa_next = self.poa_vertices[poa_idx_next]
            if poa_idx_next==poa_idx_to:
                self.align_end.append(poa_idx)
                continue
            if len(poa_next.walks)==0 or len(set(poa_next.walks.keys()).intersection(w_align_to))>0:
                self.update_alignment_score(poa_idx, poa_idx_next, seq)
                self.forward_steps(poa_idx_next, poa_idx_to, w_align_to, seq)

    def backtracking(self, w_idx, seq, poa_idx_to, w_len_aligned, collinear_from, collinear_to, graph):
        # w_len_aligned is already after adding len(seq)

        g_path = graph.genomes[block.collinear_walks[w_idx].genome].path
        w_orient = block.collinear_walks[w_idx].orient
        
        row = -1
        first_iteration = True
        for collinear_pos in range(collinear_to+w_orient, collinear_from, -w_orient):
            v_idx = g_path[collinear_pos].vertex
            v_len = graph.vertices[v_idx].v_length
            for poa_pos_on_v in range(v_len-1, -1, -1):
                if first_iteration:
                    poa_idx = poa_idx_to
                    # self.poa_vertices[poa_idx].walks[w_idx] = w_len_aligned
                    print(f'{self.align_end=}')
                    pred_idx = self.align_end[np.argmax([self.score_columns[poa_idx][-1] for poa_idx in self.align_end])]
                    last_poa_idx_aligned = pred_idx
                    first_iteration = False
                else:
                    poa_idx = pred_idx
                    pred_idx = self.score_sources[poa_idx][row][-1]

                source = self.score_sources[poa_idx][row][:-1]
                if source==(0, -1): # insertion in the new walk
                    self.poa_vertices.append(PoaVertex(seq[row], v_idx, poa_pos_on_v))
                    predecessor = self.poa_vertices[-1]
                    predecessor.next_poa.add(poa_idx)
                    predecessor.walks[w_idx] = w_len_aligned + row
                elif source==(-1, 0): # gap
                    predecessor = self.poa_vertices[pred_idx]
                else: # match or mismatch
                    if seq[row]==self.poa_vertices[poa_idx].sequence: # match
                        predecessor = self.poa_vertices[pred_idx]
                    else: # mismatch
                        self.poa_vertices.append(PoaVertex(seq[row], v_idx, poa_pos_on_v))
                        predecessor = self.poa_vertices[-1]
                    predecessor.next_poa.add(poa_idx)
                    predecessor.walks[w_idx] = w_len_aligned + row
                
                row += source[1]
        print(f'{first_iteration=}')
        print(f'{collinear_to+w_orient=}, {collinear_from=}, {w_orient}')
        print(f'{collinear_pos=}')
        print(f'{v_len-1=}')
        return last_poa_idx_aligned

                    
                
def save_blocks(blocks:list, graph_name, graph):
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

os.chdir(sys.path[0])
assert 'data' in os.listdir()
for folder in ['blocks', 'vertex_name_to_idx', 'genome_name_to_idx', 'vertex_sequences']:
    if not os.path.exists(folder):
        os.mkdir(folder)

nr_blocks = []
for graph_file_path in os.listdir('data'):
    for SORT_SEEDS in ['nr_occurrences', 'length', 'no'][:2]:
        print(f'{graph_file_path.upper()}, {SORT_SEEDS=}')
        graph = Graph('data/'+graph_file_path)
        blocks = graph.find_collinear_blocks()
        nr_blocks.append(len(blocks))
        save_blocks(blocks, graph_file_path.split('.')[0], graph)
        graph_name = re.split(r'[\.\/]', graph_file_path)[-2]
        with open(f'vertex_sequences/{graph_name}.txt') as f:
            sequences = f.readlines()

        for block in blocks:
            for walk in block.collinear_walks:
                assert graph.genomes[walk.genome].path[walk.start].vertex in block.carrying_path, f'{graph.genomes[walk.genome].path[walk.start].vertex=}'
                assert graph.genomes[walk.genome].path[walk.end].vertex in block.carrying_path, f'{graph.genomes[walk.genome].path[walk.end].vertex=}'
            PoaGraph(block, graph, sequences, match=5, gap=-1, mismatch=-2)
print(f'Number of blocks for consecutive options:\n{nr_blocks}')

# additional check
SORT_SEEDS = 'no'
for graph_file_path in os.listdir('data'):
    for i in range(10):
        print(f'{graph_file_path.upper()}, {SORT_SEEDS=}')
        graph = Graph('data/'+graph_file_path)
        blocks = graph.find_collinear_blocks()
        nr_blocks.append(len(blocks))
        save_blocks(blocks, graph_file_path.split('.')[0], graph)
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
    - Range z odpowiednim krokiem!!!!
    - Po zakończeniu dużego kroku: zapamiętać rozszerzenia i z nich korzystać.
    - Po dużym kroku zapamiętać scory poszczególnych ścieżek do końca łańcucha i później doliczać tylko odtąd.
    - Budując shortest_walk w przedłużeniach: chyyyba trzeba wziąć pod uwagę 
    wszystkie wystąpienia w0 na przedłużeniu. Wtedy nie przegapimy najkrótszej drogi.

Not necessary:
    Pamiętać po dużym kroku, które wierzchołki się końćzą w starym t (nowym w0).
'''