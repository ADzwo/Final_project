import math
from config import *
from tuples import *
from datetime import date
import pandas as pd
import numpy as np

def walk_length(walk, graph, start=None, end=None):
    '''
    Function returns length of a walk in a graph (in nucleotides).
    Arguments:
     - walk - CollinearWalk, which we measure,
     - graph - Graph, in which we measure walk,
     - start, end - integers - the postions in walk's genome path, 
     between which we measure, both default to None.
    If start / end is None, start / end of the walk is used instead (accordingly).
    Assumption: start<=end (also after replacing None with walk.start / walk.end).
    '''
    if start is None:
        start = walk.start
    if end is None:
        end = walk.end
    if start > end:
        raise ValueError(f'Function walk_length got start > end. {start=}, {end=}')
    if start<0:
        raise ValueError(f'Got negative start position, {start=}.')
    
    g_idx = walk.genome
    genome_path = graph.genomes[g_idx].path
    if end>=len(genome_path):
        raise ValueError(f'Got end value > genome path length. {end=}, {len(genome_path)=}')
    
    if start==0:
        return genome_path[end].p_length
    else:
        return genome_path[end].p_length - genome_path[start-1].p_length
    

def find_vertex_on_path_till_b(graph, walk, v_to_find, proximal, distal):
    '''
    Function searches for an occurance of vertex with index v_to_find on walk walk.
    The search begins at index proximal and ends at index distal of walk's genome's path.
    Function returns the first index found or None, if no index has been found.
    '''
    genome_path = graph.genomes[walk.genome].path
    length = 0
    assert (distal-proximal)*walk.orient>=0
    for i in range(proximal, distal+walk.orient, walk.orient):
        v_idx = genome_path[i].vertex
        if genome_path[i].used==True:
            break
        if v_idx==v_to_find:
            return i
        length += graph.vertices[v_idx].v_length
        if length>=PARAM_b:
            break
    return None

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

def remove_short_walks(block, graph):
    '''
    Function removes walks shorter than PARAM_m from collinear blocks.
    '''
    for w_idx in range(len(block.collinear_walks)-1, -1, -1):
        if walk_length(block.collinear_walks[w_idx], graph)<PARAM_m:
            block.collinear_walks.pop(w_idx)
            block.match_carrying.pop(w_idx)
            # block.scores.pop(w_idx)
    for walk in block.collinear_walks:
        assert walk_length(walk, graph)>=PARAM_m

def scoring_function(block, graph):
    '''
    Function returns score for a given block in a graph.
    Block's score is a sum of scores of collinear walks.
    Walk's score is p - (q1 + q3)**2, where
     - p is walk length,
     - q1 and q3 are lengths of carrying path's hanging ends
     (fragments which do not form a chain with the walk).
    '''
    score = 0
    for w_idx, walk in enumerate(block.collinear_walks):
        p = walk_length(walk, graph)
        if p>=PARAM_m:
            s = block.scores[w_idx]
            if s.q1>PARAM_b or s.q3>PARAM_b:
                return -math.inf
            score += p - (s.q1 + s.q3)**2
    return score

def save_blocks_to_gff(blocks:list, graph, graph_name=None):
    '''
    Function saves collinear blocks in a .gff file.
    Consecutive columns represent following information about collinear walks.
    - seqname: genome id
    - source: software name
    - feature: '.'
    - start: start position in the genome
    - end: end position in the genome
    - score: '.'
    - strand: orientation of the walk relative to the genome
    - frame: '.'
    - attribute: block_nr (in form 'ID=block_nr')
    Name of the .gff file consists of graph.name, today's date and the sorting seed mode.
    
    Example line
    3	final_project	.	212	283	.	+	.	ID=1
    '''
    gff_cols = ['seqname', 'source', 'feature', 'start', 
                'end', 'score', 'strand', 'frame', 'attribute']
    df_all = pd.DataFrame()
    walk_starts = []
    walk_ends = []
    for b_nr, block in enumerate(blocks):
        for walk in block.collinear_walks:
            walk_starts.append(walk_length(walk, graph, end=walk.start)-1)
            walk_ends.append(walk_length(walk, graph)-1)
            assert walk_starts[-1]>=0 and walk_ends[-1]>=0, f'{walk_starts[-1]=}, {walk_ends[-1]=}'
        df = pd.DataFrame.from_records(block.collinear_walks, columns=CollinearWalk._fields)
        df['attribute'] = f'ID={b_nr}'
        df_all = pd.concat([df_all, df])
    df_all['start'] = walk_starts
    df_all['end'] = walk_ends
    df_all['source'] = 'final_project'
    df_all['feature'] = '.'
    df_all['orient'] = np.where(df_all['orient']>0, '+', '-')
    df_all['score'] = '.'
    df_all['frame'] = '.'
    df_all.rename(columns={'genome':'seqname', 'orient':'strand'}, inplace=True)
    df_all = df_all[gff_cols]
    today = str(date.today()).replace('-', '_')

    if graph_name is None:
        graph_name = graph.name
    if SORT_SEEDS=='nr_occurrences':
        name = f'{graph_name}_{today}_sort_by_nr_occurrences.gff'
    elif SORT_SEEDS=='length':
        name = f'{graph_name}_{today}_sort_by_length.gff'
    else:
        name = f'{graph_name}_{today}.gff'
    print(f'{name=}')
    df_all.to_csv(f'blocks/{name}', index=False)