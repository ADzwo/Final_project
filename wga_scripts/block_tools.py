import math, random
from tuples import *
from datetime import date
import pandas as pd
import numpy as np

def sort_v_idx(graph, SORT_SEEDS):
    vertex_indices_order = list(range(len(graph.vertices)))
    if SORT_SEEDS=='nr_occurrences':
        vertex_indices_order = sorted(vertex_indices_order, key=lambda x: len(graph.vertices[x].occurrences), reverse=True)
    elif SORT_SEEDS=='length':
        vertex_indices_order = sorted(vertex_indices_order, key=lambda x: graph.vertices[x].v_length, reverse=True)
    else:
        random.shuffle(vertex_indices_order)
    return vertex_indices_order

    
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
    

def find_vertex_on_path_till_b(graph, walk, v_to_find, v_to_find_orentation, 
                               proximal, distal, PARAM_b):
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
            v_orientation = genome_path[i].orientation*walk.orient
            if v_orientation==v_to_find_orentation:
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
            graph.genomes[g_idx].path[i] = Path(*(g_path_pos[:-2]), True, g_path_pos[-1])

def remove_short_walks(block, graph, PARAM_m):
    '''
    Function removes walks shorter than PARAM_m from collinear blocks.
    '''
    for w_idx in range(len(block.collinear_walks)-1, -1, -1):
        if walk_length(block.collinear_walks[w_idx], graph)<PARAM_m:
            block.collinear_walks.pop(w_idx)
            block.match_carrying.pop(w_idx)

def scoring_function(block, graph, PARAM_m, PARAM_b):
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

def get_file_name(graph_name, SORT_SEEDS, extension):
    today = str(date.today()).replace('-', '_')
    if SORT_SEEDS=='nr_occurrences':
        return f'{graph_name}_{today}_sort_by_nr_occurrences.{extension}'
    elif SORT_SEEDS=='length':
        return f'{graph_name}_{today}_sort_by_length.{extension}'
    else:
        return f'{graph_name}_{today}.{extension}'

def save_block_to_gff(block, graph, block_nr, file):
    '''
    Function saves a collinear block in a .gff file.
    The block is appended to an existing .csv file
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
    walk_starts = []
    walk_ends = []
    for walk in block.collinear_walks:
        g_idx = walk.genome
        genome_path = graph.genomes[g_idx].path
        if walk.start==0:
            walk_starts.append(0)
        else:
            walk_starts.append(genome_path[walk.start-1].p_length)
        walk_ends.append(genome_path[walk.end].p_length-1)
        assert walk_starts[-1]>=0 and walk_ends[-1]>=0, f'{walk_starts[-1]=}, {walk_ends[-1]=}'
    df = pd.DataFrame.from_records(block.collinear_walks, columns=CollinearWalk._fields)
    if len(df)!=len(df.drop_duplicates()):
        print(f'Dataframe of block {block_nr} contains duplicates!')
    df['attribute'] = f'ID={block_nr}'
    df['start'] = walk_starts
    df['end'] = walk_ends
    df['orient'] = np.where(df['orient']>0, '+', '-')
    df.rename(columns={'genome':'seqname', 'orient':'strand'}, inplace=True)
    
    df['source'] = 'final_project'
    df['feature'] = '.'
    df['score'] = '.'
    df['frame'] = '.'
    df = df[gff_cols]

    df.to_csv(file, index=False, mode='w')
    return df[['start', 'end', 'strand']]

def save_maf(alignment, maf_file, block_df, var_graph, walks, genome_idx_to_name):
    maf_file.write('a\n')
    # The first line of the alignment represents the carrying path (we add it later).
    if len(block_df)+1!=len(alignment):
        raise ValueError(f'Block size + 1 and alignment size should be equal.\
                          Got {len(block_df)+1=}, {len(alignment)=}')

    block_df['first'] = 's'
    block_df[['label', 'alignstring']] = pd.DataFrame(alignment[1:])
    block_df['label'] = block_df['label'].apply(lambda x: walks[x].genome)
    block_df['size'] = block_df['end'] - block_df['start'] + 1
    block_df['srcSize'] = block_df['label'].apply(lambda x: var_graph.genome_lengths[x])
    block_df['label'] = block_df['label'].apply(lambda x: genome_idx_to_name[str(x)])
    block_df = block_df[['first', 'label', 'start', 'size', 'strand', 'srcSize', 'alignstring']]
    
    # add carrying path
    maf_file.write(f's carrying {var_graph.carrying_len_so_far} {len(alignment[0][1])+len(alignment[0][1])} + {var_graph.carrying_len} {alignment[0][1]}\n')
    block_df.to_csv(maf_file, sep=' ', index=None, header=None)
    maf_file.write('\n')
    var_graph.carrying_len_so_far += len(alignment[0][1])