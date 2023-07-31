import os
import sys
import argparse
import json
import pandas as pd
from tuples import *
from graph import Graph
from block_tools import get_file_name, sort_v_idx, save_block_to_gff
from find_blocks import find_collinear_block

os.chdir(sys.path[0])
os.chdir('../')
src = os.getcwd()
print(f'PWD={src}') # Reading_graphs
assert 'data' in os.listdir()
maf_path = f'{src}/maf/'

assert os.path.exists(os.getcwd()+'/../poapy/'), 'Move to the root directory (the one containing README).'
sys.path.append(os.getcwd()+'/../poapy/')
from poa import poa_align

for folder in ['blocks', 'vertex_name_to_idx', 'genome_idx_to_name', 'vertex_sequences']:
    if not os.path.exists(folder):
        os.mkdir(folder)

def save_maf(alignment, maf_file, block_df, genome_lengths, walks, genome_idx_to_name):
    maf_file.write('a\n')
    # carrying_label, carrying_string = alignment[0]
    # maf_file.write(f's\t{carrying_label}\t0\t{len(carrying_string)-1}\t+\t{len(carrying_string)-1}\t{carrying_string}\n')

    if len(block_df)+1!=len(alignment):
        raise ValueError(f'Block size + 1 and alignment size should be equal. Got {len(block_df)+1=}, {len(alignment)=}')

    block_df['first'] = 's'
    block_df[['label', 'alignstring']] = pd.DataFrame(alignment[1:])
    block_df['label'] = block_df['label'].apply(lambda x: walks[x].genome)
    block_df['size'] = block_df['end'] - block_df['start'] + 1
    block_df['srcSize'] = block_df['label'].apply(lambda x: genome_lengths[x])
    block_df['label'] = block_df['label'].apply(lambda x: genome_idx_to_name[str(x)])
    block_df = block_df[['first', 'label', 'start', 'size', 'strand', 'srcSize', 'alignstring']]
    
    block_df.to_csv(maf_file, sep=' ', index=None, header=None)
    maf_file.write('\n')

def wga(graph_file_path, SORT_SEEDS, align, _match, _mismatch, _gap, a, b, m):
    var_graph = Graph(graph_file_path)
    name = get_file_name(var_graph.name, SORT_SEEDS, 'maf')
    maf_file = open(f'{maf_path}{name}', 'w')
    print(f'maf_file = {maf_path}{name}')
    maf_file.write('##maf version=1 scoring=tba.v8\n\n')
    vertex_indices_order = sort_v_idx(var_graph, SORT_SEEDS)
    block_nr = 0
    with open(f'vertex_sequences/{var_graph.name}.txt') as f:
        sequences = f.readlines()
    with open(f'genome_idx_to_name/{var_graph.name}.json', 'r') as f:
            genome_idx_to_name = json.load(f)
    for v_idx in vertex_indices_order: # select a vertex --- seed of a new CollinearBlock
        block = find_collinear_block(var_graph, v_idx, PARAM_a=a, PARAM_b=b, PARAM_m=m)
        if block is not None: # save block and its alignment
            block_nr += 1
            block_df = save_block_to_gff(block, graph=var_graph, SORT_SEEDS=SORT_SEEDS, block_nr=block_nr)
            if align==True:
                print(f'\n ------------ ALIGN BLOCK NR {block_nr} ------------')
                alignment = poa_align(block, var_graph, sequences, _match=_match, _mismatch=_mismatch, _gap=_gap)
                save_maf(alignment, maf_file, block_df, var_graph.genome_lengths, block.collinear_walks, genome_idx_to_name)
                del block_df, block, alignment
    maf_file.close()
    print(f'Found {block_nr} blocks for graph {graph_file_path}.')


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='Realative path to input .gfa file, from /Reading_graphs/ folder.')
    parser.add_argument('-a', default=150, required=False, 
                        help='Aundance pruning parameter, default to 150.\
                        Vertices occurring more than a times are not considered as collinear walk seeds.')
    parser.add_argument('-b', default=200, required=False, 
                        help='Maximal bubble size (for each walk), default to 200 residues.')
    parser.add_argument('-m', default=50, required=False,
                        help='Minimal collinear walk length, default to 50.')
    parser.add_argument('-s', default='no', required=False,
                        help="Seed sorting mode, default to 'no'. Possible values: \
                        'no' (random vertex order), 'nr_occurrences' (most occurrences first), \
                        'length' (longest labels first).")
    parser.add_argument('--align', action='store_true',
                        help='Use to align sequences within blocks. Default.')
    parser.add_argument('--no-align', dest='align', action='store_false',
                        help='Use not to align sequences within blocks. Tool aligns by default.')
    parser.set_defaults(align=True)
    parser.add_argument('--match', default=2, required=False,
                        help='Match score in alignment (used if --align is True). Default to 2.')
    parser.add_argument('--mismatch', default=-2, required=False,
                        help='Mismatch penalty in alignment (used if --align is True). Default to -2.')
    parser.add_argument('--gap', default=-1, required=False,
                        help='Gap penalty in alignment (used if --align is True). Default to -1.')
    args = parser.parse_args()
    
    if args.align==True:
        if not os.path.exists(maf_path):
            os.mkdir(maf_path)
    wga(args.i, args.s, args.align, args.match, args.mismatch, args.gap, 
        a=int(args.a), b=int(args.b), m=int(args.m))


