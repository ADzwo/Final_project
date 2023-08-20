import os
import sys
import argparse
import json
from tuples import *
from graph import Graph
from block_tools import get_file_name, sort_v_idx, save_block_to_gff, save_maf
from find_blocks import find_collinear_block

os.chdir(sys.path[0])
os.chdir('../')
src = os.getcwd()
# print(f'PWD={src}') # Reading_graphs
maf_path = f'{src}/maf/' # folder to save the .maf files (WGA)
block_path = f'{src}/blocks/' # folder to save the .gff files (block coordinates)

assert os.path.exists(os.getcwd()+'/../poapy/'), 'Move to the root directory (the one containing README).'
sys.path.append(os.getcwd()+'/../poapy/')
from poa import poa_align

for folder in ['blocks', 'vertex_name_to_idx', 'genome_idx_to_name', 'vertex_sequences']:
    if not os.path.exists(folder):
        os.mkdir(folder)

def wga(graph_file_path, SORT_SEEDS, align, _match, _mismatch, _gap, a, b, m):
    var_graph = Graph(graph_file_path)
    print('Built graph')
    vertex_indices_order = sort_v_idx(var_graph, SORT_SEEDS)
    with open(f'vertex_sequences/{var_graph.name}.txt') as f:
        sequences = f.readlines()
    with open(f'genome_idx_to_name/{var_graph.name}.json', 'r') as f:
        genome_idx_to_name = json.load(f)

    maf_file_name = get_file_name(var_graph.name, SORT_SEEDS, 'maf')
    maf_file = open(f'{maf_path}{maf_file_name}', 'w')
    maf_file.write('##maf version=1 scoring=tba.v8\n\n')

    block_nr = 0
    block_file_name = get_file_name(var_graph.name, SORT_SEEDS, 'gff')
    block_file = open(f'{block_path}{block_file_name}', 'w')

    for v_idx in vertex_indices_order: # select a vertex --- seed of a new CollinearBlock
        block = find_collinear_block(var_graph, v_idx, PARAM_a=a, PARAM_b=b, PARAM_m=m)
        if block is not None: # save block and its alignment
            block_nr += 1
            block_df = save_block_to_gff(block, graph=var_graph, block_nr=block_nr, file=block_file)
            # print(f'Found block nr {block_nr}.')
            if align==True:
                # print(f'\tAligning.')
                alignment = poa_align(block, var_graph, sequences, _match=_match, _mismatch=_mismatch, _gap=_gap)
                save_maf(alignment, maf_file, block_df, var_graph, block.collinear_walks, genome_idx_to_name)
    maf_file.close()
    block_file.close()
    # print(f'Found {block_nr} blocks for graph {graph_file_path}.')


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
    
    if args.s not in {'no', 'length', 'nr_occurrences'}:
        raise ValueError("Argument -s must be one of 'no' (default), 'length', 'nr_occurrences'.")
    wga(args.i, args.s, args.align, args.match, args.mismatch, args.gap, 
        a=int(args.a), b=int(args.b), m=int(args.m))


