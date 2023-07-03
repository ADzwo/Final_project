import os
import sys
import argparse
import pandas as pd
from tuples import *
from graph import Graph
from block_tools import get_file_name, save_blocks_to_gff
from find_blocks import find_collinear_blocks

os.chdir(sys.path[0])
os.chdir('../')
src = os.getcwd()
print(f'PWD={src}') # Reading_graphs
assert 'data' in os.listdir()
maf_path = f'{src}/maf/'

assert os.path.exists(os.getcwd()+'/../poapy/'), 'Move to the root directory (the one containing README).'
sys.path.append(os.getcwd()+'/../poapy/')
from poa import poa_align

for folder in ['blocks', 'vertex_name_to_idx', 'genome_name_to_idx', 'vertex_sequences']:
    if not os.path.exists(folder):
        os.mkdir(folder)

def set_config(a, b, m, sort_seeds):
    with open(f'{src}/wga_scripts/config.py', 'w') as f:
        f.write(f'PARAM_a = {a}\n') # abundance pruning parameter; 150 used for mice dataset in SibeliaZ paper
        f.write(f'PARAM_b = {b}\n')
        f.write(f'PARAM_m = {m}\n')
        f.write(f"SORT_SEEDS = '{sort_seeds}'\n")

def save_maf(alignment, maf_file, block_df, genome_lengths):
    maf_file.write('a\n')
    carrying_label, carrying_string = alignment[0]
    maf_file.write(f's\t{carrying_label}\t0\t{len(carrying_string)-1}\t+\t{len(carrying_string)-1}\t{carrying_string}\n')

    if len(block_df)+1!=len(alignment):
        raise ValueError(f'Block size + 1 and alignment size should be equal. Got {len(block_df)+1=}, {len(alignment)=}')

    block_df['first'] = 's'
    block_df[['label', 'alignstring']] = pd.DataFrame(alignment[1:])
    block_df['size'] = block_df['end'] - block_df['start'] + 1
    block_df['srcSize'] = block_df['label'].apply(lambda x: genome_lengths[int(x)])
    block_df = block_df[['first', 'label', 'start', 'size', 'strand', 'srcSize', 'alignstring']]
    
    block_df.to_csv(maf_file, sep='\t', index=None, header=None)
    maf_file.write('\n')

def get_walk_genome_lengths(var_graph, collinear_walks):
    g_lens = [g.path[-1].p_length for g in var_graph.genomes]
    genome_lengths = []
    for walk in collinear_walks:
        genome_lengths.append(g_lens[walk.genome])
    return genome_lengths
    

def wga(graph_file_path, SORT_SEEDS, align, _match, _mismatch, _gap, a, b, m):
    var_graph = Graph(graph_file_path)
    blocks = find_collinear_blocks(var_graph, SORT_SEEDS, PARAM_a=a, PARAM_b=b, PARAM_m=m)
    print(f'Found {len(blocks)} blocks for graph {graph_file_path}.')
    blocks_df = save_blocks_to_gff(blocks, graph=var_graph, SORT_SEEDS=SORT_SEEDS)
    with open(f'vertex_sequences/{var_graph.name}.txt') as f:
        sequences = f.readlines()
    name = get_file_name(var_graph.name, SORT_SEEDS, 'maf')
    
    maf_file = open(f'{maf_path}{name}', 'w')
    for block_nr, block in enumerate(blocks):
        print(f'\n ------------ BLOCK NR {block_nr} ------------')
        if align==True:
            alignment = poa_align(block, var_graph, sequences, _match=_match, _mismatch=_mismatch, _gap=_gap)
            genome_lengths = get_walk_genome_lengths(var_graph, block.collinear_walks)
            save_maf(alignment, maf_file, blocks_df[block_nr], genome_lengths)
    maf_file.close()

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
    parser.set_defaults(feature=True)
    parser.add_argument('--match', default=2, required=False,
                        help='Match score in alignment (used if --align is True). Default to 2.')
    parser.add_argument('--mismatch', default=-2, required=False,
                        help='Mismatch penalty in alignment (used if --align is True). Default to -2.')
    parser.add_argument('--gap', default=-1, required=False,
                        help='Gap penalty in alignment (used if --align is True). Default to -1.')
    args = parser.parse_args()
    # set_config(args.a, args.b, args.m, args.s)
    
    if args.align==True:
        if not os.path.exists(maf_path):
            os.mkdir(maf_path)
    print(f'{args.align=}')
    wga(args.i, args.s, args.align, args.match, args.mismatch, args.gap, 
        a=int(args.a), b=int(args.b), m=int(args.m))


