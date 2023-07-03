import os
import sys
import argparse
import pandas as pd

os.chdir(sys.path[0])
os.chdir('../')
src = os.getcwd()
print(f'PWD={src}') # Reading_graphs
assert 'data' in os.listdir()

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

def save_maf(alignment, maf_file, block_df):
    maf_file.write('a\n')
    carrying_label, carrying_string = alignment[0]
    maf_file.write(f's\t{carrying_label}\t0\t{len(carrying_string)-1}\t+\t{carrying_string}\n')

    assert len(block_df)+1==len(alignment), f'{len(block_df)=}, {len(alignment)=}'
    block_df['first'] = 's'
    block_df[['label', 'alignstring']] = pd.DataFrame(alignment[1:])
    block_df = block_df[['first', 'label', 'start', 'end', 'strand', 'alignstring']]
    
    block_df.to_csv(maf_file, sep='\t', index=None, header=None)
    maf_file.write('\n')

def wga(graph_file_path, SORT_SEEDS, align, _match, _mismatch, _gap):
    var_graph = Graph(graph_file_path)
    blocks = find_collinear_blocks(var_graph)
    print(f'Found {len(blocks)} blocks for graph {graph_file_path}.')
    blocks_df = save_blocks_to_gff(blocks, graph=var_graph, SORT_SEEDS=SORT_SEEDS) # , graph_file_path.split('.')[0]
    with open(f'vertex_sequences/{var_graph.name}.txt') as f:
        sequences = f.readlines()
    maf_file = open(f'{src}/maf/{var_graph.name}.maf', 'w')
    for block_nr, block in enumerate(blocks):
        print(f'\n ------------ BLOCK NR {block_nr} ------------')
        if align==True:
            alignment = poa_align(block, var_graph, sequences, _match=_match, _mismatch=_mismatch, _gap=_gap)
            save_maf(alignment, maf_file, blocks_df[block_nr])
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
    parser.add_argument('--align', default=True, required=False,
                        help='Indicates whether to align sequences within blocks (bool). Default to True.')
    parser.add_argument('--match', default=2, required=False,
                        help='Match score in alignment (used if --align is True). Default to 2.')
    parser.add_argument('--mismatch', default=-2, required=False,
                        help='Mismatch penalty in alignment (used if --align is True). Default to -2.')
    parser.add_argument('--gap', default=-1, required=False,
                        help='Gap penalty in alignment (used if --align is True). Default to -1.')
    args = parser.parse_args()
    set_config(args.a, args.b, args.m, args.s)
    from tuples import *
    from graph import Graph
    from block_tools import save_blocks_to_gff
    from find_blocks import find_collinear_blocks
    wga(args.i, args.s, args.align, args.match, args.mismatch, args.gap)


