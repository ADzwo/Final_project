import os
import sys
import argparse

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

def wga(graph_file_path):
    var_graph = Graph(graph_file_path)
    blocks = find_collinear_blocks(var_graph)
    print(f'Found {len(blocks)} blocks for graph {graph_file_path}.')
    save_blocks_to_gff(blocks, graph=var_graph) # , graph_file_path.split('.')[0]
    with open(f'vertex_sequences/{var_graph.name}.txt') as f:
        sequences = f.readlines()
    for block_nr, block in enumerate(blocks):
        print(f'\n ------------ BLOCK NR {block_nr} ------------')
        poa_align(block, var_graph, sequences, _match=2, _mismatch=-2, _gap=-1)

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
    args = parser.parse_args()
    set_config(args.a, args.b, args.m, args.s)
    from tuples import *
    from graph import Graph
    from block_tools import save_blocks_to_gff
    from find_blocks import find_collinear_blocks
    wga(args.i)


# nr_blocks = []
# for graph_file_path in os.listdir('data'):
#     for SORT_SEEDS in ['nr_occurrences', 'length', 'no'][:2]:
#         print(f'{graph_file_path.upper()}, {SORT_SEEDS=}')
#         var_graph = Graph('data/'+graph_file_path)
#         blocks = find_collinear_blocks(var_graph)
#         nr_blocks.append(len(blocks))
#         save_blocks_to_gff(blocks, graph_file_path.split('.')[0], var_graph)
#         with open(f'vertex_sequences/{var_graph.name}.txt') as f:
#             sequences = f.readlines()

#         for block_nr, block in enumerate(blocks):
#             print(f'\n ------------ BLOCK NR {block_nr} ------------')
#             for walk in block.collinear_walks:
#                 assert var_graph.genomes[walk.genome].path[walk.start].vertex in block.carrying_path, f'{var_graph.genomes[walk.genome].path[walk.start].vertex=}'
#                 assert var_graph.genomes[walk.genome].path[walk.end].vertex in block.carrying_path, f'{var_graph.genomes[walk.genome].path[walk.end].vertex=}'
            
#             poa_align(block, var_graph, sequences, _match=2, _mismatch=-2, _gap=-1)
# print(f'Number of blocks for consecutive options:\n{nr_blocks}')

# # additional check
# SORT_SEEDS = 'no'
# for graph_file_path in os.listdir('data'):
#     for i in range(10):
#         print(f'{graph_file_path.upper()}, {SORT_SEEDS=}')
#         var_graph = Graph('data/'+graph_file_path)
#         blocks = find_collinear_blocks(var_graph)
#         nr_blocks.append(len(blocks))
#         save_blocks_to_gff(blocks, graph_file_path.split('.')[0], var_graph)
#         with open(f'vertex_sequences/{var_graph.name}.txt') as f:
#             sequences = f.readlines()

#         for block_nr, block in enumerate(blocks):
#             print(f'BLOCK NR {block_nr}')
#             for walk in block.collinear_walks:
#                 assert var_graph.genomes[walk.genome].path[walk.start].vertex in block.carrying_path, f'{var_graph.genomes[walk.genome].path[walk.start].vertex=}'
#                 assert var_graph.genomes[walk.genome].path[walk.end].vertex in block.carrying_path, f'{var_graph.genomes[walk.genome].path[walk.end].vertex=}'
#             poa_align(block, var_graph, sequences, _match=2, _mismatch=-2, _gap=-1)
# print(f'Number of blocks for consecutive options:\n{nr_blocks}')