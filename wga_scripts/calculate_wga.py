import os
import sys
from tuples import *
from params import *
from graph import Graph
from save import save_blocks_to_gff

os.chdir(sys.path[0])
os.chdir('../')
print(f'PWD={os.getcwd()}')
assert 'data' in os.listdir()

assert os.path.exists(os.getcwd()+'/../poapy/'), 'Move to the root directory (the one containing README).'
sys.path.append(os.getcwd()+'/../poapy/')
from poa import poa_align

for folder in ['blocks', 'vertex_name_to_idx', 'genome_name_to_idx', 'vertex_sequences', 'poa_visualized']:
    if not os.path.exists(folder):
        os.mkdir(folder)

nr_blocks = []
for graph_file_path in os.listdir('data'):
    for SORT_SEEDS in ['nr_occurrences', 'length', 'no'][:2]:
        print(f'{graph_file_path.upper()}, {SORT_SEEDS=}')
        var_graph = Graph('data/'+graph_file_path)
        blocks = var_graph.find_collinear_blocks()
        nr_blocks.append(len(blocks))
        save_blocks_to_gff(blocks, graph_file_path.split('.')[0], var_graph)
        # save_blocks(blocks, graph_file_path.split('.')[0], var_graph)
        # graph_name = re.split(r'[\.\/]', graph_file_path)[-2]
        with open(f'vertex_sequences/{var_graph.name}.txt') as f:
            sequences = f.readlines()

        for block_nr, block in enumerate(blocks):
            print(f'\n ------------ BLOCK NR {block_nr} ------------')
            for walk in block.collinear_walks:
                assert var_graph.genomes[walk.genome].path[walk.start].vertex in block.carrying_path, f'{var_graph.genomes[walk.genome].path[walk.start].vertex=}'
                assert var_graph.genomes[walk.genome].path[walk.end].vertex in block.carrying_path, f'{var_graph.genomes[walk.genome].path[walk.end].vertex=}'
            
            poa_align(block, var_graph, sequences, _match=2, _mismatch=-2, _gap=-1)
print(f'Number of blocks for consecutive options:\n{nr_blocks}')

# additional check
SORT_SEEDS = 'no'
for graph_file_path in os.listdir('data'):
    for i in range(10):
        print(f'{graph_file_path.upper()}, {SORT_SEEDS=}')
        var_graph = Graph('data/'+graph_file_path)
        blocks = var_graph.find_collinear_blocks()
        nr_blocks.append(len(blocks))
        save_blocks_to_gff(blocks, graph_file_path.split('.')[0], var_graph)
        # save_blocks(blocks, graph_file_path.split('.')[0], var_graph)
        # graph_name = re.split(r'[\.\/]', graph_file_path)[-2]
        with open(f'vertex_sequences/{var_graph.name}.txt') as f:
            sequences = f.readlines()

        for block_nr, block in enumerate(blocks):
            print(f'BLOCK NR {block_nr}')
            for walk in block.collinear_walks:
                assert var_graph.genomes[walk.genome].path[walk.start].vertex in block.carrying_path, f'{var_graph.genomes[walk.genome].path[walk.start].vertex=}'
                assert var_graph.genomes[walk.genome].path[walk.end].vertex in block.carrying_path, f'{var_graph.genomes[walk.genome].path[walk.end].vertex=}'
            poa_align(block, var_graph, sequences, _match=2, _mismatch=-2, _gap=-1)
print(f'Number of blocks for consecutive options:\n{nr_blocks}')