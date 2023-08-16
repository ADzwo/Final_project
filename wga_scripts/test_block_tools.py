import unittest
import os, sys, math, json
from collections import namedtuple
import pandas as pd
from datetime import date
from tuples import *
from block_tools import *
from graph import Graph

os.chdir(sys.path[0])
os.chdir('../')
src = os.getcwd() + '/'
print(f'{src=}')
assert 'test_data' in os.listdir()
data_path = src+'test_data/'

today = str(date.today()).replace('-', '_')

class TestGraph(unittest.TestCase):
    def setUp(self):
        self.graph1 = Graph(data_path+'test_data1.gfa')
        self.graph2 = Graph(data_path+'test_data2.gfa')
        self.vertex_dict_path = f'{src}/vertex_name_to_idx/'
        self.genome_dict_path = f'{src}/genome_idx_to_name/'
        self.vertex_seq_path = f'{src}/vertex_sequences/'
    
    def tearDown(self):
        for graph in [self.graph1, self.graph2]:
            assert os.path.exists(f'{self.vertex_dict_path}/{graph.name}.json')
            os.remove(f'{self.vertex_dict_path}/{graph.name}.json')
            assert os.path.exists(f'{self.genome_dict_path}/{graph.name}.json')
            os.remove(f'{self.genome_dict_path}/{graph.name}.json')
            assert os.path.exists(f'{self.vertex_seq_path}/{graph.name}.txt')
            os.remove(f'{self.vertex_seq_path}/{graph.name}.txt')

    def test_sort_v_idx(self):
        vertices_sorted = sort_v_idx(graph=self.graph1, SORT_SEEDS='nr_occurrences')
        self.assertEqual(vertices_sorted[0], 0)
        
        vertices_sorted = sort_v_idx(graph=self.graph2, SORT_SEEDS='length')
        self.assertEqual(vertices_sorted[0], 0)
        self.assertEqual(vertices_sorted[-1], 2) # vertex with name '4' has idx 2 (checked below)
        with open(f'{self.vertex_dict_path}/{self.graph2.name}.json', 'r') as f:
            vertex_dict = json.load(f)
        self.assertEqual(vertex_dict['4'], 2)

    def test_walk_length(self):
        # one element walk
        walk1 = CollinearWalk(0, 0, 0, 1)
        len1 = walk_length(walk1, self.graph1)
        self.assertEqual(len1, 32)
        len2 = walk_length(walk1, self.graph2)
        self.assertEqual(len2, 23)
        # longer walk
        walk2 = CollinearWalk(1, 1, 2, -1)
        len1 = walk_length(walk2, self.graph1)
        self.assertEqual(len1, 36)
        len2 = walk_length(walk2, self.graph2)
        self.assertEqual(len2, 4)
        # start, end not None
        walk3 = CollinearWalk(1, 0, 0, -1)
        len1 = walk_length(walk3, self.graph1, 0, 3)
        self.assertEqual(len1, 44)
        len2 = walk_length(walk3, self.graph2, 0, 0)
        self.assertEqual(len2, 23)

        walk4 = CollinearWalk(1, 1, 3, -1)
        with self.assertRaises(ValueError):
            len1 = walk_length(walk4, self.graph1, 2, 1)
        
        with self.assertRaises(ValueError):
            len1 = walk_length(walk4, self.graph2)

        walk5 = CollinearWalk(1, 1, 3, -1)
        with self.assertRaises(ValueError):
            len1 = walk_length(walk5, self.graph1, -1, 1)

    def test_find_vertex_on_path_till_b(self):
        walk1 = CollinearWalk(0, 0, 9, 1)

        pos = find_vertex_on_path_till_b(self.graph1, walk1, 9, 1, 7, 9, PARAM_b=200)
        self.assertEqual(pos, 9)
        
        pos = find_vertex_on_path_till_b(self.graph1, walk1, 9, 1, 0, 9, PARAM_b=200)
        self.assertIsNone(pos) # vertex 9 is not reachable within PARAM_b from 0

        walk1 = CollinearWalk(0, 0, 9, -1)
        pos = find_vertex_on_path_till_b(self.graph1, walk1, 9, -1, 9, 0, PARAM_b=200)
        self.assertEqual(pos, 9)

    def test_mark_vertices_as_used(self):
        MockBlock = namedtuple('MockBlock', ['collinear_walks'])

        # one-element block
        block = MockBlock([CollinearWalk(0, 0, 9, -1)])
        mark_vertices_as_used(self.graph1, block)
        for pos in self.graph1.genomes[0].path:
            self.assertTrue(pos.used)

        # 2-element block
        block = MockBlock([CollinearWalk(0, 0, 6, -1), CollinearWalk(1, 1, 3, 1)])
        mark_vertices_as_used(self.graph1, block)
        for pos in self.graph1.genomes[0].path[:7]:
            self.assertTrue(pos.used)
        for pos in self.graph1.genomes[1].path[1:]:
            self.assertTrue(pos.used)
        # check that no more vertices are marked as used
        self.assertFalse(self.graph1.genomes[1].path[0].used)
        for i in range(2,5):
            for pos in self.graph1.genomes[i].path:
                self.assertFalse(pos.used)

        # 2 walks on the same genome
        block = MockBlock([CollinearWalk(1, 0, 1, -1), CollinearWalk(1, 1, 2, 1)])
        mark_vertices_as_used(self.graph2, block)
        for pos in self.graph2.genomes[1].path:
            self.assertTrue(pos.used)

    def test_remove_short_walks(self):
        MockBlock = namedtuple('MockBlock', 
                               ['collinear_walks', 'match_carrying'])
        
        # no need to remove walks
        block = MockBlock([CollinearWalk(0, 0, 9, -1)], [None])
        remove_short_walks(block, self.graph1, PARAM_m=50)
        self.assertEqual(len(block.collinear_walks), 1)
        self.assertEqual(len(block.match_carrying), 1)

        # 2 walks, one to remove
        block = MockBlock([CollinearWalk(0, 1, 8, -1), CollinearWalk(1, 1, 2, 1)],
                          [0, 1])
        remove_short_walks(block, self.graph1, PARAM_m=50)
        self.assertEqual(len(block.collinear_walks), 1)
        self.assertEqual(block.match_carrying, [0])

        # 2 short walks
        block = MockBlock([CollinearWalk(0, 0, 1, -1), CollinearWalk(1, 1, 2, 1)],
                          [0, 1])
        remove_short_walks(block, self.graph2, PARAM_m=50)
        self.assertEqual(len(block.collinear_walks), 0)
        self.assertEqual(len(block.match_carrying), 0)
    
    def test_scoring_function_big_param_m(self):
        MockBlock = namedtuple('MockBlock', 
                               ['collinear_walks', 'match_carrying', 'scores'])
        
        block = MockBlock([CollinearWalk(0, 0, 0, -1)], [None], Score(0, 0))
        score = scoring_function(block, self.graph1, PARAM_b=200, PARAM_m=50)
        self.assertEqual(score, 0)
        
        block = MockBlock([CollinearWalk(0, 0, 0, -1)], [None], [Score(201, 0)])
        score = scoring_function(block, self.graph1, PARAM_b=200, PARAM_m=50)
        self.assertEqual(score, 0)

        block = MockBlock([CollinearWalk(0, 0, 0, -1)], [None], [Score(0, 201)])
        score = scoring_function(block, self.graph1, PARAM_b=200, PARAM_m=50)
        self.assertEqual(score, 0)
        
    def test_scoring_function_small_params(self):
        MockBlock = namedtuple('MockBlock',
                               ['collinear_walks', 'match_carrying', 'scores'])
        # one-element block with different q1 and q3 values
        block = MockBlock([CollinearWalk(0, 0, 0, -1)], [None], [Score(0, 0)])
        score = scoring_function(block, self.graph1, PARAM_b=10, PARAM_m=5)
        self.assertEqual(score, 32)

        block = MockBlock([CollinearWalk(0, 0, 0, 1)], [None], [Score(1, 2)])
        score = scoring_function(block, self.graph1, PARAM_b=10, PARAM_m=5)
        self.assertEqual(score, 32-9)

        block = MockBlock([CollinearWalk(0, 0, 0, -1)], [None], [Score(11, 0)])
        score = scoring_function(block, self.graph1, PARAM_b=10, PARAM_m=5)
        self.assertTrue(math.isinf(score))

        block = MockBlock([CollinearWalk(0, 0, 0, -1)], [None], [Score(0, 11)])
        score = scoring_function(block, self.graph1, PARAM_b=10, PARAM_m=5)
        self.assertTrue(math.isinf(score))

        # check if block score is a sum of walk's scores
        block = MockBlock([CollinearWalk(0, 0, 0, -1), CollinearWalk(1, 0, 1, 1)],
                           [None, None], [Score(0, 0), Score(0, 10)])
        score = scoring_function(block, self.graph1, PARAM_b=10, PARAM_m=5)
        self.assertEqual(score, -38) # -38 = 32+26+4-100

        block = MockBlock([CollinearWalk(0, 0, 0, -1), CollinearWalk(1, 0, 1, 1)], 
                          [None, None], [Score(0, 0), Score(0, 11)])
        score = scoring_function(block, self.graph1, PARAM_b=10, PARAM_m=5)
        self.assertTrue(math.isinf(score))
    
class TestSave(unittest.TestCase):
    def setUp(self):
        self.graph1 = Graph(data_path+'test_data1.gfa')
        self.graph2 = Graph(data_path+'test_data2.gfa')
        self.vertex_dict_path = f'{src}/vertex_name_to_idx/'
        self.genome_dict_path = f'{src}/genome_idx_to_name/'
        self.vertex_seq_path = f'{src}/vertex_sequences/'
        self.blocks_path = f'{src}/blocks/'

        MockBlock = namedtuple('MockBlock', ['collinear_walks'])
        self.empty_block = MockBlock([])
        self.block1 = MockBlock([CollinearWalk(0, 0, 0, -1), CollinearWalk(1, 0, 1, 1)])

    
    def tearDown(self):
        for graph in [self.graph1, self.graph2]:
            assert os.path.exists(f'{self.vertex_dict_path}/{graph.name}.json')
            os.remove(f'{self.vertex_dict_path}/{graph.name}.json')
            assert os.path.exists(f'{self.genome_dict_path}/{graph.name}.json')
            os.remove(f'{self.genome_dict_path}/{graph.name}.json')
            assert os.path.exists(f'{self.vertex_seq_path}/{graph.name}.txt')
            os.remove(f'{self.vertex_seq_path}/{graph.name}.txt')
            for file in os.listdir(self.blocks_path):
                if file.startswith(graph.name):
                    if file.endswith(f'{today}.gff'):
                        os.remove(f'{self.blocks_path}/{graph.name}_{today}.gff')
                    if file.endswith('length.gff'):
                        os.remove(f'{self.blocks_path}/{graph.name}_{today}_sort_by_length.gff')
                    if file.endswith('ences.gff'):
                        os.remove(f'{self.blocks_path}/{graph.name}_{today}_sort_by_nr_occurrences.gff')

    def test_save_block_to_gff(self):
        for sort_seeds, suffix in [('no', ''), ('nr_occurrences', '_sort_by_nr_occurrences'), ('length', '_sort_by_length')]:
            file_name = get_file_name(self.graph2.name, SORT_SEEDS=sort_seeds, extension='gff')
            with open(src+'blocks/'+file_name, 'w') as file:
                save_block_to_gff(self.empty_block, graph=self.graph2, block_nr=0, file=file)
            block_df = pd.read_csv(f'{self.blocks_path}{self.graph2.name}_{today}{suffix}.gff', sep=',')
            self.assertEqual(len(block_df), 0)

            file_name = get_file_name(self.graph1.name, SORT_SEEDS=sort_seeds, extension='gff')
            with open(src+'blocks/'+file_name, 'w') as file:
                save_block_to_gff(self.block1, graph=self.graph1, block_nr=0, file=file)
            block_df = pd.read_csv(f'{self.blocks_path}{self.graph1.name}_{today}{suffix}.gff', sep=',')
            self.assertEqual(len(block_df), 2)
            self.assertTrue(block_df['seqname'].tolist(), ['0', '1'])
            self.assertTrue(block_df['source'].unique(), ['final_project'])
            self.assertTrue(block_df['feature'].unique(), ['.'])
            self.assertTrue(block_df['frame'].unique(), ['.'])
            self.assertTrue(block_df['score'].unique(), ['.'])
            self.assertTrue(block_df['attribute'].unique(), ['ID=0'])
            self.assertTrue(block_df['strand'].tolist(), ['-', '+'])
            self.assertTrue(block_df['start'].tolist(), ['0', '0'])
            self.assertTrue(block_df['end'].tolist(), ['32', '4'])

    def test_save_maf(self):
        alignment = [(-1, -1), (0, 0), (1, 1)]
        maf_file = open('tmp.maf', 'w')

        block_file_name = get_file_name(self.graph1.name, SORT_SEEDS='no', extension='gff')
        block_file = open(src+'blocks/'+block_file_name, 'w')

        block_df = save_block_to_gff(self.block1, graph=self.graph1, block_nr=0, file=block_file)
        with open(f'{self.genome_dict_path}/{self.graph1.name}.json', 'r') as f:
            genome_idx_to_name = json.load(f)
        save_maf(alignment, maf_file, block_df, self.graph1.genome_lengths, self.block1.collinear_walks, genome_idx_to_name)
        maf_file.close()
        os.remove('tmp.maf')

class TestGetFileName:

    def test_get_file_name(self):
        name = get_file_name('', SORT_SEEDS='no', extension='a')
        self.assertEqual(name, f'_{today}.a')
        name = get_file_name('graph_name', SORT_SEEDS='no', extension='a')
        self.assertEqual(name, f'graph_name_{today}.a')

        name = get_file_name('', SORT_SEEDS='length', extension='a')
        self.assertEqual(name, f'_{today}_sort_by_length.a')

        name = get_file_name('', SORT_SEEDS='nr_occurrences', extension='a')
        self.assertEqual(name, f'_{today}_sort_by_nr_occurrences.a')
        

if __name__=='__main__':
    unittest.main()
     