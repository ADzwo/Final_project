import unittest
from unittest.mock import patch
import os, sys
from tuples import *
import block
from block import CollinearBlock, BlockExtensions
from graph import Graph
from find_blocks import *

os.chdir(sys.path[0])
os.chdir('../')
src = os.getcwd()
print(f'{src=}')
assert 'test_data' in os.listdir()
data_path = src+'/test_data/'

def walk_start_end_check(walk, genome_length):
    assert walk.end<genome_length, f'end >= genome_length; {walk=}, {genome_length=}'
    assert walk.start>=0, f'start < 0; {walk=}, {genome_length-1=}'
    assert walk.start<=walk.end, f'end < start; {walk=}'

def match_carrying_check(walk, matches):
    last_match_carrying = 0
    for match_carrying, match_walk in matches:
        assert walk.start<=match_walk<=walk.end, f'Matching positions on walk should be between walk.start and walk.end {walk=}, {match_walk=}'
        assert match_carrying>=last_match_carrying, f'Carrying path indices for consecutive matches should not decrease. Got {last_match_carrying=}, {match_carrying=}.'
        last_match_carrying = match_carrying
    assert matches[0][1] in {walk.start, walk.end}, f'Walk must start on carrying path. Got {walk=}, {matches[0]=}.'
    assert matches[-1][1] in {walk.start, walk.end}, f'Walk must end on carrying path. Got {walk=}, {matches[-1]=}.'

def check_block(graph, block):
    for w_idx, (walk, matches) in enumerate(zip(block.collinear_walks, block.match_carrying)):
        g_idx = walk.genome
        g_path = graph.genomes[g_idx].path
        walk_start_end_check(walk, len(g_path))
        match_carrying_check(walk, matches)
        for match_carrying, match_walk in matches:
            g_pos = g_path[match_walk]
            carrying_orientation = block.carrying_path_orientations[match_carrying]
            assert g_pos.orientation*walk.orient==carrying_orientation, f'Orientation inconsistency for {w_idx=}.'


class TestFunctions(unittest.TestCase):
    def setUp(self):
        self.block1 = CollinearBlock(seed_idx=1, 
                                     seed_occurrences=[CollinearWalk(0, 1, 1, 1),
                                                       CollinearWalk(2, 2, 2, -1),
                                                       CollinearWalk(3, 1, 1, 1)],
                                     carrying_seed_orientation=1)
        self.block1.collinear_walks = [CollinearWalk(0, 1, 4, 1),
                                       CollinearWalk(2, 2, 3, -1),
                                       CollinearWalk(3, 1, 2, 1),
                                       CollinearWalk(1, 3, 3, -1)]
        self.block1.match_carrying = [[(0, 1), (1, 4)], [(0, 2), (1, 3)],
                                      [(0, 1), (2, 2)], [(4, 3)]]
        self.block1.scores.append(Score(0, 0))
        self.block1.carrying_path = [1, 4, 2, 9]
        self.block1.carrying_path_orientations = [1, 1, 1, 1]
        self.block2 = CollinearBlock(seed_idx=1, 
                                     seed_occurrences=[CollinearWalk(0, 1, 1, -1),
                                                       CollinearWalk(2, 2, 2, 1),
                                                       CollinearWalk(3, 1, 1, -1)],
                                     carrying_seed_orientation=1)
        self.block2.collinear_walks.append(CollinearWalk(4, 1, 1, 1))
        self.block2.match_carrying = [[(0, 1)], [(0, 2)],
                                      [(0, 1)], [(1, 1)]]
        self.block2.scores.append(Score(0, 0))
        self.block2.carrying_path = [1, 0]
        self.block2.carrying_path_orientations = [-1, 1]

    def test_merge_match_carrying(self):
        matches_f = [(i, j) for i, j in zip(range(5), range(5))]
        matches_b = [(i, j) for i, j in zip(range(5), range(5))]
        matches_final = merge_match_carrying(5, matches_f=matches_f, matches_b=matches_b)
        self.assertEqual(len(matches_final), 9)
        self.assertEqual([m[0] for m in matches_final], list(range(9)))
        self.assertEqual([m[1] for m in matches_final], list(range(4, 0, -1))+list(range(5)))

        matches_final = merge_match_carrying(5, matches_f=matches_f, matches_b=None)
        self.assertEqual(len(matches_final), 5)
        self.assertEqual([m[0] for m in matches_final], list(range(4,9)))
        self.assertEqual([m[1] for m in matches_final], list(range(5)))

        matches_final = merge_match_carrying(5, matches_f=None, matches_b=matches_b)
        self.assertEqual(len(matches_final), 5)
        self.assertEqual([m[0] for m in matches_final], list(range(5)))
        self.assertEqual([m[1] for m in matches_final], list(range(4, -1, -1)))


class TestFunctionsGraph(unittest.TestCase):
    def setUp(self):
        self.graph1 = Graph(data_path+'test_data1.gfa')
        self.vertex_dict_path = f'{src}/vertex_name_to_idx/'
        self.genome_dict_path = f'{src}/genome_idx_to_name/'
        self.vertex_seq_path = f'{src}/vertex_sequences/'
        
        self.block1 = CollinearBlock(seed_idx=1, 
                                     seed_occurrences=[CollinearWalk(0, 1, 1, 1),
                                                       CollinearWalk(2, 2, 2, -1),
                                                       CollinearWalk(3, 1, 1, 1)],
                                     carrying_seed_orientation=1)
        self.block1.collinear_walks = [CollinearWalk(0, 1, 4, 1),
                                       CollinearWalk(2, 2, 3, -1),
                                       CollinearWalk(3, 1, 2, 1),
                                       CollinearWalk(1, 3, 3, -1)]
        self.block1.match_carrying = [[(0, 1), (1, 4)], [(0, 2), (1, 3)],
                                      [(0, 1), (2, 2)], [(3, 3)]]
        self.block1.scores.append(Score(0, 0))
        self.block1.carrying_path = [1, 4, 2, 9]
        self.block1.carrying_path_orientations = [1, 1, 1, 1]
        self.block2 = CollinearBlock(seed_idx=1, 
                                     seed_occurrences=[CollinearWalk(0, 1, 1, -1),
                                                       CollinearWalk(2, 2, 2, 1),
                                                       CollinearWalk(3, 1, 1, -1)],
                                     carrying_seed_orientation=1)
        self.block2.collinear_walks.append(CollinearWalk(4, 1, 1, 1))
        self.block2.match_carrying = [[(0, 1)], [(0, 2)],
                                      [(0, 1)], [(1, 1)]]
        self.block2.scores.append(Score(0, 0))
        self.block2.carrying_path = [1, 0]
        self.block2.carrying_path_orientations = [-1, 1]

        self.block3 = CollinearBlock(seed_idx=4, 
                                     seed_occurrences=[CollinearWalk(0, 4, 4, 1)],
                                     carrying_seed_orientation=1)
    
    def tearDown(self):
        name = self.graph1.name
        assert os.path.exists(f'{self.vertex_dict_path}/{name}.json')
        os.remove(f'{self.vertex_dict_path}/{name}.json')
        assert os.path.exists(f'{self.genome_dict_path}/{name}.json')
        os.remove(f'{self.genome_dict_path}/{name}.json')
        assert os.path.exists(f'{self.vertex_seq_path}/{name}.txt')
        os.remove(f'{self.vertex_seq_path}/{name}.txt')

    def test_find_collinear_blocks(self):
        for v_idx in range(len(self.graph1.vertices)):
            final_block = find_collinear_block(self.graph1, v_idx=v_idx, PARAM_a=300, PARAM_b=300, PARAM_m=50)
            if final_block is not None:
                check_block(self.graph1, final_block)        

    def test_merge_forward_backward_blocks(self):
        nr_seeds = 3
        walk_nr = len(self.block1.match_carrying) + len(self.block2.match_carrying) - nr_seeds
        block_dict = {1:self.block1, -1:self.block2}
        final_block = merge_forward_backward_blocks(block_dict, nr_seeds)
        for walk, matches in zip(final_block.collinear_walks, final_block.match_carrying):
            match_carrying_check(walk, matches)
        self.assertEqual(len(final_block.match_carrying), walk_nr,
                         msg=f'Number of match_carrying elements should be equal to {walk_nr}.\
                            Got {len(final_block.match_carrying)=}.')
        self.assertEqual(len(final_block.collinear_walks), walk_nr,
                         msg=f'Number of collinear walks should be equal to {walk_nr}.\
                            Got {len(final_block.collinear_walks)=}')
        carrying_len = len(self.block1.carrying_path) + len(self.block2.carrying_path) - 1
        self.assertEqual(carrying_len, len(final_block.carrying_path_orientations),
                         msg=f'Carrying path orientations length should be equal to {carrying_len}.\
                            Got {len(final_block.carrying_path_orientations)=}.')
        self.assertEqual(len(final_block.carrying_path), carrying_len,
                         msg=f'Carrying path length should be equal to {carrying_len}.\
                            Got {len(final_block.carrying_path)=}.')

    def test_merge_blocks_and_postprocess(self):
        block_dict = {1:self.block1, -1:self.block2}
        final_block = merge_blocks_and_postprocess(block_dict, self.graph1, nr_seeds=3, PARAM_m=50)
        for matches in final_block.match_carrying:
            for i in range(len(matches)-1):
                self.assertTrue(matches[i][0]<matches[i+1][0], 
                                msg=f'Final match_carrying should have non-descending first elements.\
                                    Got {matches[i]=}, {matches[i+1]=} for {i=}.')
    
    def test_find_walk_to_extend(self):
        # no walk found
        e = BlockExtensions(self.block1.collinear_walks, self.graph1, 9, 1, PARAM_b=200)
        g_path = self.graph1.genomes[3].path
        walk_to_extend = find_walk_to_extend(self.block1, e, g_idx=3, g_path=g_path,
                                             o_nr_on_path=0, carrying_orient=-1)
        self.assertIsNone(walk_to_extend)

        # occurrence already in a walk
        e = BlockExtensions(self.block1.collinear_walks, self.graph1, 9, 1, PARAM_b=200)
        g_path = self.graph1.genomes[0].path
        walk_to_extend = find_walk_to_extend(self.block1, e, g_idx=0, g_path=g_path,
                                             o_nr_on_path=1, carrying_orient=1)
        self.assertEqual(walk_to_extend, -1)

        # one walk found
        e = BlockExtensions(self.block2.collinear_walks, self.graph1, 2, -1, PARAM_b=200)
        g_path = self.graph1.genomes[0].path
        walk_to_extend = find_walk_to_extend(self.block2, e, g_idx=0, g_path=g_path,
                                             o_nr_on_path=0, carrying_orient=-1)
        self.assertEqual(walk_to_extend, 0) # number of walk containing element nr 0 of genome 0

    def test_update_collinear_walks(self): # vertex_info - tuple(vertex index, orientation)
        self.block1.carrying_path_orientations[-1] = -1
        e = BlockExtensions(self.block1.collinear_walks, self.graph1, 9, 1, PARAM_b=200)
        update_collinear_walks(self.block1, CarryingPathExtension(9, -1), e, 
                               self.graph1, PARAM_b=200)
        self.assertIn(CollinearWalk(1, 0, 3, -1), self.block1.collinear_walks)
        self.assertIn([(3, 3), (3, 0)], self.block1.match_carrying)
        for walk, matches in zip(self.block1.collinear_walks, self.block1.match_carrying):
            match_carrying_check(walk, matches)
        check_block(self.graph1, self.block1)

        e = BlockExtensions(self.block3.collinear_walks, self.graph1, 2, -1, PARAM_b=200)
        # mark one occurrence of v_idx 5 as used (genome 2, position 0)
        g_path_pos = self.graph1.genomes[2].path[0]
        self.graph1.genomes[2].path[0] = Path(*(g_path_pos[:-2]), True, g_path_pos[-1])
        update_collinear_walks(self.block3, CarryingPathExtension(5, 1), e, 
                               self.graph1, PARAM_b=200)
        self.assertEqual(self.block3.collinear_walks, [CollinearWalk(0, 4, 5, 1)])
        self.assertEqual(self.block3.scores, [(0, 0)])
        for walk, matches in zip(self.block3.collinear_walks, self.block3.match_carrying):
            match_carrying_check(walk, matches)
        check_block(self.graph1, self.block3)

        
if __name__=='__main__':
    unittest.main()