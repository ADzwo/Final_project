import unittest
import os, sys
from tuples import *
from block import CollinearBlock, BlockExtensions
from graph import Graph

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

def match_carrying_check_block(block):
    for w_idx, walk in enumerate(block.collinear_walks):
            matches = block.match_carrying[w_idx]
            match_carrying_check(walk, matches)

def check_walk_orientations(graph, block):
    for walk, matches in zip(block.collinear_walks, block.match_carrying):
        g_idx = walk.genome
        g_path = graph.genomes[g_idx].path
        for match_carrying, match_walk in matches:
            g_pos = g_path[match_walk]
            carrying_orientation = block.carrying_path_orientations[match_carrying]
            assert g_pos.orientation*walk.orient==carrying_orientation


class TestCollinearBlock(unittest.TestCase):
    def setUp(self):
        seed_occurrences = [CollinearWalk(0, 1, 1, -1)]
        self.block1 = CollinearBlock(seed_idx=0, 
                                    seed_occurrences=seed_occurrences, 
                                    carrying_seed_orientation=1)

    def test_init(self):        
        self.assertEqual(self.block1.collinear_walks[0].orient, -1)
        self.assertEqual(self.block1.carrying_path_orientations[0], 1)
        self.assertEqual(len(self.block1.collinear_walks), 1)
        self.assertEqual(self.block1.match_carrying, [[(0, 1)]])
        self.assertEqual(self.block1.scores, [(0, 0)])
        

    def test_remove_doubled_matches(self):
        # example where no matches are removed
        self.block1.carrying_path = [0, 1, 2, 3, 0]
        self.block1.carrying_path_orientations = [1, -1, 1, 1, -1]
        self.block1.collinear_walks.append(CollinearWalk(1, 10, 20, 1))
        self.block1.match_carrying.append([(1,10), (3,12), (4,20)])
        match_carrying_check_block(self.block1)
        self.block1.remove_doubled_matches()
        match_carrying_check_block(self.block1)
            # check that there is no change
        self.assertEqual(len(self.block1.collinear_walks), 2)
        self.assertEqual(len(self.block1.match_carrying), 2)
        self.assertEqual(self.block1.match_carrying[0], [(0,1)])
        self.assertEqual(self.block1.match_carrying[1], [(1,10), (3,12), (4,20)])
        
        # example where matches are removed and walk.orient is 1
        self.block1.match_carrying[1] = [(1,10), (1,11), (3,12), (3,13), (3,14), (3,15), (4,18), (4,20)]
        match_carrying_check_block(self.block1)
        self.block1.remove_doubled_matches()
        match_carrying_check_block(self.block1)
            # check that there is no change
        self.assertEqual(len(self.block1.collinear_walks), 2)
        self.assertEqual(len(self.block1.match_carrying), 2)
        self.assertEqual(self.block1.match_carrying[0], [(0,1)])
        self.assertEqual(self.block1.match_carrying[1], [(1,10), (3,15), (4,20)])

        # example where matches are removed and walk.orient is -1
        self.block1.collinear_walks[1] = CollinearWalk(1, 10, 20, -1)
        self.block1.match_carrying[1] = [(1,20), (1,14), (1,13), (3,11), (4,12), (4,10)]
        match_carrying_check_block(self.block1)
        self.block1.remove_doubled_matches()
        match_carrying_check_block(self.block1)
            # check that there is no change
        self.assertEqual(len(self.block1.collinear_walks), 2)
        self.assertEqual(len(self.block1.match_carrying), 2)
        self.assertEqual(self.block1.match_carrying[0], [(0,1)])
        self.assertEqual(self.block1.match_carrying[1], [(1,20), (3,11), (4,10)])


class TestBlockExtensions(unittest.TestCase):
    def setUp(self):
        self.walks1 = [CollinearWalk(0, 0, 0, 1), CollinearWalk(1, 2, 3, -1)]
        self.walks2 = [CollinearWalk(0, 0, 1, -1), CollinearWalk(1, 1, 2, 1)]
        self.graph1 = Graph(data_path+'test_data1.gfa')
        self.graph2 = Graph(data_path+'test_data2.gfa')
        self.vertex_dict_path = f'{src}/vertex_name_to_idx/'
        self.genome_dict_path = f'{src}/genome_name_to_idx/'
        self.vertex_seq_path = f'{src}/vertex_sequences/'

    def tearDown(self):
        for graph in [self.graph1, self.graph2]:
            assert os.path.exists(f'{self.vertex_dict_path}/{graph.name}.json')
            os.remove(f'{self.vertex_dict_path}/{graph.name}.json')
            assert os.path.exists(f'{self.genome_dict_path}/{graph.name}.json')
            os.remove(f'{self.genome_dict_path}/{graph.name}.json')
            assert os.path.exists(f'{self.vertex_seq_path}/{graph.name}.txt')
            os.remove(f'{self.vertex_seq_path}/{graph.name}.txt')
    
    def test_init(self):
        # extensions = BlockExtensions(self.walks1, self.graph1, w0_idx=9)
        pass

if __name__=='__main__':
    unittest.main()