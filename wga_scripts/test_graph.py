import unittest
import os, sys
from tuples import *
from config import *
from graph import *

os.chdir(sys.path[0])
os.chdir('../')
src = os.getcwd()
print(f'{src=}')
assert 'test_data' in os.listdir()
data_path = src+'/test_data/'

def check_v_idx_in_genome_path(v_idx, v, graph):
    for g_idx, i in v.occurrences: # occurrences of the vertex
        genome = graph.genomes[g_idx]
        g_path_pos = genome.path[i] # vertex index: int, orientation: int, used: bool
        assert g_path_pos.vertex==v_idx, f'For occurrences of the same vertex, v_idx should be the same in genome path.\
                                            Got primary {v_idx=}, path_v_idx={g_path_pos.vertex}.'


class TestGraph(unittest.TestCase):
    def setUp(self):
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

    def test_v_idx_in_genome_path(self):
        # test_data 1 and 2
        for graph in [self.graph1, self.graph2]:
            for v_idx, v in enumerate(graph.vertices):
                check_v_idx_in_genome_path(v_idx, v, graph)
        # test_data 3
        with self.assertRaises(ValueError):
            Graph(data_path+'test_data3.gfa')
        os.remove(f'{self.vertex_seq_path}/test_data3.txt')
    
    def test_add_vertex(self):
        nr_vertices_old = len(self.graph1.vertices)
        # check the increase of the number of vertices and check the sequence length retrieved
        v_name, v_sequence = self.graph1.add_vertex(line='S 11 aaagtg')
        self.assertEqual(nr_vertices_old+1, len(self.graph1.vertices))
        self.assertEqual(self.graph1.vertices[-1].v_length, 6)
        self.assertEqual(v_name, '11')
        self.assertEqual(v_sequence.strip(), 'AAAGTG')
        # check that nothing is changed after trying to add an empty vertex
        v_name, v_sequence = self.graph1.add_vertex(line='S 12 ')
        self.assertEqual(nr_vertices_old+1, len(self.graph1.vertices))
        self.assertEqual(self.graph1.vertices[-1].v_length, 6)
        self.assertIsNone(v_name)
        self.assertIsNone(v_sequence)

        nr_vertices_old = len(self.graph2.vertices)
        # check the increase of the number of vertices and check the sequence length retrieved
        v_name, v_sequence = self.graph2.add_vertex(line='S agtc A')
        self.assertEqual(nr_vertices_old+1, len(self.graph2.vertices))
        self.assertEqual(self.graph2.vertices[-1].v_length, 1)
        self.assertEqual(v_name, 'agtc')
        self.assertEqual(v_sequence.strip(), 'A')
        # check that nothing is changed after trying to add an empty vertex
        v_name, v_sequence = self.graph2.add_vertex(line='S C')
        self.assertEqual(nr_vertices_old+1, len(self.graph2.vertices))
        self.assertEqual(self.graph2.vertices[-1].v_length, 1)
        self.assertIsNone(v_name)
        self.assertIsNone(v_sequence)

    def test_vertex_name_to_idx(self):
        for graph in [self.graph1, self.graph2]:
            with open(f'{self.vertex_dict_path}/{graph.name}.json') as f:
                vertex_name_to_idx = json.load(f)
            self.assertEqual(len(set(vertex_name_to_idx.values())), len(vertex_name_to_idx.values()),
                             msg=f'vertex_name_to_idx dict should have different values for different keys.\
                                  Check {self.vertex_dict_path}/{graph.name}.json')

            with open(f'{self.genome_dict_path}/{graph.name}.json') as f:
                genome_name_to_idx = json.load(f)
            self.assertEqual(len(set(genome_name_to_idx.values())), len(genome_name_to_idx.values()),
                             msg=f'genome_name_to_idx dict should have different values for different keys.\
                                  Check {self.genome_dict_path}/{graph.name}.json')

            with open(f'{self.vertex_seq_path}/{graph.name}.txt') as f:
                sequences = list(f.readlines())
            self.assertEqual(len(sequences), len(graph.vertices),
                             msg=f'Number of sequences should be equal to the number of vertices. \
                                Got {len(sequences)=}, {len(graph.vertices)=}')

    def test_find_seeds(self):
        collinear_seeds, carrying_seed_orientation = self.graph1.find_seeds(4)
        self.assertEqual(len(collinear_seeds), 2)
        self.assertEqual(collinear_seeds, [CollinearWalk(0, 4, 4, 1), CollinearWalk(2, 3, 3, -1)])
        self.assertEqual(carrying_seed_orientation, 1)

        collinear_seeds, carrying_seed_orientation = self.graph2.find_seeds(1)
        self.assertEqual(len(collinear_seeds), 1)
        self.assertEqual(collinear_seeds, [CollinearWalk(1, 1, 1, 1)])
        self.assertEqual(carrying_seed_orientation, -1)


if __name__=='__main__':
    unittest.main()
     
