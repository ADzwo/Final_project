import unittest
from params import *

class TestParams(unittest.TestCase):
    def test_numeric_params(self):
        self.assertTrue(PARAM_m>0, msg=f'PARAM_m must be positive, got {PARAM_m=}')
        self.assertTrue(PARAM_b>0, msg=f'PARAM_b must be positive, got {PARAM_b=}')
        self.assertTrue(PARAM_a>=0, msg=f'PARAM_a must be nonnegative, got {PARAM_a=}')

    def test_SORT_SEEDS(self):
        self.assertIn(SORT_SEEDS, {'no', 'nr_occurrences', 'length'}, 
                      msg=f'SORT_SEEDS must be one of "no", "nr_occurrences", "length", got {SORT_SEEDS}')



if __name__=='__main__':
    unittest.main()