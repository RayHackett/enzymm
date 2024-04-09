import unittest

from matcher import utils

class TestCommonFunctions(unittest.TestCase):

    def test_ranked_argsort(self):
        self.assertEqual(common_functions.ranked_argsort([0, 4, 8, 6]), [1, 2, 4, 3])
        self.assertEqual(common_functions.ranked_argsort([2, 3, 20, 9]), [1, 2, 4, 3])
        self.assertEqual(common_functions.ranked_argsort([-3, 3, 20, 9]), [1, 2, 4, 3])
        self.assertEqual(common_functions.ranked_argsort([2, 3, 20, 20, 9]), [1, 2, 4, 4, 3])
        self.assertEqual(common_functions.ranked_argsort([2, 20, 3, 20, 9]), [1, 4, 2, 4, 3])
