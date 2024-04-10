import unittest

from matcher import utils

class TestUtils(unittest.TestCase):

    def test_ranked_argsort(self):
        self.assertEqual(utils.ranked_argsort([0, 4, 8, 6]), [1, 2, 4, 3])
        self.assertEqual(utils.ranked_argsort([2, 3, 20, 9]), [1, 2, 4, 3])
        self.assertEqual(utils.ranked_argsort([-3, 3, 20, 9]), [1, 2, 4, 3])
        self.assertEqual(utils.ranked_argsort([2, 3, 20, 20, 9]), [1, 2, 4, 4, 3])
        self.assertEqual(utils.ranked_argsort([2, 20, 3, 20, 9]), [1, 4, 2, 4, 3])

    def test_chunks(self):
        self.assertEqual(list(utils.chunks(['a', 'b', 'c'], 2)), [['a', 'b'], ['c']])
        self.assertEqual(list(utils.chunks(['a', 'b', 'c'], 1)), [['a'],['b'],['c']])
        self.assertEqual(len(list(utils.chunks(set(['a', 'b', 'b', 'c']), 2))), 2)
        self.assertEqual(list(utils.chunks({1: 'a', 2: 'b', 3: 'c'}, 2)), [[1,2], [3]])

    def test_convert_sets_to_lists(self):
        self.assertEqual(utils.convert_sets_to_lists('abc'), 'abc')
        self.assertEqual(utils.convert_sets_to_lists(12), 12)
        self.assertEqual(utils.convert_sets_to_lists([1, True, 'a']), [1, True, 'a'])
        self.assertEqual(utils.convert_sets_to_lists({1, 2, 3}), [1, 2, 3])
        self.assertEqual(utils.convert_sets_to_lists(('a', 'b', 'b')), ['a', 'b', 'b'])
        self.assertEqual(utils.convert_sets_to_lists([{1, 2, 3}, ['a', {0, 'b'}]]), [[1, 2, 3], ['a', [0, 'b']]])
        self.assertEqual(utils.convert_sets_to_lists(([1, 3, 4], (1, 3, 4), {1, 3, 4})), [[1, 3, 4], [1, 3, 4], [1, 3, 4]])
        self.assertEqual(utils.convert_sets_to_lists({1: {1: {0, 1}, 2: [0, 1]}, 2: [{2, 4}, (1, 3, 'a')]}), {1: {1: [0, 1], 2: [0, 1]}, 2: [[2,4], [1, 3, 'a']]})
        
    # def test_request_url(self): no idea how to test this tbh

    def test_json_extract(self):
        self.assertEqual(utils.json_extract({1: 'a', 2: 'b', 3: {1: 'c', 2: 'd'}}, 2), ['b', 'd'])