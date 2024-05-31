import unittest
import json

from matcher import utils

class TestUtils(unittest.TestCase):

    def test_ranked_argsort(self):
        self.assertEqual(utils.ranked_argsort([0, 4, 8, 6]), [1, 2, 4, 3])
        self.assertEqual(utils.ranked_argsort([2, 3, 20, 9]), [1, 2, 4, 3])
        self.assertEqual(utils.ranked_argsort([-3, 3, 20, 9]), [1, 2, 4, 3])
        self.assertEqual(utils.ranked_argsort([2, 3, 20, 20, 9]), [1, 2, 4, 4, 3])
        self.assertEqual(utils.ranked_argsort([2, 20, 3, 20, 9]), [1, 4, 2, 4, 3])

    def test_chunks(self):
        self.assertEqual(list(utils.chunks(['a', 'b', 'c'], 2)), [('a', 'b'), ('c', )])
        self.assertEqual(list(utils.chunks(['a', 'b', 'c'], 1)), [('a',), ('b',), ('c',)])
        self.assertEqual(len(list(utils.chunks(set(['a', 'b', 'b', 'c']), 2))), 2)
        self.assertEqual(list(utils.chunks({1: 'a', 2: 'b', 3: 'c'}, 2)), [(1,2), (3,)])
        
    # def test_request_url(self): no idea how to test this tbh

    def test_json_extract(self):
        self.assertEqual(utils.json_extract({1: 'a', 2: 'b', 3: {1: 'c', 2: 'd'}}, 2), ['b', 'd'])
        self.assertEqual(utils.json_extract([{1: 'a', 2: 'b'}, {1: 'c', 2: 'd'}], 2), ['b', 'd'])

class TestSetEncoder(unittest.TestCase):
    def test_encode_set(self):
        result = json.dumps({'numbers': {1, 2, 3, 4, 5}}, cls=utils.SetEncoder)
        self.assertEqual(result, '{"numbers": [1, 2, 3, 4, 5]}')