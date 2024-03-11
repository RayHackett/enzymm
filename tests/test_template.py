import unittest

from matcher import template

class TestTemplate(unittest.TestCase):

    def test_load_template_missing_folder(self):
        with self.assertRaises(FileNotFoundError):
            templates = template.load_templates("/some/bogus/folder")
