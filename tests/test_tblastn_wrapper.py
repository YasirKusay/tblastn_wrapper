import unittest

from tblastn_wrapper.tblastn_wrapper import add

class TestHelperFunctions(unittest.TestCase):
    def test_add(self):
        self.assertEqual(add(1, 1), 2)
