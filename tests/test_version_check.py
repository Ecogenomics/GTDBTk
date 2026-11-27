import unittest

from gtdbtk.tools import is_newer_version


class TestVersionCheck(unittest.TestCase):

    def test_is_newer_version(self):
        self.assertTrue(is_newer_version('3.0.1', '3.0.0'))
        self.assertTrue(is_newer_version('2.10.0', '2.9.0'))     # not a string comparison
        self.assertFalse(is_newer_version('2.8.0', '3.0.0'))     # local ahead of PyPI
        self.assertFalse(is_newer_version('3.0.0', '3.0.0'))


if __name__ == '__main__':
    unittest.main()
