import unittest
import os
import sys
from PrimerDesigner.Primer import design_primers


class DesignPrimers(unittest.TestCase):

    def setUp(self):
        try:
            os.remove(os.path.join(os.getcwd(), 'temp', 'tmp_blast.jobs.db'))
        except FileNotFoundError:
            pass

    def tearDown(self):
        try:
            os.remove(os.path.join(os.getcwd(), 'temp', 'tmp_blast.jobs.db'))
        except FileNotFoundError:
            pass

    @unittest.skip
    def test_design_primers(self):
        p = design_primers(os.path.join(os.getcwd(), 'data', 'random_sequence_0.fa'), 5,
                           database=os.path.join(os.getcwd(), 'data', 'random.fa'),
                           primer_pairs_to_screen=100)
        print(p)
        self.assertEqual(len(p), 5)


if __name__ == '__main__':
    unittest.main()
