import unittest
import PrimerDesigner

class Job(unittest.TestCase):

    def test_basic_function(self):
        job = PrimerDesigner.Job.Job()
        self.assertNotEqual(str(job), '')
        self.assertNotEqual(job.__repr__(), '')

    def test_basic_attributes(self):
        job = PrimerDesigner.Job.Job()
        self.assertEqual(job.stdout, '')
        self.assertEqual(job.stderr, '')


if __name__ == '__main__':
    unittest.main()
