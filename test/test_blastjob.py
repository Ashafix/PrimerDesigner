import unittest
import os
import sys
import sqlite3
from PrimerDesigner import FlaskJob


class Blast(unittest.TestCase):

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

    def test_basic_blast_positive(self):
        job = FlaskJob.BlastJob(conf_file=os.path.join(os.getcwd(), 'data', 'blast.conf'),
                                blast_db=os.path.join(os.getcwd(), 'data', 'random.fa'))
        job.set_database(os.path.join(os.getcwd(), 'data', 'random.fa'))
        job.result_db = os.path.join(os.getcwd(), 'temp', 'tmp_blast.jobs.db')
        self.assertIsNotNone(job.blast_db)
        self.assertIsNotNone(job.blast_executable)
        self.assertEqual(job.stdout, '')
        self.assertEqual(job.stderr, '')
        self.assertNotEqual(len(job.defaults), 0)

    def test_basic_blast_negative(self):
        self.assertRaises(ValueError, FlaskJob.BlastJob, 'I should not exist.txt')

    def test_basic_blast_run(self):
        job = FlaskJob.BlastJob(conf_file=os.path.join(os.getcwd(), 'data', 'blast.conf'),
                                blast_db=os.path.join(os.getcwd(), 'data', 'random.fa'))
        job.set_database(os.path.join(os.getcwd(), 'data', 'random.fa'))
        job.result_db = os.path.join(os.getcwd(), 'temp', 'tmp_blast.jobs.db')
        job.directory_database = os.path.join(os.getcwd(), 'test')
        with open(os.path.join(os.getcwd(), 'data', 'random_sequence_0.fa'), 'r') as f:
            seq = f.read()
        parameters = {'sequence': seq,
                      'num_threads': 1}
        job.run(parameters)
        self.assertNotEqual(job.stdout, '')
        self.assertEqual(job.stderr, '')

    def test_basic_blast_run_from_filename(self):
        job = FlaskJob.BlastJob(conf_file=os.path.join(os.getcwd(), 'data', 'blast.conf'),
                                blast_db=os.path.join(os.getcwd(), 'data', 'random.fa'))
        job.set_database(os.path.join(os.getcwd(), 'data', 'random.fa'))
        job.result_db = os.path.join(os.getcwd(), 'temp', 'tmp_blast.jobs.db')
        job.directory_database = os.path.join(os.getcwd(), 'data')
        parameters = {'sequence': os.path.join(os.getcwd(), 'data', 'random_sequence_0.fa'),
                      'num_threads': 1}
        job.run(parameters, query_is_file=True)
        self.assertEqual(job.stderr, '')
        self.assertNotEqual(job.stdout, '')

    def test_basic_blast_run_negative(self):
        job = FlaskJob.BlastJob(conf_file=os.path.join(os.getcwd(), 'data', 'blast.conf'),
                                blast_db=os.path.join(os.getcwd(), 'data', 'random.fa'))
        job.set_database(os.path.join(os.getcwd(), 'data', 'random.fa'))
        job.result_db = os.path.join(os.getcwd(), 'temp', 'tmp_blast.jobs.db')
        job.directory_database = os.path.join(os.getcwd(), 'temp')
        with open(os.path.join(os.getcwd(), 'data', 'random_sequence_0.fa'), 'r') as f:
            seq = f.read()
        parameters = {'sequence': seq,
                      'num_threads': 1024}
        job.run(parameters)
        self.assertEqual(job.stdout, '')
        self.assertNotEqual(job.stderr, '')

    def test_blast_cache(self):
        job = FlaskJob.BlastJob(conf_file=os.path.join(os.getcwd(), 'data', 'blast.conf'),
                                blast_db=os.path.join(os.getcwd(), 'data', 'random.fa'))
        job.result_db = os.path.join(os.getcwd(), 'temp', 'tmp_blast.jobs.db')
        job.directory_database = os.path.join(os.getcwd(), 'data')
        with open(os.path.join(os.getcwd(), 'data', 'random_sequence_0.fa'), 'r') as f:
            seq = f.read()
        parameters = {'sequence': seq,
                      'num_threads': 1}
        job.run(parameters)
        conn = sqlite3.connect(job.result_db)
        c = conn.cursor()
        c.execute('SELECT COUNT (id) FROM jobs')
        r0 = c.fetchall()
        parameters = {'sequence': seq + 'A',
                      'num_threads': 1}
        job.run(parameters)

        c.execute('SELECT COUNT (id) FROM jobs')
        r1 = c.fetchall()
        self.assertLess(r0[0], r1[0])

        parameters = {'sequence': seq + 'AT',
                      'num_threads': 1}
        job.run(parameters, cache=False)

        c.execute('SELECT COUNT (id) FROM jobs')
        r2 = c.fetchall()
        self.assertEqual(r2[0], r1[0])




if __name__ == '__main__':
    unittest.main()
