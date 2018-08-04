import sys
import io
import os
import yaml
import uuid
import subprocess
import hashlib
import time
import sqlite3
from Bio.Blast import NCBIXML


class Job:
    def __init__(self, status='submitted', error=False, finished=False, future=None, executor=None):
        self.status = status
        self.error = error
        self.finished = finished
        self.future = future
        self.stdout = ''
        self.stderr = ''
        self.executor = executor

    def __str__(self):
        return str({'status': self.status, 'error': self.error, 'finished': self.finished,
                    'output': self.stdout, 'stderr': self.stderr})

    def __repr__(self):
        return '{}(status={}, error={}, finished={}, future={}'.format(self.__class__,
                                                                       self.status,
                                                                       self.error,
                                                                       self.finished,
                                                                       self.future)


class BlastJob(Job):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.blast_executable = ''
        self.directory_db = ''
        self.directory_tmp = '/tmp/'
        self.directory_query = ''
        self.get_locations()
        self.defaults = self._get_defaults()

    def get_locations(self):

        conf_filename = 'blast.conf'

        if os.path.isfile(conf_filename):
            with open(conf_filename, 'r') as f:
                config = yaml.load(f)
        else:
            config = {'blast': os.environ.get('BLAST_EXECUTABLE'),
                      'blast_db': os.environ.get('BLAST_DATABASE')}
        self.blast_executable = config.get('blast')
        self.directory_db = config.get('blast_db')
        if self.blast_executable is None or self.directory_db is None:
            raise ValueError('Blast executable and database dir need to be present '
                             'in config file or environment variable')
        if self.directory_query is None or len(self.directory_query) == 0:
            self.directory_query = os.path.join(self.directory_db, '..', 'query')

        if not os.path.isdir(self.directory_query):
            os.makedirs(self.directory_query)
        return self.blast_executable, self.directory_db

    def run(self, parameters=None, cache=True):
        parameters = self._clean_parameters(parameters)
        filename_query = os.path.join(self.directory_db, 'blast_{}.fa'.format(parameters['job_id']))
        seq = parameters['sequence']
        with open(filename_query, 'w') as f:
            f.write(seq)
        call = [self.blast_executable,
                '-db', '{}/nt'.format(self.directory_db),
                '-outfmt', str(parameters['outfmt'])]
        if len(seq.split('\n', 1)[-1]) < self.defaults['short_sequence']:
            call.append('-task')
            call.append('blastn-short')

        if cache:
            h = hashlib.md5(('_'.join(call) + '_' + seq.split('\n', 1)[-1]).encode('utf-8')).hexdigest()
            conn = sqlite3.connect('blast_jobs.db')
            c = conn.cursor()
            cmd = "SELECT * FROM jobs WHERE id='{}'".format(h)
            c.execute(cmd)
            rows = c.fetchall()

        if cache and len(rows) > 0:
            self.stdout = rows[0][5]
            self.stderr = rows[0][6]
        else:
            call.append('-query')
            call.append(filename_query)
            call.append('-num_threads')
            call.append(str(parameters['num_threads']))
            proc = subprocess.Popen(call,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

            self.status = 'running'
            self.stdout, self.stderr = proc.communicate()
            self.stdout, self.stderr = self.stdout.decode('utf-8'), self.stderr.decode('utf-8')
            if cache:
                cmd = ("INSERT INTO jobs VALUES ('{}', '{}', '{}', '{}', '{}', '{}', '{}')".format(h,
                                                                                                   seq,
                                                                                                   0000,
                                                                                                   time.time(),
                                                                                                   'finished',
                                                                                                   self.stdout,
                                                                                                   self.stderr))
                c.execute(cmd)

        if self.stderr is None or len(self.stderr) > 0:
            self.error = True
        self.finished = True
        self.status = 'finished'
        if cache:
            conn.commit()
            conn.close()

        return parameters['job_id']

    def get_accession(self, accession):

        filename = None
        while filename is None or os.path.exists(filename):
            filename = os.path.join(self.directory_tmp, "blastdbcmd_{}".format(BlastJob.get_job_id()))
        if isinstance(accession, str):
            accession = accession.replace(';', '\n')
        elif isinstance(accession, list):
            accession = '\n'.join(accession)
        else:
            raise ValueError('accession must be either str or list')

        with open(filename, 'w') as f:
            f.write(accession)
        call = [self.blast_executable.replace('blastn', 'blastdbcmd'),
                '-db', '{}/nt'.format(self.directory_db),
                '-entry_batch', filename]
        proc = subprocess.Popen(call,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True)
        stdout, stderr = proc.communicate()
        if stderr is None or len(stderr) > 0:
            print(stderr, file=sys.stderr)
            return {'error': stderr}
        return stdout

    def get_accessions_from_list(self, accessions):
        acc = list()
        for accession in accessions:
            acc.append(self.get_accession(accession))

        acc = list(set(acc))
        return acc

    @staticmethod
    def extract_hits_from_blast(blast_xml):
        if isinstance(blast_xml, str):
            blast_xml = io.StringIO(blast_xml)
        blast_records = NCBIXML.read(blast_xml)
        hits = []
        for alignment in blast_records.alignments:
            hits.append(alignment.accession)
        return hits

    @staticmethod
    def get_job_id():
        return uuid.uuid4().hex[0:8]

    def _clean_parameters(self, parameters):
        seq = parameters.get('sequence')
        if seq is None:
            raise ValueError('No sequence provided')
        job_id = parameters.get('job_id', self.get_job_id())
        seq = seq.strip()
        if len(seq) == 0:
            raise ValueError('No sequence provided')

        if not seq.startswith('>'):
            seq = '>{}\n{}'.format(job_id, seq)
        
        parameters['sequence'] = seq
        parameters['job_id'] = job_id

        num_threads = parameters.get('num_threads', self.defaults['num_threads'])
        if not isinstance(num_threads, int):
            try:
                num_threads = int(num_threads)
            except ValueError:
                raise ValueError('num_threads must be an integer')
        if num_threads < 1:
            raise ValueError('num_threads needs to be 1 or higher')
        parameters['num_threads'] = num_threads

        outfmt = parameters.get('outfmt', self.defaults['outfmt'])
        if not isinstance(outfmt, int):
            try:
                num_threads = int(num_threads)
            except ValueError:
                raise ValueError('outfmt must be an integer')
        valid_outfmt = (1, 2, 3, 4, 5, 6)
        if outfmt not in valid_outfmt:
            raise ValueError('outfmt needs to be in {}'.format(valid_outfmt))
        parameters['outfmt'] = outfmt

        return parameters

    def _get_defaults(self):
        filename = 'blast_defaults'
        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                defaults = yaml.loads(f)
        else:
            defaults = {}
        defaults['short_sequence'] = int(defaults.get('short_sequence', 25))
        defaults['num_threads'] = int(defaults.get('num_threads', 6))
        defaults['outfmt'] = int(defaults.get('outfmt', 5))
        self.defaults = defaults
        return defaults