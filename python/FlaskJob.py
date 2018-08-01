import subprocess
import sys
import io
import os
import yaml
import uuid
from Bio.Blast import NCBIXML


class FlaskJob:
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


class BlastJob(FlaskJob):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.blast_executable = ''
        self.directory_db = ''
        self.directory_tmp = '/tmp/'
        self.directory_query = ''
        self.get_locations()


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

    def run(self, parameters=None):
        if parameters is None:
            parameters = {}
        job_id = parameters.get('job_id', 'job_id')
        num_threads = parameters.get('num_threads', 4)
        output_format = parameters.get('outfmt', 5)
        filename_query = '{}/{}.fa'.format(self.directory_db, job_id)
        seq = seq_from_fasta(parameters['sequence'])

        with open(filename_query, 'w') as f:
            f.write(seq)
        call = [self.blast_executable,
                '-db', '{}/nt'.format(self.directory_db),
                '-query', filename_query,
                '-num_threads', str(num_threads),
                '-outfmt', str(output_format)]
        proc = subprocess.Popen(call,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        self.status = 'running'
        self.stdout, self.stderr = proc.communicate()
        self.stdout = self.stdout.decode('utf-8')
        self.stderr = self.stderr.decode('utf-8')

        if self.stderr is None or len(self.stderr) > 0:
            self.error = True
        self.finished = True
        self.status = 'finished'

    def get_accession(self, accession):

        filename = None
        while filename is None or os.path.exists(filename):
            filename = os.path.join(self.directory_tmp, "blastdbcmd_{}".format(get_job_id()))
        accession = '\n'.join(accession.split(';'))
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
            return {'error': stderr}
        return stdout


def get_job_id():
    return uuid.uuid4().hex[0:8]


def seq_from_fasta(seq):
    return seq
    seq = seq.strip()
    lines = seq.splitlines()
    if lines[0].startswith('>'):
        start = 1
    else:
        start = 0
    return ''.join(seq[start:])


def extract_hits_from_blast(blast_xml):
    if isinstance(blast_xml, str):
        blast_xml = io.StringIO(blast_xml)
    blast_records = NCBIXML.read(blast_xml)
    hits = []
    for alignment in blast_records.alignments:
        hits.append(alignment.accession)
    return hits