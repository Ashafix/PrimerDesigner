import sys
import io
import os
import yaml
import uuid
import subprocess
import hashlib
import time
import sqlite3
import multiprocessing
from Bio.Blast import NCBIXML
from PrimerDesigner.tools import tools
#from . import tools


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
    def __init__(self, conf_file='blast.conf', result_db='blast_jobs.db', blast_db=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.blast_executable = ''
        self.directory_database = ''
        self.directory_tmp = '/tmp/'
        self.directory_query = ''
        self.conf_file = conf_file
        self.blast_db = None
        self.set_database(blast_db)
        self.get_locations()
        self.defaults = self._get_defaults()
        self.result_db = result_db
        self.run_hash = None

    def get_locations(self):
        if os.path.isfile(self.conf_file):
            with open('/tmp/tmp.txt', 'a') as f:
                f.write('\nfound conf file')

            with open(self.conf_file, 'r') as f:
                config = yaml.load(f)
        else:
            config = {'blast': os.environ.get('BLAST_EXECUTABLE'),
                      'blast_dir': os.environ.get('BLAST_DIRECTORY')}
        self.blast_executable = config.get('blast')
        self.directory_database = config.get('blast_dir')
        with open('/tmp/tmp.txt', 'a') as f:
            f.write('\n')
            f.write(os.getcwd())
            f.write(str(self.blast_executable))
            f.write('\n')
            f.write(str(self.directory_tmp))
            f.write('\n')
            f.write(str(self.blast_db))
        if self.blast_executable is None or (self.directory_tmp is None and not os.path.isfile(self.blast_db)):
            raise ValueError('Blast executable and database dir need to be present '
                             'in config file or environment variable')
        if self.directory_query is None or len(self.directory_query) == 0:
            self.directory_query = os.path.join(os.path.dirname(self.directory_database), '..', 'query')

        os.makedirs(self.directory_database, exist_ok=True)
        os.makedirs(self.directory_query, exist_ok=True)
        os.makedirs(self.directory_tmp, exist_ok=True)

        return self.blast_executable, self.blast_db

    def set_database(self, database):
        if database is None:
            return None
        elif not isinstance(database, str):
            raise TypeError('database needs to be a str, got {}'.format(type(database)))
        if os.path.isfile(database):
            self.blast_db = database
        elif os.path.join(self.directory_database, database):
            self.blast_db = os.path.join(self.directory_database, database)
        else:
            raise ValueError('could not find database or in directory {}'.format(self.directory_database))
        return self.blast_db

    def run(self, parameters, cache=True, query_is_file=False, delete_query_file=False):
        parameters = self._clean_parameters(parameters, query_is_file=query_is_file)

        if query_is_file:
            filename_query = parameters['sequence']
            with open(filename_query, 'r') as f:
                seq = f.read()
        else:
            filename_query = os.path.join(self.directory_query, 'blast_{}.fa'.format(parameters['job_id']))
            seq = parameters['sequence']
            with open(filename_query, 'w') as f:
                f.write(seq)
        call = [self.blast_executable,
                '-db', '{}'.format(self.blast_db),
                '-outfmt', str(parameters['outfmt'])]

        valid_blast_parameters = ('word_size', 'word_size', 'word_size', 'word_size', 'gapopen', 'gapextend',
                                  'gapopen', 'gapextend', 'reward', 'penalty', 'reward', 'penalty', 'reward',
                                  'penalty', 'strand', 'dust', 'filtering_db', 'window_masker_taxid',
                                  'window_masker_db', 'soft_masking', 'lcase_masking', 'db_soft_mask', 'db_hard_mask',
                                  'perc_identity', 'template_type', 'template_length', 'use_index', 'index_name',
                                  'xdrop_ungap', 'xdrop_gap', 'xdrop_gap_final', 'no_greedy', 'min_raw_gapped_score',
                                  'ungapped', 'window_size', 'evalue')
        other_blast_parameters = ('sequence', 'forward', 'reverse', 'job_id', 'num_threads', 'job_id', 'outfmt')
        for param_k, param_v in parameters.items():
            if param_k in valid_blast_parameters:
                call.append('-{}'.format(param_k))
                call.append(param_v)
            elif param_k not in other_blast_parameters:
                print('encountered invalid parameters: {}'.format(param_k), file=sys.stderr)

        if len(seq.split('\n', 1)[-1]) < self.defaults['short_sequence']:
            call.append('-task')
            call.append('blastn-short')

        if cache:
            self.run_hash = hashlib.md5(('_'.join(call) + '_' + seq.split('\n', 1)[-1]).encode('utf-8')).hexdigest()
            rows = self.get_cached_results()
        else:
            rows = None

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
        if self.stderr is None or len(self.stderr) > 0:
            self.error = True
        self.finished = True
        self.status = 'finished'
        if cache:
            self.write_cached_results(seq)
        if query_is_file and delete_query_file:
            try:
                os.remove(filename_query)
            except (FileNotFoundError, PermissionError):
                pass
        return parameters['job_id']

    def get_cached_results(self):
        conn = sqlite3.connect(self.result_db)
        c = conn.cursor()
        cmd = "SELECT * FROM jobs WHERE id='{}'".format(self.run_hash)
        try:
            c.execute(cmd)
        except sqlite3.OperationalError as e:
            if 'no such table' in str(e):
                tools.create_empty_database(conn)
            else:
                raise e
            c.execute(cmd)
        rows = c.fetchall()
        conn.close()
        return rows

    def write_cached_results(self, seq):
        conn = sqlite3.connect(self.result_db)
        c = conn.cursor()
        cmd = ("INSERT INTO jobs VALUES ('{}', '{}', '{}', '{}', '{}', '{}', '{}')".format(self.run_hash,
                                                                                           seq,
                                                                                           0000,
                                                                                           time.time(),
                                                                                           'finished',
                                                                                           self.stdout.replace("'", "''"),
                                                                                           self.stderr.replace("'", "''")))
        c.execute(cmd)
        conn.commit()
        conn.close()

    def set_arguments_for_primer_blast(self, parameters=None, cache=True):
        if parameters is None:
            parameters = {}

        # based on the idea described here:
        # https://eu.idtdna.com/pages/education/decoded/article/tips-for-using-blast-to-locate-pcr-primers
        parameters['word_size'] = 7
        parameters['evalue'] = 1000
        # TODO check if work as expected
        parameters['dust'] = 'no'

        parameters['sequence'] = parameters['forward'] + 'N' * 10 + parameters['reverse']
        return parameters
    
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
        #TODO: make this replace prettier
        call = [self.blast_executable.replace('blastn', 'blastdbcmd'),
                '-db', self.blast_db,
                '-entry_batch', filename]
        proc = subprocess.Popen(call,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                universal_newlines=True)
        stdout, stderr = proc.communicate()
        if (stderr is not None and len(stderr.strip()) > 0) or stdout.startswith('Error:'):
            error = '\n'.join([stderr, stdout])
            raise ValueError('blastdbcmd failed with error: {}'.format(error))
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

    def _clean_parameters(self, parameters, query_is_file=False):
        seq = parameters.get('sequence')
        if seq is None:
            raise ValueError('No sequence provided')
        job_id = parameters.get('job_id', self.get_job_id())
        seq = seq.strip()
        if len(seq) == 0:
            raise ValueError('No sequence provided')

        if not seq.startswith('>') and not query_is_file:
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

        outfmt = str(parameters.get('outfmt', self.defaults['outfmt'])).lower()
        valid_outfmt = ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11')
        outfmt_parts = outfmt.split(' ')
        if outfmt_parts[0] not in valid_outfmt:
            raise ValueError('outfmt needs to be in {}'.format(valid_outfmt))
        if len(outfmt_parts) > 1:
            valid_fmtext = ('qseqid', 'qgi', 'qacc', 'sseqid', 'sallseqid', 'sgi', 'sallgi', 'sacc', 'sallacc',
                            'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'evalue', 'bitscore', 'score', 'length',
                            'pident', 'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'frames', 'qframe',
                            'sframe', 'btop', 'staxids', 'sscinames', 'scomnames', 'sblastnames', 'vsskingdoms',
                            'stitle', 'salltitles', 'sstrand', 'qcovs', 'qcovhsp', 'qcovus')
            for p in outfmt_parts[1:]:
                if p not in valid_fmtext:
                    raise ValueError('found invalid format extension: {}'.format(p))
        parameters['outfmt'] = outfmt

        return parameters

    def _get_defaults(self):
        filename = 'blast_defaults'
        if os.path.isfile(filename):
            with open(filename, 'r') as f:
                defaults = yaml.load(f)
        else:
            defaults = {}
        defaults['short_sequence'] = int(defaults.get('short_sequence', 25))
        defaults['num_threads'] = int(defaults.get('num_threads', min(6, multiprocessing.cpu_count())))
        defaults['outfmt'] = int(defaults.get('outfmt', 5))
        self.defaults = defaults
        return defaults
