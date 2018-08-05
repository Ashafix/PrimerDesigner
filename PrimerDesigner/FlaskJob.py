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
        with open('touch.me', 'w') as f:
            f.write(str(parameters))

        print(parameters, file=sys.stderr)
        filename_query = os.path.join(self.directory_db, 'blast_{}.fa'.format(parameters['job_id']))
        print(filename_query, file=sys.stderr)
        seq = parameters['sequence']
        with open(filename_query, 'w') as f:
            f.write(seq)
        call = [self.blast_executable,
                '-db', '{}/nt'.format(self.directory_db),
                '-outfmt', str(parameters['outfmt'])]

        if False:
            valid_blast_parameters = ('word_size', 'word_size', 'word_size', 'word_size', 'gapopen', 'gapextend',
                                      'gapopen', 'gapextend', 'reward', 'penalty', 'reward', 'penalty', 'reward',
                                      'penalty', 'strand', 'dust', 'filtering_db', 'window_masker_taxid',
                                      'window_masker_db', 'soft_masking', 'lcase_masking', 'db_soft_mask', 'db_hard_mask',
                                      'perc_identity', 'template_type', 'template_length', 'use_index', 'index_name',
                                      'xdrop_ungap', 'xdrop_gap', 'xdrop_gap_final', 'no_greedy', 'min_raw_gapped_score',
                                      'ungapped', 'window_size', 'outfmt', 'evalue')
            other_blast_parameters = ('sequence', 'forward', 'reverse', 'job_id', 'num_threads', 'job_id')
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
            h = hashlib.md5(('_'.join(call) + '_' + seq.split('\n', 1)[-1]).encode('utf-8')).hexdigest()
            conn = sqlite3.connect('blast_jobs.db')
            c = conn.cursor()
            cmd = "SELECT * FROM jobs WHERE id='{}'".format(h)
            c.execute(cmd)
            rows = c.fetchall()

        print(' '.join(call), file=sys.stderr)

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

        print(' '.join(call), file=sys.stderr)
        if self.stderr is None or len(self.stderr) > 0:
            self.error = True
        self.finished = True
        self.status = 'finished'
        if cache:
            conn.commit()
            conn.close()

        return parameters['job_id']

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
        print(parameters, file=sys.stderr)
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
                defaults = yaml.loads(f)
        else:
            defaults = {}
        defaults['short_sequence'] = int(defaults.get('short_sequence', 25))
        defaults['num_threads'] = int(defaults.get('num_threads', 6))
        defaults['outfmt'] = int(defaults.get('outfmt', 5))
        self.defaults = defaults
        return defaults