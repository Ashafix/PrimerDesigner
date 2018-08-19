import sys
import os
import subprocess
import time
import primer3
import concurrent.futures
import functools
import yaml

from Bio import SeqIO
from PrimerDesigner.FlaskJob import BlastJob


class Primer:
    def __init__(self):
        self.seq = ''
        self.gc = None

    def __hash__(self):
        return hash(self.seq)

    def __eq__(self, other):
        return self.seq == other.seq

    def __repr__(self):
        repr = ''
        repr += 'Sequence: {}\n'.format(self.seq)
        repr += 'GC      : {}\n'.format(self.gc)
        return repr

    @staticmethod
    def parse_primer3(primer3_output, index=0, forward=True):
        if forward:
            pos = 'LEFT'
        else:
            pos = 'RIGHT'
        primer = Primer()
        primer.seq = primer3_output['PRIMER_{}_{}_SEQUENCE'.format(pos, index)]
        primer.gc = float(primer3_output['PRIMER_{}_{}_GC_PERCENT'.format(pos, index)])
        return primer


class PrimerPair:

    def __init__(self):
        self.forward = Primer()
        self.reverse = Primer()

    def __repr__(self):
        repr = ''
        repr += 'Forward: {}\n'.format(self.forward.seq)
        repr += 'Reverse: {}'.format(self.reverse.seq)
        return repr
    def __hash__(self):
        return hash(self.forward.seq + "|" + self.reverse.seq)

    def __eq__(self, other):
        return (self.forward.seq == other.forward.seq and self.reverse.seq == other.reverse.sequence) or (
                    self.forward.seq == other.reverse.seq and self.reverse.seq == other.forward.seq)

    @staticmethod
    def parse_primer3(primer3_output, index=0):
        pp = PrimerPair()
        pp.forward = Primer.parse_primer3(primer3_output, index=index, forward=True)
        pp.reverse = Primer.parse_primer3(primer3_output, index=index, forward=False)
        return pp


class GfServer:
    def __init__(self, port=12345, executable=None, file_2bit=None, file_fasta=None):
        self.port = port
        self.executable = executable
        self.process = None
        self.file_2bit = file_2bit
        self.file_fasta = file_fasta
        if executable is None:
            self.get_location()

        if file_fasta is not None and file_2bit is None:
            self.convert_fasta_to_2bit()

    def get_location(self):
        conf_filename = 'blast.conf'

        if os.path.isfile(conf_filename):
            with open(conf_filename, 'r') as f:
                self.executable = yaml.load(f)['gfserver']
        else:
            self.executable = os.environ.get('GFSERVER')
        if not os.path.isfile(self.executable):
            raise ValueError('gfServer executable not found in location: {}'.format(self.executable))
        return self.executable

    def start(self):
        self.process = subprocess.Popen(
            [self.executable, '-canStop', '-stepSize=5', 'start', 'localhost', str(self.port), self.file_2bit],
            stdout=subprocess.PIPE)
        return self.process

    def convert_fasta_to_2bit(self):
        """
        Converts a FASTA file to 2bit format
        :return: str, the full path of the 2bit file
        """
        exec_2bit = self.executable[0:self.executable.rfind('gfServer')] + 'faToTwoBit'
        if self.file_fasta.endswith('.fa'):
            self.file_2bit = self.file_fasta[0:self.file_fasta.rfind('.fa')] + '.2bit'
        else:
            self.file_2bit = self.file_fasta + '.2bit'
        p = subprocess.Popen([exec_2bit, self.file_fasta, self.file_2bit])
        p.communicate()
        return self.file_2bit

    @staticmethod
    def parse_response(response):
        if isinstance(response, bytes):
            response = response.decode()
        response = response.strip()

        return len(response.splitlines())

    def call(self, primer_pair, max_distance=1500, trials=100):
        sub = None
        while trials > 0 and (sub is None or sub.stderr != b''):
            sub = subprocess.run(
                [self.executable, 'pcr', 'localhost', str(self.port), primer_pair.forward.seq, primer_pair.reverse.seq,
                 str(max_distance)],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            trials -= 1
        if trials == 0:
            print('error: {}'.format(sub.stderr))
        return sub

    def stop(self):
        subprocess.run([self.executable, 'stop', 'localhost', str(self.port)])


def call_ispcr(primer_pair):
    with open('primerPair.txt', 'w') as f:
        f.write('{} {} {}'.format('bla', primer_pair.forward.seq, primer_pair.reverse.seq))
    x = subprocess.run(['./isPcr', 'top500.fa', 'primerPair.txt', 'stdout'], stdout=subprocess.PIPE)
    return x.stdout


def create_primers(record, number_of_primers=1000):
    return primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': record.id,
            'SEQUENCE_TEMPLATE': str(record.seq),
        },
        {
            'PRIMER_OPT_SIZE': 21,
            'PRIMER_PICK_INTERNAL_OLIGO': 1,
            'PRIMER_INTERNAL_MAX_SELF_END': 8,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 57.0,
            'PRIMER_MAX_TM': 65.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 80.0,
            'PRIMER_MAX_POLY_X': 100,
            'PRIMER_INTERNAL_MAX_POLY_X': 100,
            'PRIMER_SALT_MONOVALENT': 50.0,
            'PRIMER_DNA_CONC': 50.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_MAX_SELF_ANY': 12,
            'PRIMER_MAX_SELF_END': 8,
            'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            'PRIMER_PAIR_MAX_COMPL_END': 8,
            'PRIMER_PRODUCT_SIZE_RANGE': [[450, 1000]],
            'PRIMER_NUM_RETURN': number_of_primers
        })


def validate_primerpairs(primer_pairs, filename=None):
    gfserver = GfServer(file_fasta=filename)
    gfserver.start()
    validated = []
    for pp in primer_pairs:
        r = gfserver.call(pp)
        no_r = (GfServer.parse_response(r.stdout))
        if no_r == 1:
            validated.append(pp)

    gfserver.stop()
    return validated


def make_directories():
    os.makedirs(os.path.join(os.path.dirname(__file__), 'data', 'input'), exist_ok=True)


def write_sequence_to_file(sequence):

    filename = os.path.join('/tmp/', 'file.fa')
    with open(filename, 'w') as f:
        f.write(sequence)
    return filename


def design_primers(filename, number_of_primers, database='nt', primer_pairs_to_screen=3200):

    make_directories()
    # get target sequence

    if filename.startswith('>') and not os.path.isfile(filename):
        filename = write_sequence_to_file(filename)

    try:
        record = SeqIO.read(filename, 'fasta')
    except ValueError as e:
        raise e
    # run BLAST in the background
    blast = BlastJob(blast_db=database)
    blast.run(parameters={'sequence': record.format('fasta')})

    while not blast.finished:
        time.sleep(0.1)
    if blast.stderr is None or blast.stderr != '':
        raise RuntimeError('BLAST failed with error: {}'.format(blast.stderr))

    # get BLAST sequences
    executor = concurrent.futures.ThreadPoolExecutor(4)
    future_blast = executor.submit(functools.partial(blast.extract_hits_from_blast, blast.stdout))

    acc_hits = future_blast.result(timeout=120)
    # TODO default
    filename_hits = os.path.join(os.path.dirname(__file__), 'data', 'input', blast.get_job_id() + '.fa')
    with open(filename_hits, 'w') as f:
        f.write(blast.get_accession(acc_hits))

    primer_pairs = []
    valid_pairs = []
    primers = {}
    while len(valid_pairs) < number_of_primers:
        old_len = primers.get('PRIMER_LEFT_NUM_RETURNED', 0)
        future_primers = executor.submit(functools.partial(create_primers, record,
                                                           number_of_primers=primer_pairs_to_screen))
        primers = future_primers.result(timeout=120)
        for i in range(old_len, primers['PRIMER_LEFT_NUM_RETURNED']):
            pp = PrimerPair.parse_primer3(primers, index=i)
            primer_pairs.append(pp)
        valid_pairs = list(set(validate_primerpairs(primer_pairs, filename=filename_hits)))
        primer_pairs_to_screen = primer_pairs_to_screen * 2

    with open('optimal_pairs.txt', 'a') as f:
        f.write(str(primer_pairs_to_screen))
        f.write('\n')
    blast_outputs = []
    for v, valid in enumerate(valid_pairs):
        # run against all nucleotides
        for orientation in ('forward', 'reverse'):
            sequence = '>{}_{}\n{}'.format(orientation, v, valid.__getattribute__(orientation).seq)
            blast.run(parameters={'sequence': sequence})
            while not blast.finished:
                time.sleep(0.1)
            blast_outputs.extend(blast.extract_hits_from_blast(blast.stdout))
        # collect new sequences
        # print(blast_outputs, file=sys.stderr)
        # add new sequences to initial
        pass

    #acc_hits = blast.get_accessions_from_list(blast_outputs)
    #print(acc_hits, file=sys.stderr)

    return valid_pairs[0:number_of_primers]