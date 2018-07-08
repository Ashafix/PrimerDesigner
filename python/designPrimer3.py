import sys
import subprocess
import primer3
from Bio import SeqIO
from isPcrParser import *


class Primer:
	def __init__(self):
		seq = ''
		gc = None
	
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
		forward = Primer()
		reverse = Primer()
	
	def __repr__(self):
		repr = ''
		repr += 'Forward: {}\n'.format(self.forward.seq)
		repr += 'Reverse: {}'.format(self.reverse.seq)
		return repr
	
	def __hash__(self):
		return hash(self.forward.seq + "|" + self.reverse.seq)
	
	def __eq__(self, other):
		return (self.forward.seq == other.forward.seq and self.reverse.seq == other.reverse.sequence) or (self.forward.seq == other.reverse.seq and self.reverse.seq == other.forward.seq)

	@staticmethod
	def parse_primer3(primer3_output, index=0):
		pp = PrimerPair()
		pp.forward = Primer.parse_primer3(primer3_output, index=index, forward=True)
		pp.reverse = Primer.parse_primer3(primer3_output, index=index, forward=False)
		return pp

class GfServer:
	def __init__(self, port=12345, executable='./gfServer'):
		self.port = port
		self.executable = executable
		self.process = None
		
	def start(self):
		self.process = subprocess.Popen([self.executable, '-canStop', '-stepSize=5', 'start', 'localhost', str(self.port), 'seq100.2bit'], stdout=subprocess.PIPE)
		return self.process
	
	@staticmethod
	def parse_response(response):
		if isinstance(response, bytes):
			response = response.decode()
		response = response.strip()
	
		return(len(response.splitlines()))
	
	def call(self, primer_pair, max_distance=1500, trials=100):
		sub = None
		while trials > 0 and (sub is None or sub.stderr != b''):
			sub = subprocess.run([self.executable, 'pcr', 'localhost', str(self.port), primer_pair.forward.seq, primer_pair.reverse.seq, str(max_distance)], 
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
        'PRIMER_PRODUCT_SIZE_RANGE': [[450,1000]],
        'PRIMER_NUM_RETURN': number_of_primers
    })

def validate_primerpairs(primer_pairs):
	gfserver = GfServer()
	gfserver.start()
	validated = list()
	for pp in primer_pairs:
		r = gfserver.call(pp)
		no_r = (GfServer.parse_response(r.stdout))
		if no_r == 1:
			validated.append(pp)
			
	gfserver.stop()
	return validated

def main(filename, number_of_primers):
	#get target sequence
	record = SeqIO.read(filename, 'fasta')

	#run BLAST in the background

	#get BLAST sequences
	primers = create_primers(record, number_of_primers=number_of_primers)

	primer_pairs = list()
	for i in range(number_of_primers):
		pp = PrimerPair.parse_primer3(primers, index=i)
		primer_pairs.append(pp)

	valid_pairs = validate_primerpairs(primer_pairs)

	for valid in valid_pairs:
		#run against all nucleotides
		
		#collect new sequences
		
		#add new sequences to initial
		print(valid)
		pass

if __name__ == '__main__':
	main(sys.argv[1], int(sys.argv[2]))
