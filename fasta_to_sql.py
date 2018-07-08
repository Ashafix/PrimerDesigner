import sys
import argparse


def fasta_to_csv(filename, output_location='stdout', output_header=True, output_sequence=False, delim="\t"):
	
	if output_sequence:
		output_format = ('{}' + delim) * 4 + '{}\n'
	else:
		output_format = ('{}' + delim) * 3 + '{}\n'
	if output_location == 'stdout':
		output_location = sys.stdout
	else:
		output_location = open(output_location, 'w')

	if output_header:
		if output_sequence:
			output = 'AccessionNumber{delim}Description{delim}Start{delim}Stop{delim}Sequence\n'.format(delim=delim)
		else:
			output = 'AccessionNumber{delim}Description{delim}Start{delim}Stop\n'.format(delim=delim)
	else:
		output = ''
	output_location.write(output)
	with open(filename, 'r') as f:
		for line_no, line in enumerate(f):
			if line.startswith('>'):
				if line_no != 0:
					output_location.write(output_format.format(accession, description, start, end, seq))
				start = line_no + 1
				header = line[1:].strip()
				accession = header[0:header.find(' ')]
				description = header[header.find(' ') + 1:]
				seq = ''
			else:
				if output_sequence:
					seq += line.strip()
				end = line_no + 1
		output_location.write(output_format.format(accession, description, start, end, seq))
	output_location.close()
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("filename", type=str, help="The FASTA file which should be translated")
	parser.add_argument("-o", "--output", type=str, help="The output location, default=stdout, i.e. printing", default='stdout')
	parser.add_argument("-s", "--sequence", type=bool, help="Output sequence?, default=False", default=False)
	parser.add_argument("-p", "--header", type=bool, help="Output header?, default=True", default=True)
	parser.add_argument("-d", "--delimiter", type=str, help="The delimiter used in the output", default="\t")
	args = parser.parse_args()
	print(args)
	fasta_to_csv(args.filename, output_location=args.output, output_header=args.header, output_sequence=args.sequence, delim=args.delimiter)
