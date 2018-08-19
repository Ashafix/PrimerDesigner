import sys
import os
import io


class Amplicon:
    def __init__(self):
        self.accession = None
        self.forward = ''
        self.reverse = ''
        self.size = ''
        self.pos_start = None
        self.pos_end = None
        self.sequence = ''

    def __repr__(self):
        repr = ''
        repr += 'Accession: {}\n'.format(self.accession)
        repr += 'Forward  : {}\n'.format(self.forward)
        repr += 'Reverse  : {}\n'.format(self.reverse)
        repr += 'Size     : {} bp\n'.format(self.size)
        repr += 'Sequence : {}\n'.format(self.sequence)
        return repr

    @staticmethod
    def parse(header):
        amplicon = Amplicon()
        cells = header.strip().split(' ')
        amplicon.forward = cells[3]
        amplicon.reverse = cells[4]
        amplicon.size = int(cells[2][0:-2])
        cells = cells[0].split(':')
        amplicon.accession = cells[0][1:]
        cells = cells[1].split('+')
        amplicon.pos_start = int(cells[0])
        amplicon.pos_end = int(cells[1])
        assert (amplicon.size == abs(amplicon.pos_start - amplicon.pos_end) + 1)
        return amplicon


class IsPcrParser:
    def __init__(self):
        self.amplicons = []
        self.number_of_amplicons = None

    def parse(self, filename):
        if os.path.exists(filename):
            f = open(filename, 'r')
        else:
            f = filename.splitlines()
        for line in f:
            if line.startswith('>'):
                self.amplicons.append(Amplicon.parse(line))
            else:
                self.amplicons[-1].sequence += line.strip()

        self.number_of_amplicons = len(self.amplicons)
        if isinstance(f, io.TextIOWrapper):
            f.close()


if __name__ == '__main__':
    parser = IsPcrParser()
    parser.parse(sys.argv[1])
