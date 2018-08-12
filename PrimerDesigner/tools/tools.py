import os
import random
import ftplib
import tarfile
import sqlite3


def download_from_ftp(url, direc='.'):
    """
    Downloads a file from a FTP server and stores it in a directory
    :param url: str, the url of the file to download
    :param direc: str, the target directory
    :return: int, size of the downloaded file
    """

    prefix = 'ftp://'
    if url.startswith(prefix):
        url = url[len(prefix) + 1:]

    cells = url.split('/')
    if len(cells) < 2:
        raise ValueError()

    ftp = ftplib.FTP(cells[0])
    ftp.login()
    if len(cells) > 2:
        ftp.cwd('/'.join(cells[1:-1]))
    filename = cells[-1]
    size = ftp.size(filename)
    if cells[-1].endswith('.tar.gz'):
        with tarfile.open(filename, 'r:gz') as tar:
            tar.extractall(direc)
    return size


def random_sequence(seq_len=1000, seq_type='nuc'):
    seq_type = seq_type.lower()
    if seq_type in ('nuc', 'nucleotide'):
        seq_char = 'ATGC'
    elif seq_type == 'protein':
        seq_char = 'ACDEFGHIKLMNPQRSTVWY'
    else:
        raise ValueError("seq_type must either be 'nuc' or 'protein'")

    return ''.join(random.choice(seq_char) for _ in range(seq_len))


def create_random_fasta(num_seq=10, seq_type='nuc'):
    seq = []
    for i in range(num_seq):
        seq.append('>random_sequence_{}\n{}'.format(i, random_sequence(seq_type=seq_type)))
    return '\n'.join(seq)


def create_empty_database(conn=None):
    close = False
    if conn is None:
        conn = sqlite3.connect('blast_jobs.db')
        close = True
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE jobs
                         (id text, sequence text, parameters text, date text, status text, stdout text, stderr text)''')
    except sqlite3.OperationalError as e:
        if 'table jobs already exists' not in str(e):
            raise e

    conn.commit()
    if close:
        conn.close()


def make_dirs(dirs=None):
    """
    Creates a list of directories needed to store files for PrimerDesigner
    :param dirs: list, list of directories to create
    :return: None
    """

    if dirs is None:
        dirs = ['bin', 'data', 'db']
    for direc in dirs:
        direc = os.path.join(os.getcwd(), '..', direc)
        os.makedirs(direc, exist_ok=True)
