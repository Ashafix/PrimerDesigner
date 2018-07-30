import subprocess
import sys
import os
import yaml


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
        return 'FlaskJob(status={}, error={}, finished={}, future={}'.format(self.status,
                                                                             self.error,
                                                                             self.finished,
                                                                             self.future)


class BlastJob(FlaskJob):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.blast_executable = ''
        self.db_directory = ''
        self.query_directory = ''
        self.get_locations()

    def get_locations(self):

        conf_filename = 'blast.conf'

        if os.path.isfile(conf_filename):
            with open(conf_filename, 'r') as f:
                config = yaml.load(f)
        else:
            config = {'blast': os.environ.get('BLAST_EXECUTABLE'),
                      'blast_db': os.environ.get('BLAST_DATABASE')}
        print(config, file=sys.stderr)
        self.blast_executable = config.get('blast')
        self.db_directory = config.get('blast_db')
        if self.blast_executable is None or self.db_directory is None:
            raise ValueError('Blast executable and database dir need to be present '
                             'in config file or environment variable')
        self.query_directory = os.path.join(self.db_directory, '..', 'query')
        if not os.path.isdir(self.query_directory):
            os.makedirs(self.query_directory)
        return self.blast_executable, self.db_directory

    def run(self, parameters=None):
        if parameters is None:
            parameters = {}
        job_id = parameters.get('job_id', 'job_id')
        num_threads = parameters.get('num_threads', 4)
        output_format = parameters.get('outfmt', 5)
        call = [self.blast_executable,
                '-db', '{}/nt'.format(self.db_directory),
                '-query', '{}/{}.fa'.format(self.db_directory, job_id),
                '-num_threads', str(num_threads),
                '-outfmt', str(output_format)]
        proc = subprocess.Popen(call,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        self.status = 'running'
        self.stdout, self.stderr = proc.communicate()
        self.stdout = self.stdout.decode('utf-8')
        self.stderr = self.stderr.decode('utf-8')

        print(type(self.stdout), file=sys.stderr)
        if self.stderr is None or len(self.stderr) > 0:
            self.error = True
        self.finished = True
        self.status = 'finished'
