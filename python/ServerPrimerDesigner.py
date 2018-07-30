from flask import Flask
from flask_restful import Resource, Api
from flask_restful import reqparse
import subprocess
import sys
import sqlite3
import time
import uuid
import concurrent.futures
import functools
import uuid
from FlaskJob import BlastJob

app = Flask(__name__)
api = Api(app)
jobs = dict()
executor = concurrent.futures.ThreadPoolExecutor(10)


class RestBlast(Resource):
    def get(self, blast_id):
        if blast_id is None:
            return 'noooone'
        if blast_id in jobs:
            return str(jobs[blast_id].stdout)

        return [str(j) for j in jobs.values()]

    def post(self, blast_id):
        parser = reqparse.RequestParser()
        parser.add_argument('sequence', type=str, help='Sequence')
        args = parser.parse_args()
        #conn = sqlite3.connect('blast_jobs.db')
        #c = conn.cursor()

        job_id = uuid.uuid4().hex
        cmd = ("INSERT INTO jobs VALUES ('{}', '{}', '{}', '{}', '{}')".format(job_id, args['sequence'], 0000, time.time(), 'submitted'))
        #with open('job_id.fa', 'w') as f:
        #    f.write('>{}\n{}'.format('bla', args['sequence']))
        #with open('f.txt', 'w') as f:
        #    f.write(cmd)
            #c.execute(cmd)
            
            #call = ['/media/ashafix/WesternDigital/ncbi-blast-2.7.1+/bin/blastn', '-db', '/media/ashafix/WesternDigital/db/nt', '-query', '/media/ashafix/WesternDigital/job_id.fa', '>>', job_id + 'txt']
            #f.write('\n'.join(call))
            #subprocess.run(' '.join(call), shell=True)
        job_id = get_job_id()
        jobs[job_id] = BlastJob()
        print('starting', file=sys.stderr)
        jobs[job_id].future = executor.submit(jobs[job_id].run)
        print('started', file=sys.stderr)

        return job_id
        #conn.commit()
        #conn.close()

class RestNucleotide(Resource):
    def get(self):
        return {'id': 'values'}

def get_job_id():
    return uuid.uuid4().hex[0:8]


api.add_resource(RestBlast, '/blast/<blast_id>')
api.add_resource(RestNucleotide, '/nucleotide')


if __name__ == '__main__':
    conn = sqlite3.connect('blast_jobs.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE jobs
                 (id text, sequence text, parameters text, date text, status text)''')
    except:
        pass
    conn.commit()
    conn.close()
    app.run(host='0.0.0.0', debug=True)

