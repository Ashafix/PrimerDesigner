from flask import Flask
from flask_restful import Resource, Api
from flask_restful import reqparse
from flask import jsonify
import subprocess
import sys
import sqlite3
import time
import uuid
import concurrent.futures
import functools
import uuid
from FlaskJob import BlastJob, get_job_id, extract_hits_from_blast
import json
import flask


app = Flask(__name__)
api = Api(app)
jobs = dict()
executor = concurrent.futures.ThreadPoolExecutor(10)

parser = reqparse.RequestParser()
parser.add_argument('accession', action='append')
parser.add_argument('format', required=False, default='txt', choices=['txt', 'json'])


class RestBlastMinimal(Resource):

    def get(self, blast_id):
        if blast_id is None or len(blast_id) == 0:
            return flask.abort(404)
        if blast_id in jobs:
            return str(jobs[blast_id].stdout)

        return [str(j) for j in jobs.values()]

class RestBlast(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('blast_id', type=str, help='BLAST ID')
        args = parser.parse_args()
        return RestNucleotideMinimal().get(args['blast_id'])

    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('sequence', type=str, help='Sequence')
        args = parser.parse_args()

        job_id = get_job_id()
        cmd = ("INSERT INTO jobs VALUES ('{}', '{}', '{}', '{}', '{}')".format(job_id, args['sequence'], 0000, time.time(), 'submitted'))

        conn = sqlite3.connect('blast_jobs.db')
        c = conn.cursor()
        c.execute(cmd)
        conn.commit()
        conn.close()
        job = BlastJob()
        args['job_id'] = job_id
        job.future = executor.submit(functools.partial(job.run, parameters=args))
        jobs[job_id] = job

        return job_id


class RestNucleotide(Resource):
    def get(self):
        args = parser.parse_args()
        args['accession'] = args['accession'][0]
        if args['format'] == 'txt':
            return RestNucleotideMinimal().get(args['accession'])
        else:
            job = BlastJob()
            acc = job.get_accession(args['accession'])
            return jsonify({'response': acc})

    def post(self):
        args = parser.parse_args()
        job = BlastJob()
        resp = {}
        for acc in args['accession']:
            resp[acc] = job.get_accession(acc)
        if args['format'] == 'txt':
            return resp
        else:
            return jsonify(resp)


class RestNucleotideMinimal(Resource):
    def get(self, accession):
        job = BlastJob()
        acc = job.get_accession(accession)
        return flask.Response(acc, mimetype='txt')


class RestBlastHits(Resource):
    def get(self, blast_id):
        if blast_id is None or len(blast_id) == 0 or blast_id not in jobs:
            return flask.abort(404)
        if not jobs[blast_id].finished:
            return flask.abort(404)
        return extract_hits_from_blast(jobs[blast_id].stdout)

api.add_resource(RestBlast, '/blast/')
api.add_resource(RestBlastMinimal, '/blast/<blast_id>')
api.add_resource(RestBlastHits, '/blast/hits/<blast_id>')
api.add_resource(RestNucleotide, '/nucleotide/')
api.add_resource(RestNucleotideMinimal, '/nucleotide/<accession>')


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

