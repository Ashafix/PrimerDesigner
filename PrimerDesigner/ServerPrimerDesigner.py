import flask
from flask import Flask
from flask_restful import Resource, Api
from flask_restful import reqparse
from flask import jsonify
import sys
import os
import sqlite3
import time
import concurrent.futures
import functools
from PrimerDesigner.FlaskJob import BlastJob
import PrimerDesigner.design_primers
from PrimerDesigner.tools import tools

app = Flask(__name__)
api = Api(app)
jobs = {}
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

        job = BlastJob()
        job_id = job.get_job_id()
        cmd = ("INSERT INTO jobs VALUES ('{}', '{}', '{}', '{}', '{}', '{}', '{}')".format(job_id,
                                                                                           args['sequence'],
                                                                                           0000,
                                                                                           time.time(),
                                                                                           'submitted',
                                                                                           '',
                                                                                           ''))

        conn = sqlite3.connect(job.result_db)
        c = conn.cursor()
        c.execute(cmd)
        conn.commit()
        conn.close()
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
        return BlastJob.extract_hits_from_blast(jobs[blast_id].stdout)

class RestDesignPrimers(Resource):
    def get(self):
        pass

    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('sequence', type=str, help='Sequence')
        parser.add_argument('number_of_pairs', type=int, required=False, default=5, help='Number of pairs to design')

        args = parser.parse_args()
        # TODO default
        filename = os.path.join(os.getcwd(), '..', 'input', BlastJob.get_job_id() + '.fa')
        resp = {'error': False, 'primers': [], 'message': ''}
        with open(filename, 'w') as f:
            f.write(args['sequence'])
        try:
            primers = design_primers(filename, args['number_of_pairs'])
        except ValueError as e:
            resp['error'] = True
            resp['message'] = str(e)
            return jsonify(resp)
        for p, primer in enumerate(primers):
            primers[p] = str(primer)
        resp['primers'] = primers
        if len(primers) == 0:
            resp['message'] = 'Failed to design primers for the target sequence'
        else:
            resp['message'] = 'Successfully designed primers'
        return jsonify(resp)

class RestBlastPrimers(Resource):
    def get(self):
        RestBlast().get()

    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('forward', type=str, help='The sequence of the forward primer')
        parser.add_argument('reverse', type=str, help='The sequence of the reverse primer')

        args = parser.parse_args()
        args['sequence'] = args['forward'] + '_' + args['reverse']
        job = BlastJob()
        job_id = job.get_job_id()
        cmd = (
            "INSERT INTO jobs VALUES ('{}', '{}', '{}', '{}', '{}', '{}', '{}')".format(job_id,
                                                                                        args['sequence'],
                                                                                        0000,
                                                                                        time.time(),
                                                                                        'submitted',
                                                                                        '',
                                                                                        '',
                                                                                        ''))

        conn = sqlite3.connect(job.result_db)
        c = conn.cursor()
        c.execute(cmd)
        conn.commit()
        conn.close()
        args['job_id'] = job_id
        print(args, file=sys.stderr)
        args = job.set_arguments_for_primer_blast(args)
        print(str(args) + '#' * 20, file=sys.stderr)
        job.future = executor.submit(functools.partial(job.run, parameters=args))
        jobs[job_id] = job

        return job_id


class RestShutdown(Resource):

    def get(self):
        func = flask.request.environ.get('werkzeug.server.shutdown')
        if func is None:
            raise RuntimeError('Not running with the Werkzeug Server')
        func()

    def post(self):
        self.get()


api.add_resource(RestBlast, '/blast/')
api.add_resource(RestBlastMinimal, '/blast/<blast_id>')
api.add_resource(RestBlastHits, '/blast/hits/<blast_id>')
api.add_resource(RestBlastPrimers, '/blast_primers/')
api.add_resource(RestNucleotide, '/nucleotide/')
api.add_resource(RestNucleotideMinimal, '/nucleotide/<accession>')
api.add_resource(RestDesignPrimers, '/design/')
api.add_resource(RestShutdown, '/shutdown/')


def start(host='0.0.0.0', debug=False):
    tools.create_empty_database()
    app.run(host=host, debug=debug)


if __name__ == '__main__':
    start()
