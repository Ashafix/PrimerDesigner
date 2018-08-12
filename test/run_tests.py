import unittest
import sys
import os
import io
import time


class TestRunner(unittest.TextTestRunner):

    if __name__ == '__main__':
        file = sys.argv[0]
    else:
        file = __file__
    testdir = os.path.abspath(os.path.dirname(file) or os.curdir)

    def __init__(self, tests=None, verbosity=0):
        if tests is None:
            self.tests = []
        else:
            self.tests = tests
        if not self.tests:
            for filename in os.listdir(TestRunner.testdir):
                if filename.startswith('test_') and filename.endswith('py'):
                    self.tests.append(filename[:-3])
            self.tests.sort()

        stream = io.StringIO()
        unittest.TextTestRunner.__init__(self, stream,
                                         verbosity=verbosity)

    def runTest(self, name):
        result = self._makeResult()
        output = io.StringIO()
        os.chdir(self.testdir)

        sys.stdout = output

        sys.stderr.write('{} ... '.format(name))
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromName(name)
        if hasattr(loader, 'errors') and loader.errors:
            sys.stderr.write('loading tests failed:\n')
            for msg in loader.errors:
                sys.stderr.write('{}\n'.format(msg))
            return False

        suite.run(result)
        if self.testdir != os.path.abspath('.'):
            sys.stderr.write('FAIL\n')
            result.stream.write(result.separator1 + '\n')
            result.stream.write('ERROR: %s\n' % name)
            result.stream.write(result.separator2 + '\n')
            result.stream.write('Current directory changed\n')
            result.stream.write('Was: %s\n' % self.testdir)
            result.stream.write('Now: %s\n' % os.path.abspath('.'))
            os.chdir(self.testdir)
            if not result.wasSuccessful():
                result.printErrors()
            return False
        elif result.wasSuccessful():
            sys.stderr.write('ok\n')
            return True
        else:
            sys.stderr.write('FAIL\n')
            result.printErrors()
        return False

    def run(self):
        failures = 0
        start_time = time.time()
        for test in self.tests:
            ok = self.runTest(test)
            if not ok:
                failures += 1
        total = len(self.tests)
        stop_time = time.time()
        time_taken = stop_time - start_time
        sys.stderr.write(self.stream.getvalue())
        sys.stderr.write('-' * 70 + '\n')
        sys.stderr.write('Ran %d test%s in %.3f seconds\n' %
                         (total, total != 1 and 's' or '', time_taken))
        sys.stderr.write('\n')
        if failures:
            sys.stderr.write('FAILED (failures = %d)\n' % failures)
        return failures


def main(args, verbosity=0):
    runner = TestRunner(args, verbosity)
    return runner.run()


if __name__ == '__main__':
    errors = main(sys.argv[1:])
    if errors:
        # Doing a sys.exit(...) isn't nice if run from IDLE...
        sys.exit(1)