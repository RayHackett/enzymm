from . import test_jess_run, test_template, test_utils


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_jess_run))
    suite.addTests(loader.loadTestsFromModule(test_template))
    suite.addTests(loader.loadTestsFromModule(test_utils))
    return suite
