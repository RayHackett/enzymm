from . import test_jess_run, test_template, test_utils, test_cli


def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_jess_run))
    suite.addTests(loader.loadTestsFromModule(test_template))
    suite.addTests(loader.loadTestsFromModule(test_utils))
    suite.addTests(loader.loadTestsFromModule(test_cli))
    return suite
