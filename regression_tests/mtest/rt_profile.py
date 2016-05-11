#!/usr/bin/env python
"""
Profile the regression_test manager:


cd ${PFLOTRAN_DIR}/regression_tests
python \
    -m cProfile \
    -o rt.profile \
    regression_tests.py \
        -e ../src/pflotran/pflotran \
        --mpiexec /opt/local/bin/mpiexec \
        --test fail-discrete-mean \
        --config-files mtest/expected-fail.cfg
"""

import pstats
p = pstats.Stats("rt.profile")
#p.strip_dirs().sort_stats(-1).print_stats()

# identify long algorithms
p.sort_stats('cumulative').print_stats(10)
# identify long functions
p.sort_stats('time').print_stats(10)

