Reproducing the ant-torus computations
=====================================

Run the commands below from this directory with all the files together.
The scripts need only Python 3.10 or later and Node.js. They were checked
with Python 3.14.2 and Node.js 24.13.0.

    python -B count_blank_class.py 10 10 --check
    python -B count_blank_class.py 12 12
    python -B verify_ant.py
    python -B orbit_distribution.py

The first command checks the counts against direct enumeration on ten
small tori before printing the exact scale for 10 by 10. The 12 by 12
count took about two seconds in the recorded check. The verifier reports
151 checks. The distribution script reproduces the small-board period
histograms and tests the width-and-back walks.

    node ant.js blank 10 10 150192928160

The last command repeats the saved 10 by 10 computation in
blank-10x10.json. The saved run took about 26 minutes. It returns the
same value for the reduced and full periods. The command below gives
a short check that finishes almost immediately.

    node ant.js blank 4 6 696
