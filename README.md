# vw-Jemez

Functions for running our modified CASiMiR model and D-FLOW coupled
model on a supercomputer.

# Installation

### 1. Clone the repo and cd in to the root directory

```bash
git clone https://github.com/VirtualWatershed/vw-jemez && cd vw-jemez
```

### 2. Use a virtual environment and install dependencies

```bash
virtualenv venv
```

```bash
source venv/bin/activate
```

```bash
pip install -r requirements.txt
```

## Unit tests

To check that all is well, try running the unit tests:

```bash
nosetests -v
```


The output should be

```
test_casimir (test.test_dflow_casimir.TestDflow) ... ok
test_vegmap_properly_read (test.test_dflow_casimir.TestDflow) ... ok

----------------------------------------------------------------------
Ran 2 tests in 0.127s

OK
```
