from distutils.core import setup

setup(
    name="dmcalc",
    version="0.2.0",
    description="Package for estimating dispersion measures from radio pulsar observations",
    author="MA Krishnakumar",
    author_email="kkambalappat@gmail.com",
    url="http://github.com/abhisrkckl/dmcalc",
    scripts=["scripts/dmcalc.py", "scripts/run_dmcalc.py"],
)
