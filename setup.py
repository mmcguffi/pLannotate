from glob import glob
import os
import re
from setuptools import setup, find_packages
from setuptools import Distribution
from setuptools.command.install import install
import shutil
import sys

__pkg_name__ = 'plannotate'
__author__ = 'Matt McGuffie'
__author_email__ = 'mmcguffi@gmail.com'
__description__ = 'Webserver and command line tool for annotating engineered plasmids'

# Use readme as long description and say its github-flavour markdown
from os import path
this_directory = path.abspath(path.dirname(__file__))
kwargs = {'encoding':'utf-8'} if sys.version_info.major == 3 else {}
with open(path.join(this_directory, 'README.md'), **kwargs) as f:
    __long_description__ = f.read()
__long_description_content_type__ = 'text/markdown'

__path__ = os.path.dirname(__file__)
__pkg_path__ = os.path.join(os.path.join(__path__, __pkg_name__))

# Get the version number from __init__.py, and exe_path
verstrline = open(os.path.join(__pkg_name__, '__init__.py'), 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}/__init__.py".'.format(__pkg_name__))

dir_path = os.path.dirname(__file__)
install_requires = []
with open(os.path.join(dir_path, 'requirements.txt')) as fh:
    reqs = (
        r.split('#')[0].strip()
        for r in fh.read().splitlines() if not r.startswith('#')
    )
    for req in reqs:
        if req == '':
            continue
        if req.startswith('git+https'):
            req = req.split('/')[-1].split('@')[0]
        install_requires.append(req)
print(install_requires)

setup(
    name=__pkg_name__,
    version=__version__,
    author=__author__,
    author_email=__author_email__,
    description=__description__,
    long_description=__long_description__,
    long_description_content_type=__long_description_content_type__,
    install_requires=install_requires,
    tests_require=[].extend(install_requires),
    python_requires='==3.7.*',
    packages=find_packages(exclude=['*.test', '*.test.*', 'test.*', 'test']),
    package_data={__pkg_name__:["data/**/*"]},
    zip_safe=False,
    test_suite='discover_tests',
    entry_points={
        'console_scripts': [
            'plannotate = {}.pLannotate:main'.format(__pkg_name__),
        ] # leave this here
    },
)
