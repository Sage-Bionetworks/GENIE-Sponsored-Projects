"""genie package setup"""
import os
from setuptools import setup, find_packages

# figure out the version
# about = {}
# here = os.path.abspath(os.path.dirname(__file__))
# with open(os.path.join(here, "genie", "__version__.py")) as f:
#     exec(f.read(), about)

# Add readme
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='geniesp',
      version="0.0.1",
      description='Export cBioPortal files for BPC',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/Sage-Bionetworks/GENIE-Sponsored-Projects',
      author='Thomas Yu',
      author_email='thomas.yu@sagebionetworks.org',
      license='MIT',
      packages=find_packages(),
      zip_safe=False,
      python_requires='>=3.6',
      entry_points={'console_scripts': ['geniesp = geniesp.__main__:main']},
      install_requires=['aacrgenie>=12.5.0',
                        'synapseclient>=2.3.0',
                        'pandas'])
