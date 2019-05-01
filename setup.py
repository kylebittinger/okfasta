from setuptools import setup

setup(name='okfasta',
      version='0.0.1',
      description='FASTA utilities',
      author='Kyle Bittinger',
      author_email='kylebittinger@gmail.com',
      url='https://github.com/kylebittinger/okfasta',
      packages=['okfastalib'],
      entry_points = {
          'console_scripts': [
              'okfasta=okfastalib.command:main',
          ],
      }
)
