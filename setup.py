from setuptools import setup, Extension, Command
import os,sys,glob

fast_miscentering_module = Extension('_fast_miscentering',
                            extra_compile_args=[os.path.expandvars("-I${GSLI}")],
                            extra_link_args=[os.path.expandvars("-L${GSLL}"),"-lm","-lgsl"],
                            sources=["fast_miscentering.c"])

def read(fname):
    """Quickly read in the README.md file."""
    return open(os.path.join(os.path.dirname(__file__),fname)).read()

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

setup(name='fast_miscentering',
      install_requires=['numpy'],
      version='1.0',
      py_modules=['fast_miscentering'],
      description='A fast calculation of a miscentered Sigma(R).',
      long_description=read('README.md'),
      author='Tom McClintock',
      author_email='tmcclintock@email.arizona.edu',
      url='https://github.com/tmcclintock/FastMiscentering',
      ext_modules=[fast_miscentering_module],
      cmdclass={'clean': CleanCommand})
