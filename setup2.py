import setuptools
from setuptools.command.install import install

def post_install(func):
    def wrapped(cls):
        orig_run = cls.run

        def new_run(self):
            orig_run(self)
            func(self)
        cls.run = new_run
        return cls

    return wrapped

@post_install
def run_julia_setup(self):
    import julia
    julia.install()

@run_julia_setup
class CustomInstallCommand(install):
    pass

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='cv19gm',
    version='0.0.1',
    author='Samuel Ropert',
    author_email='sropert@dlab.cl',
    description='Epidemic dynamic modeling and simulation tools',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/DLab/covid19geomodeller',
    project_urls={
        "Bug Tracker": "https://github.com/DLab/covid19geomodeller/issues"
    },
    license='MIT',
    packages=setuptools.find_packages(),

    install_requires=['requests==2.28.1', 'numpy==1.23.4', 'matplotlib==3.7.1', 'pandas==2.0.0rc0', 'scipy==1.10.1', 'datetime==5.0', 'toml==0.10.2', 'pygmo', 'julia'],

    include_package_data=True,
    cmdclass={'install': CustomInstallCommand}
)
