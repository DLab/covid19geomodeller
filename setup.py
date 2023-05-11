import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='cv19gm',
    version='0.0.2',
    author='Samuel Ropert',
    author_email='sropert@dlab.cl',
    description='Epidemic dynamic modeling and simulation tools',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/DLab/covid19geomodeller',
    project_urls = {
        "Bug Tracker": "https://github.com/DLab/covid19geomodeller/issues"
    },
    license='MIT',
    packages = setuptools.find_packages(),

    install_requires=['requests==2.28.1','numpy==1.23.4','matplotlib==3.7.1','pandas==2.0.0rc0','scipy==1.10.1','datetime==5.1','toml==0.10.2','pygmo','julia'],

    include_package_data=True, 
    #cmdclass={'install': CustomInstallCommand}  # Use the custom install command
    )