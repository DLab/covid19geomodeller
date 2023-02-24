import setuptools

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
    project_urls = {
        "Bug Tracker": "https://github.com/DLab/covid19geomodeller/issues"
    },
    license='MIT',
    packages = setuptools.find_packages(),
    install_requires=['requests','numpy','matplotlib','pandas','scipy','datetime','toml'], # 'pygmo'
    include_package_data=True 
)

#packages = ['cv19gm','cv19gm.models']