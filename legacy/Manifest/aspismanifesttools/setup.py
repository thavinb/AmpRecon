from setuptools import setup, find_packages
setup(
    name="aspismanifesttools",
    version="0.1.10",
    packages=find_packages(),
    
    install_requires=['mysql-connector-python==8.0.14'],
    
    author='Jon Keatley',
    author_email='jk23@sanger.ac.uk',
    description='Provides a set of aspis manifest tools',
    project_urls={'Source Code':'https://github.com/malariagen/aspis-pipeline'},
    
    entry_points={'console_scripts':[
     'aspis-tools=aspismanifesttools:execCmd']}
)
