from setuptools import setup, find_packages
setup(
    name="aspisdb",
    version="0.1.1",
    packages=find_packages(),

    install_requires=['mysql-connector-python==8.0.14'],

    author='Jon Keatley',
    author_email='jk23@sanger.ac.uk',
    description='Provides a set of tools for interacting with aspis related databases',
    project_urls={'Source Code':'https://gitlab.com/malariagen-aspis/aspis-store'}
)
