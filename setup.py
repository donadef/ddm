#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=6.0', 'pip==9.0.1', 'flake8==3.5.0', 'coverage==4.5.1', 'Sphinx==1.7.1', 'alchemlyb==0.1.0']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', ]

setup(
    author="Dona de Francquen",
    author_email='donat.defrancquen@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="ddm allows binding free energy for a complex, using the double decoupling method.",
    entry_points={
        'console_scripts': [
            'ddm=ddm.ddm_cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='ddm',
    name='ddm',
    packages=find_packages(include=['ddm']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/donadef/ddm',
    version='0.1.0',
    zip_safe=False,
)
