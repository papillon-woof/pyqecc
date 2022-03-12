from setuptools import setup
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="pyqecc",  # package name
    packages=[
        "pyqecc",
        "pyqecc.qecc",
        "pyqecc.sim",
        "pyqecc.channel",
        "pyqecc.util",
    ],
    version="0.0.2",  # version
    license="MIT",  # lisence
    install_requires=[
        "numpy>=1.21",
    ],
    author="papillon",
    author_email="ogyogugige@gmail.com",
    url="https://github.com/papillon-woof/pyqecc",
    description="Simulater for quantum error correction code.",
    long_description_content_type="text/markdown",
    keywords="pyqecc PyQecc py-qecc qecc qec quantum error correction quantum computer quantum computing",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.9",
    ],
)
