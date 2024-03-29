from setuptools import setup

with open("VERSION") as fin:
    VERSION = fin.read().strip()

with open("README.rst") as fin:
    README = fin.read()

with open("LICENSE") as fin:
    LICENSE = fin.read()

URL = "https://github.com/cjdrake/gvmagic"

DOWNLOAD_URL = ""

CLASSIFIERS = [
    "License :: OSI Approved :: Apache Software License",
    "Programming Language :: Python",
]

setup(
    name="gvmagic",
    version=VERSION,
    author="Chris Drake",
    author_email="cjdrake@gmail.com",
    description="Graphviz IPython magic commands",
    license=LICENSE,
    url=URL,
    download_url=DOWNLOAD_URL,
    classifiers=CLASSIFIERS,
    py_modules=["gvmagic"],
)