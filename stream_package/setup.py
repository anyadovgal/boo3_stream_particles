from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="stream_package",
    version="0.0.1",
    description="A place to store shortcut classes & functions",
    author="Anya Dovgal",
    author_email="anyadovgal1@gmail.com",
    packages=find_packages(),
)
