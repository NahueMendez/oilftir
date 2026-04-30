from setuptools import setup, find_packages

setup(
    name="oilftir",
    version="0.1.0",
    author="Nahuel Mendez",
    author_email="nahueldanielmendez@gmail.com",
    description="FTIR-based condition monitoring of lubricating oils — ASTM E2412 toolkit",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
        "numpy",
        "pandas",
        "scipy",
        "matplotlib",
        "pybaselines",
    ],
    extras_require={
        "specio": ["specio"],
    },
)