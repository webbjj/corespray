import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="corespray",
    version="0.2.dev1",
    author="Steffani Grondin & Jeremy J. Webb",
    author_email="grondin@astro.utoronto.ca, webb@astro.utoronto.ca",
    description="A python packaged for sampling a distribution function for stars that have been ejected from a star cluster's core",
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='MIT',
    packages=["corespray","corespray/df","corespray/potential"],
    setup_requires=['numpy>=1.8','scipy'],
    install_requires=['matplotlib','galpy'],
    )
