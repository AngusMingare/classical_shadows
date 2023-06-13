from setuptools import setup, find_packages

setup(
    name="classical_shadows",
    version="1.0.0",
    packages=find_packages(exclude=['test']),
    install_requires=["numpy", "pytest", "qiskit"]   
    )