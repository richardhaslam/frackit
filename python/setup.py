from setuptools import setup, find_namespace_packages

setup(
    name = "frackit",
    version = "1.0.0",
    author = "Dennis Glaeser",
    author_email = "dennis.glaeser@iws.uni-stuttgart.de",
    description = "Frackit Python library",
    long_description = "",
    packages = find_namespace_packages(include=["frackit.*"]),
    package_data = {"": ["*.so"]},
    zip_safe = False,
)
