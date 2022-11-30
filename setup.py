from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

setup(
    name="Newton-Polygon",
    version="0.0.5",
    author="Ivan Kuznetsov",
    author_email="maxbeyn@icloud.com",
    description="Python power geometry package.",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/Garrisonmond/Newton-Polygon/tree/main",
    packages=find_packages(),
    install_requires=['matplotlib', 'seaborn', 'wolframclient'],
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)