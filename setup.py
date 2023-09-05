from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    # Basic package information:
    name="gwlss",
    version="0.1",
    packages=find_packages(),

    # Package metadata:
    author="Richard Stiskalek",
    author_email="richard.stiskalek@protonmail.com",
    description="GW-LSS..",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Richard-Sti/GW-LSS",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=["numpy", "matplotlib", "ipympl"],
)