from setuptools import setup

setup(
    name="template_matcher",
    version="0.1.0",
    author="Raymund Hackett",
    author_email="r.e.hackett@lumc.nl",
    description="Detect catalytic enzyme residues in protein structures by matching a library of known templates",
    packages=["template_matcher", "template_matcher.config"],
    package_data={
        "template_matcher": ["config/rules.yml"],
    },
    include_package_data=True,
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=[
        "pyjess>=0.2.1",
        # List any other dependencies your module requires
    ],
)
