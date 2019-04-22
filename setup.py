import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
     name='cgstatistic',  
     version='0.1',
     scripts=['cgstatistic'] ,
     author="Dimmy Magalhaes",
     author_email="dimmyk@gmail.com",
     description="A Docker and AWS utility package",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/dimmykarson/CGStatistic",
     packages=setuptools.find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )