# Tipping point  analysis

This file is a jupyter notebook describes how to conduct the analysis and recreate the figures given in: [Shah I, Setzer RW, Jack J, Houck KA, Judson RS, Knudsen TB, Liu J, Martin MT, Reif DM, Richard AM, Thomas RS, Crofton KM, Dix DJ, Kavlock RJ. 2016. Using ToxCast™ data to reconstruct dynamic cell state trajectories and estimate toxicological points of departure. Environ Health Perspect 124:910–919](http://dx.doi.org/10.1289/ehp.1409029). 

# Installation
All requirements for running this code are open-source:-

## Operating system
This code has only been tested under [Red Hat Enterprise Linux v6](https://www.redhat.com/en/technologies/linux-platforms/enterprise-linux) and [Ubuntu Server 16.04 LTS](http://www.ubuntu.com). The Ubuntu system can be downloaded freely and installed on any Intel PC. 

## High-throughput screening database (hts-db)
The [hts-db system](https://github.com/i-shah/hts-db) contains the data and functions necessary for replicating the results in this notebook. It is freely available from GitHub. The hts-db system uses Mongo nosql database to store all high-throughput screening (HTS) data from the ToxCast HepG2 high-content imaging (HCI) study including: raw well level data, plate-level normalized data and lowest effect concentration data. The can also be used to manage the raw image data and this capability will be released in future updates. You must download the hts-db system and follow the installation instructions before running this notebook.


## Python 
This code has been tested on [Python 2.7](http://python.org). The easiest way to install Python 2.7 is via [Anaconda](https://www.continuum.io/downloads). Install the python packages given in py_requirements.txt as follows:

```
pip install -r py_requirements.txt
```

* Python library: The code depends on a python source file (traj.py) that must be installed locally under your PYTHONPATH so that it can be imported by the interpreter. The python source files required for this notebook should be organized as follows (the source files installed via the hts-db system are included): 

* Jupyter: [Jupyter notebook](http://jupyter.org/) is an interactive computing environment based and is freely available. It can be installed using Anaconda. Please the the documents about Jupyter or find a tutorial on YouTube. 
After jupyter is installed read the [quickstart instructions](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/) to create a configuration file. Set the following variables in your jupyter_notebook_config.py:
```
c.NotebookApp.notebook_dir = your_notebook_directory
c.NotebookApp.port = 7777
```
Start the notebook server from the command line using (I always run this in a [screen](https://www.gnu.org/software/screen/manual/screen.html)): 
jupyter notebook

Note: If another service is running on port 7777 you will have to change the port number. 

Once you have jupyter running load the notebook (Shah-Tipping-Point-EHP-2016-v0.2.ipynb) and follow the instructions.





