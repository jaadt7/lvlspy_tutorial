# lvlspy_tutorial
This repository contains Jupyter notebooks that demonstrate the lvlspy [API](https://lvlspy.readthedocs.io/en/latest/lvlspy.html). 

# Creating XMLs
This notebook demonstrates how to create levels and transitions, add those to species, add the species to a species collection, and then how to write those data to an XML file.  You can download the notebook and run it locally, or you can run it on Google colab by clicking on the following badge: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jaadt7/lvlspy_tutorial/blob/master/create_xml_notebook.ipynb)

# Validation and Calculation
This notebook will allow you to read in data from lvlspy appropriate XML and perform time evolution calculations. You can download the notebook and run it locally, or you can run it on Google colab by clicking on the following badge: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jaadt7/lvlspy_tutorial/blob/master/validate_and_calculate.ipynb)

# ENSDF to XML
This notebook will you to download and read in [ENSDF files](https://www.nndc.bnl.gov/ensdfarchivals/) and outputs them as a collection to xml format utilizing the libnucnet syntax. It is equipped with a routine that will fill in any downward transitions that are not found within the files according to a Weisskopf prescription. A log file is generated to indicate how many levels and transitions were read from the files, which of those were excluded and why. You can download the notebook and run it locally, or you can run it on on Google Colab by clicking on the following badge: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jaadt7/lvlspy_tutorial/blob/master/ENSDF_to_xml.ipynb)
