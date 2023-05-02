# lvlspy_tutorial
This repository contains two Jupyter notebooks that utilize the lvlspy API found [here](https://github.com/jaadt7/lvlspy). 
# Creating XMLs
This notebook will show you how to create energy levels and transitions and write them to XML format, you can run it on google colab [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jaadt7/lvlspy_tutorial/blob/master/create_xml_notebook.ipynb)

# Validation and Calculation
This notebook will allow you to read in generated xmls and perform time evolution calculations. You can load up on the notebook on google colab [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jaadt7/lvlspy_tutorial/blob/master/validate_and_calculate.ipynb)

# ENSDF to XML
The jupyter notebook does exactly that, reads in ENSDF files for specified species and outputs them as a collection to xml format utilizing the libnucnet syntax for species labeling. You would have to specify where you downloaded and extracted your database, which can be downloaded from [here](https://www.nndc.bnl.gov/ensdfarchivals/). The notebook will fill in any downward transitions not found within the files according to a Weisskopf prescription. It also outputs a log file with a txt extension to indicate how many levels and transitions were read from the files and what was excluded and why. You can load up the notebook on google colab [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jaadt7/lvlspy_tutorial/blob/master/ENSDF_to_xml.ipynb)
