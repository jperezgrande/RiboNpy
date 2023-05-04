# RiboNpy workflow
A Python-based workflow that describes the incorporation of ribozyme constraints into metabolic models


## About

The RiboNpy workflow was written and tested with Python 3.10. 

The core libraries essential for the workflow include: cobra (simulation), plotly(figures), pandas, numpy, and related packages. 


## Installation

1. create RiboNpy environment using conda:

```shell
$ conda create -n RiboNpy python=3.10
```

2. install related packages using pip:

```shell 
$ conda activate RiboNpy
$ pip install cobra
$ pip install plotly
$ pip install -U kaleido
$ pip install nbformat
$ pip install ipykernel
$ python -m ipykernel install --user --name RiboNpy --display-name "RiboNpy"
```

## Steps to reproduce the analysis in the publication

Download all data and analysis code from github (directlt download or use git clone). 

 ```shell
$ cd /file path/project save path/
$ git clone https://github.com/jperezgrande/RiboNpy.git
```

Results can be reproduced by executing the following Jupyter Python notebooks:

+ 01.metabolic_cost_iML1515.ipynb
  + Calculate the metabolic cost of amino acid biosynthesis in iML1515.

+ 02.construct_toy_model.ipynb
  + Using COBRApy, construct the toy model considered in the study.

+ 03.RiboNpy_toy_model.ipynb
  + Construct the ribozyme-constrained toy model (rcToy).

+ 04.RiboNpy_iML1515.ipynb
  + Construct the ribozyme-constrained iML1515 model (rciML1515).
  
+ 05.RiboNpy_simulations.ipynb
  + Results of the simulations included in the report.
 
