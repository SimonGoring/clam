clam
====

An R package for the age-depth modelling software clam by Maarten Blaauw.

Clam ([Blaauw, 2010](http://dx.doi.org/10.1016/j.quageo.2010.01.002); this version 2.2) was written to perform 'classic' age-depth modelling - prior to applying more sophisticated techniques such as Bayesian age-depth modelling ([Blaauw and Christen, 2011](http://ba.stat.cmu.edu/journal/2011/vol06/issue03/christen.pdf)).  The original clam files can be found on Maarten Blaauw's [clam webpage](http://chrono.qub.ac.uk/blaauw/clam.html).  This work represents an effort to put clam into a package so that it can work on any file, or on a set of vector data (or a `data.frame`), removing the necessity to have the clam source files in specific locations.

### Publications
+ [Maarten Blaauw](http://chrono.qub.ac.uk/blaauw/) - Queen's University - Belfast, School of Geography, Archaeology and Palaeoecology.

Blaauw, M., 2010. Methods and code for 'classical' age-modelling of radiocarbon sequences. *Quaternary Geochronology* **5**: 512-518

### Package Development
+ [Maarten Blaauw](http://chrono.qub.ac.uk/blaauw/) - Queen's University - Belfast, School of Geography, Archaeology and Palaeoecology.
+ [Simon Goring](http://downwithtime.wordpress.com) - University of Wisconsin-Madison, Department of Geography

### Install `clam`:

+ Development version from GitHub:

```coffee
install.packages("devtools")
require(devtools)
install_github("clam", "SimonGoring")
require(clam)
```

#  Example
This example uses a built in function in the package `netoma` to build a `clam` compatible `csv` file for use in age modeling.  Otherwise the user must generate their own based on the standard discussed in the [clam manual](http://chrono.qub.ac.uk/blaauw/clam.html).  The software package Tilia also provides an option to export clam style data files for modeling.  Alternately, you can use this example as a template.

```coffee
#  Requires the neotoma package:
require(devtools)
install_github("neotoma", "ropensci")
require(neotoma)

marion.site <- get_site(sitename='Marion Lake%')
marion.data <- get_dataset(siteid=marion.site$siteid, datasettype = 'pollen')

marion.download <- get_download(marion.data[[1]]$DatasetID)

#  You need to have a 'Cores' directory in the current working directory.

if(!'Cores' %in% list.files(include.dirs=TRUE)){
  dir.create('Cores')
}

write_agefile(marion.download, path = '.', corename = 'Marion', cal.prog = 'Clam')

#  Build a model using a smooth spline:
clam('Marion',type = 4)

#  The Cores/Marion folder now has a number of files.
```