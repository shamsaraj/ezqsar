# ezqsar
ezqsar_f {ezqsar}	
ezqsar_f
Description

This function can easily create a MLR-QSAR model from a proper set of compounds.

Usage:

ezqsar_f(SDFfile, activityfile, propertyfield = "title",
  propertyfield_newset = "title", propertyfield_newset2 = "title",
  Nofdesc = 6, correlation = 1, partition = 0.8,
  des_sel_meth = "forward", testset = 0, newdataset = 0,
  newdataset2 = 0, activity = 0, Cutoff = 3)

It can be installed by one of the following methods in R:

1- 

install.packages("devtools")

devtools::install_github("shamsaraj/ezqsar")

2-

setwd("D:/")#set a working directory

download.file(url="https://github.com/shamsaraj/ezqsar", destfile="ezqsar.zip")

unzip("ezqsar.zip")

system("R CMD build ezqsar")

install.packages("ezqsar_0.8.0.tar.gz", repos = NULL)

3-

Download ezqsar-master.zip file from "clone or download" link in https://github.com/shamsaraj/ezqsar

setwd("D:/")#set the path to the download folder

unzip("ezqsar-master.zip")

system("R CMD build ezqsar-master")

install.packages("ezqsar_0.8.0.tar.gz", repos = NULL)

