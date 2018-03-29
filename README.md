# DistributedFerroSimulation

## Setup Guide

### Scala
* To handle dependencies, we use sbt, whose website can be found here: http://www.scala-sbt.org/  
* Our dependencies can be found in the build.sbt file in the "scala" directory. We use Nd4s for linear algebra and Vegas for visualization.
* After installing sbt, clone this repository into a new directory.   
* To run the project, start a terminal within the "scala" root directory and type "sbt run". More info on how to use sbt can be found here: https://alvinalexander.com/scala/sbt-how-to-compile-run-package-scala-project 
* Note, the methods used to write data (monoc.data_write, polyc.data_write, and polyfunc.data_write) are commented out by default, but rely on paths that you must create/edit yourself if you want to use them.

### Python
* Our python script uses numpy for linear algebra and matplotlib for visualization. Make sure you have installed both of these with pip before running the script. 
* Running python is easy. From a terminal, cd into the python root directory and type "python ferroScript.py".
* As with scala, the methods to write data to file are commented out by default, but the path to write to must be changed in order to use these methods. 

## How it works
* Based on various input coefficients (found in constants.py or monoc.scala), the simulation constructs a representation of the ferroelectric crystal across 6 different domains. This work is done in the monocrystal class (monocrystal structure in python and monoc.scala in scala). 
* The polycrystal class (polycrystal structure in python and polyc.scala in scala) creates an isotropic distribution of crystal grains, represented by directional arrays. Using these arrays, the data in from the monocrystal is rotated to each grain and stored in the polycrystal class.  
* After constructing the poly crystal, the polyfunc class is used to iterate over each grain and apply different loads in the form of stress and electric field. The behavior of each individual grain is averaged to output curves describing the behavior of the entire polycrystal.
