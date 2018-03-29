# Summer-Ferro-IntelliJ

## Setup Guide
* To handle dependencies, we use sbt, whose website can be found here: http://www.scala-sbt.org/  
* After installing sbt, clone this repository into a new directory.   
* Next, go into the "build.sbt" file, and change the library dependency on line 9 to match your system architecture (for example, the current system architecture is set to "windows-x86_64"). This changes the linear algebra backend used by the Nd4s library. Info on backend dependencies can be found here: http://nd4j.org/dependencies  
* Running the project is easy, simply start a terminal within the project root directory (the repository you just cloned), and type "sbt run." More info on how to use sbt can be found here: https://alvinalexander.com/scala/sbt-how-to-compile-run-package-scala-project  

## How it works
* The only internal values that should be changed are the input coefficients found from lines 13-24 of monoc.scala. These define the input data for the crystal  
* In order to change the number of grains, set the n1, n2, and n3 variables found in lines 9-11 of polyc.scala. These define the number of different starting orientations along the x, y, and z axis.  
* From this input, the simulation generates a polycrystal and gives piezoelectric information as output.
