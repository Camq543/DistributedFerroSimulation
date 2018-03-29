name := "FerroIntelliJ"

version := "0.1"

scalaVersion := "2.11.12"

val nd4jVersion = "0.9.1"

libraryDependencies += "org.nd4j" % "nd4j-native-platform" % nd4jVersion
libraryDependencies += "org.nd4j" %% "nd4s" % nd4jVersion

libraryDependencies += "org.vegas-viz" %% "vegas" % "0.3.11"