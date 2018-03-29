import java.io.PrintWriter

import org.nd4j.linalg.api.ndarray.INDArray
import org.nd4j.linalg.factory.Nd4j
import org.nd4j.linalg.inverse

class polyc(monoCrystal:monoc) {
  //sets up the polycrystal
  //n1, n2, and n3 are used to create a distribution of orientations
  //they correspond to the x, y, and z dimensions
  //n1 = 13, n2 = 7, n3 = 6
  val n1:Int = 7
  val n2:Int = 4
  val n3:Int = 3
  val nbg:Int = n1*n2*n3
  val tetraAS:Double = monoCrystal.AS
  var angles:INDArray = Nd4j.create(nbg,1,3,1)
  anglesetter()
  //angles array is set up, now we turn it into an array of grain orientation matrixes
  var M_rot:INDArray = Nd4j.create(nbg,1,3,3)
  for(i<-0 until nbg){
    var toPut = MME.matrot(angles.getRow(i).getRow(0))
    M_rot.getRow(i).putRow(0,toPut)
  }

  //transfer over the monocrystal data, and put the grains into alignment with monocrystal orientation
  var dir100:INDArray = mtimesx.mul(M_rot,monoCrystal.dir100,'N')
  var dir110:INDArray = mtimesx.mul(M_rot,monoCrystal.dir110,'N')
  var dir111:INDArray = mtimesx.mul(M_rot,monoCrystal.dir111,'N')

  var Pol:INDArray = dir100.mul(monoCrystal.Pol0)
  var epsferro61:INDArray = MME.voigt_33_61(mtimesx.mul(mtimesx.mul(M_rot,monoCrystal.epsferro33,'N'),M_rot,'T'))
  var kappa33:INDArray = mtimesx.mul(mtimesx.mul(M_rot,monoCrystal.kappa33,'N'),M_rot,'T')
  var d36:INDArray = MME.voigt_333_36(dHelper())
  var C66:INDArray = MME.voigt_3333_66(cHelper())
  var invC66:INDArray = Nd4j.create(nbg,monoCrystal.nbdom,6,6)
  for(n<-0 until nbg){
    for(m<-0 until monoCrystal.nbdom){
      invC66.getRow(n).putRow(m,inverse.InvertMatrix.invert(C66.getRow(n).getRow(m),false))
    }
  }
  var ceffMeca = C66.getRow(0).getRow(0)
  var ceffElec = kappa33.getRow(0).getRow(0)
  var semiAxis = Array(Array(1),Array(1),Array(1))
  var eshMeca = MME.Eshelby66_muphylin(ceffMeca,1,semiAxis)
  var eshElec = MME.Eshelby33_muphylin(ceffElec,1)

//  println("ceffMeca",ceffMeca)
//  println("eshMeca",eshMeca)
//  println("ceffElec",ceffElec)
//  println("eshElec",eshElec)
  var cStarMeca = MME.matRep(ceffMeca.mmul(inverse.InvertMatrix.invert(eshMeca,false)).sub(ceffMeca),nbg,1)
  var cStarElec = MME.matRep(ceffElec.mmul(inverse.InvertMatrix.invert(eshElec,false)).sub(ceffElec),nbg,1)

  eshMeca = MME.matRep(eshMeca,nbg,1)
  eshElec = MME.matRep(eshElec,nbg,1)

  println("Polycrystal setup finished")

  def anglespace(min:Double, max:Double, len:Int): INDArray ={
    //uses linspace to give an array with evenly distributed values across a range
    //cuts out last element
    var out = Nd4j.create(len)
    var space = Nd4j.linspace(min, max, len + 1)
    for(i <- 0 until len){
      out.putScalar(i,space.getDouble(i))
    }
    out
  }
  def anglesetter(): Unit ={
    //constructs the array which holds the orientation of all polycrystal grains
    val angle1 = anglespace(0,2*math.Pi,n1)
    val angle2 = anglespace(0,1,n2)
    for(i<-0 until n2){
      angle2.putScalar(i,math.acos(angle2.getDouble(i)))
    }
    val angle3 = anglespace(0,math.Pi/2,n3)
    var count = 0
    var toPut = Nd4j.create(3)
    for(i<- 0 until n1){
      for(j<- 0 until n2){
        for(k<- 0 until n3){
          toPut = Nd4j.create(Array(angle1.getDouble(i),angle2.getDouble(j),angle3.getDouble(k)))
          this.angles.getRow(count).getRow(0).putColumn(0,toPut)
          count += 1
        }
      }
    }
  }
  def dHelper():INDArray = {
    var out: INDArray = null
    var dFill = Nd4j.create(nbg,monoCrystal.nbdom,3,3,3)
    var rotFill = Nd4j.create(nbg,monoCrystal.nbdom,3,3)
    for(i<-0 until nbg){
      dFill.putRow(i,monoCrystal.dmonoc333.getRow(0))
      for(j<-0 until monoCrystal.nbdom){
        rotFill.getRow(i).putRow(j,M_rot.getRow(i).getRow(0))
      }
    }
    out = MME.rotate_ordre3(dFill,rotFill)
    out
  }
  def cHelper(): INDArray ={
    var cFill = Nd4j.create(nbg,monoCrystal.nbdom,3,3,3,3)
    var rotFill = Nd4j.create(nbg,monoCrystal.nbdom,3,3)
    for(i<-0 until nbg){
      for(j<-0 until monoCrystal.nbdom){
        cFill.getRow(i).putRow(j,monoCrystal.Cmonoc3333.getRow(0).getRow(j))
        rotFill.getRow(i).putRow(j,M_rot.getRow(i).getRow(0))
      }
    }
    var out = MME.rotate_ordre4(cFill,rotFill)
    out
  }
  def data_write():Unit = {
    var pw = new PrintWriter("DataFiles/poly/kappa33.txt")
    for(grain<-0 until this.kappa33.shape()(0)){
      for(dir<-0 until this.kappa33.shape()(1)){
        for(i<-0 until this.kappa33.shape()(2)){
          for(j<-0 until this.kappa33.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.kappa33.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/poly/epsferro61.txt")
    for(grain<-0 until this.epsferro61.shape()(0)){
      for(dir<-0 until this.epsferro61.shape()(1)){
        for(i<-0 until this.epsferro61.shape()(2)){
          for(j<-0 until this.epsferro61.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.epsferro61.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/poly/d36.txt")
    for(grain<-0 until this.d36.shape()(0)){
      for(dir<-0 until this.d36.shape()(1)){
        for(i<-0 until this.d36.shape()(2)){
          for(j<-0 until this.d36.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.d36.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/poly/C66.txt")
    for(grain<-0 until this.C66.shape()(0)){
      for(dir<-0 until this.C66.shape()(1)){
        for(i<-0 until this.C66.shape()(2)){
          for(j<-0 until this.C66.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.C66.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/poly/eshMeca.txt")
    for(grain<-0 until this.eshMeca.shape()(0)){
      for(dir<-0 until this.eshMeca.shape()(1)){
        for(i<-0 until this.eshMeca.shape()(2)){
          for(j<-0 until this.eshMeca.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.eshMeca.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/poly/eshElec.txt")
    for(grain<-0 until this.eshElec.shape()(0)){
      for(dir<-0 until this.eshElec.shape()(1)){
        for(i<-0 until this.eshElec.shape()(2)){
          for(j<-0 until this.eshElec.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.eshElec.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/poly/cStarMeca.txt")
    for(grain<-0 until this.cStarMeca.shape()(0)){
      for(dir<-0 until this.cStarMeca.shape()(1)){
        for(i<-0 until this.cStarMeca.shape()(2)){
          for(j<-0 until this.cStarMeca.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.cStarMeca.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/poly/cStarElec.txt")
    for(grain<-0 until this.cStarElec.shape()(0)){
      for(dir<-0 until this.cStarElec.shape()(1)){
        for(i<-0 until this.cStarElec.shape()(2)){
          for(j<-0 until this.cStarElec.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.cStarElec.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/poly/M_rot.txt")
    for(grain<-0 until this.M_rot.shape()(0)){
      for(dir<-0 until this.M_rot.shape()(1)){
        for(i<-0 until this.M_rot.shape()(2)){
          for(j<-0 until this.M_rot.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.M_rot.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/poly/angles.txt")
    for(grain<-0 until this.angles.shape()(0)){
      for(dir<-0 until this.angles.shape()(1)){
        for(i<-0 until this.angles.shape()(2)){
          for(j<-0 until this.angles.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.angles.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/poly/invC66.txt")
    for(grain<-0 until this.invC66.shape()(0)){
      for(dir<-0 until this.invC66.shape()(1)){
        for(i<-0 until this.invC66.shape()(2)){
          for(j<-0 until this.invC66.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.invC66.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
  }
}