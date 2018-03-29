import java.io.PrintWriter

import org.nd4j.linalg.api.ndarray.INDArray
import org.nd4j.linalg.factory.Nd4j

//class to hold monocrystal
//performs setup from donn_monoc_tetra in constructor
class monoc {

  //constants
  val nbdom = 6
  val AS:Double = "1E-5".toDouble
  val Pol0: Double = 0.3
  //set the monocrystal direction matrixes
  //these will be used to set up the polycrystal grains
  var dir100:INDArray = Nd4j.zeros(1,6,3,1)
  dir100.getRow(0).putRow(0,Nd4j.create(Array(Array(1.0),Array(0.0),Array(0.0))))
  dir100.getRow(0).putRow(1,Nd4j.create(Array(Array(-1.0),Array(0.0),Array(0.0))))
  dir100.getRow(0).putRow(2,Nd4j.create(Array(Array(0.0),Array(1.0),Array(0.0))))
  dir100.getRow(0).putRow(3,Nd4j.create(Array(Array(0.0),Array(-1.0),Array(0.0))))
  dir100.getRow(0).putRow(4,Nd4j.create(Array(Array(0.0),Array(0.0),Array(1.0))))
  dir100.getRow(0).putRow(5,Nd4j.create(Array(Array(0.0),Array(0.0),Array(-1.0))))
  var dir110:INDArray = Nd4j.zeros(1,12,3,1)
  dir110.getRow(0).putRow(0,Nd4j.create(Array(Array(1.0),Array(1.0),Array(0.0))))
  dir110.getRow(0).putRow(1,Nd4j.create(Array(Array(1.0),Array(-1.0),Array(0.0))))
  dir110.getRow(0).putRow(2,Nd4j.create(Array(Array(-1.0),Array(1.0),Array(0.0))))
  dir110.getRow(0).putRow(3,Nd4j.create(Array(Array(-1.0),Array(-1.0),Array(0.0))))
  dir110.getRow(0).putRow(4,Nd4j.create(Array(Array(0.0),Array(1.0),Array(1.0))))
  dir110.getRow(0).putRow(5,Nd4j.create(Array(Array(0.0),Array(1.0),Array(-1.0))))
  dir110.getRow(0).putRow(6,Nd4j.create(Array(Array(0.0),Array(-1.0),Array(1.0))))
  dir110.getRow(0).putRow(7,Nd4j.create(Array(Array(0.0),Array(-1.0),Array(-1.0))))
  dir110.getRow(0).putRow(8,Nd4j.create(Array(Array(1.0),Array(0.0),Array(1.0))))
  dir110.getRow(0).putRow(9,Nd4j.create(Array(Array(1.0),Array(0.0),Array(-1.0))))
  dir110.getRow(0).putRow(10,Nd4j.create(Array(Array(-1.0),Array(0.0),Array(1.0))))
  dir110.getRow(0).putRow(11,Nd4j.create(Array(Array(-1.0),Array(0.0),Array(-1.0))))
  dir110.muli(1.0/math.sqrt(2))
  var dir111:INDArray = Nd4j.ones(1,8,3,1)
  dir111.getRow(0).putRow(0,Nd4j.create(Array(Array(1.0),Array(1.0),Array(1.0))))
  dir111.getRow(0).putRow(1,Nd4j.create(Array(Array(1.0),Array(1.0),Array(-1.0))))
  dir111.getRow(0).putRow(2,Nd4j.create(Array(Array(1.0),Array(-1.0),Array(1.0))))
  dir111.getRow(0).putRow(3,Nd4j.create(Array(Array(1.0),Array(-1.0),Array(-1.0))))
  dir111.getRow(0).putRow(4,Nd4j.create(Array(Array(-1.0),Array(1.0),Array(1.0))))
  dir111.getRow(0).putRow(5,Nd4j.create(Array(Array(-1.0),Array(1.0),Array(-1.0))))
  dir111.getRow(0).putRow(6,Nd4j.create(Array(Array(-1.0),Array(-1.0),Array(1.0))))
  dir111.getRow(0).putRow(7,Nd4j.create(Array(Array(-1.0),Array(-1.0),Array(-1.0))))
  dir111.muli(1.0/math.sqrt(3))
  //other arrays
  var dir_ref: INDArray = Nd4j.create(Array(0,0,1.0)).reshape(3,1)
  var Pol: INDArray = this.dir100.mul(this.Pol0)

  //setup Tensors
  def tensorHelper(): INDArray ={
    //finds the Rotation matrix R which is used to calculate cmonoc3333, kappa33, and dmonoc333
    //math below is to calculate the "R" matrix
    var R = Nd4j.zeros(1,nbdom,3,3)
    var cos_t = 1.0
    var dir_r = Nd4j.create(Array(0.0,0,0)).reshape(3,1)
    var mat_cross = Nd4j.zeros(3,3)
    var toPut = Nd4j.create(3,3)

    for(i <- 0 until nbdom){
      cos_t = dir_ref.transpose().mmul(this.dir100.getRow(0).getRow(i)).getDouble(0)
      dir_r.putColumn(0,Nd4j.create(Array(this.dir100.getDouble(0,i,1,0),-1*this.dir100.getDouble(0,i,0,0),0)))
      if(dir_r.norm2Number().asInstanceOf[Double] != 0){
        dir_r.divi(dir_r.norm2Number())
      }

      mat_cross.putRow(0,Nd4j.create(Array(0,0,dir_r.getDouble(1,0))))
      mat_cross.putRow(1,Nd4j.create(Array(0,0,-1*dir_r.getDouble(0,0))))
      mat_cross.putRow(2,Nd4j.create(Array(-dir_r.getDouble(1,0),dir_r.getDouble(0,0),0)))

      toPut = Nd4j.eye(3).mul(cos_t).add(dir_r.mmul(dir_r.transpose()).mul(1-cos_t)).sub(mat_cross.mul(math.sqrt(1-(cos_t * cos_t))))
      R.getRow(0).putRow(i,toPut)
    }
    R
  }
  var R:INDArray =tensorHelper()
  var epsferro33 = epsferroHelper()
  var kappa = kappaHelper()
  var dmonoc = dmonocHelper()
  var Cmonoc = CmonocHelper()

  //now, with R, we calculate Cmonoc, kappa, and dmonoc, which are used later
  var Cmonoc3333temp:INDArray = MME.voigt_66_3333(Cmonoc)
  var dmonoc333temp:INDArray = MME.voigt_36_333(dmonoc)

  //with each of the tensors calculated, we create arrays of 1 x M x T where
  //T is the tensor matrix Cmonoc, kappa, or Dmonoc
  var Cmonoc3333:INDArray = MME.matRep(Cmonoc3333temp,1,nbdom)
  this.Cmonoc3333 = MME.rotate_ordre4(this.Cmonoc3333,R)

  var kappa33:INDArray = MME.matRep(kappa,1,nbdom)
  this.kappa33 = mtimesx.mul(mtimesx.mul(R,this.kappa33,'N'),R,'T')

  var dmonoc333:INDArray = MME.matRep(dmonoc333temp,1,nbdom)
  this.dmonoc333 = MME.rotate_ordre3(this.dmonoc333,R)
  data_write()

  println("Monocrystal setup finished")



  //SETUP FUNCTIONS
  def CmonocHelper(): INDArray = {
    val C11 = "106E9".toDouble
    val C12 = "62e9".toDouble
    val C44 = "44E9".toDouble
    var Cmonoc = Nd4j.create(Array(Array(C11, C12, C12, 0, 0, 0),
      Array(C12, C11, C12, 0, 0, 0),
      Array(C12, C12, C11, 0, 0, 0),
      Array(0, 0, 0, C44, 0, 0),
      Array(0, 0, 0, 0, C44, 0),
      Array(0, 0, 0, 0, 0, C44)))
    Cmonoc
  }
  def dmonocHelper():INDArray = {
    val d31 = "-2.1E-10".toDouble
    val d33 = "4.5E-10".toDouble
    val d15 = "5.8E-10".toDouble
    var dmonoc = Nd4j.create(Array(Array(0, 0, 0, 0, d15, 0),
      Array(0, 0, 0, d15, 0, 0),
      Array(d31, d31, d33, 0, 0, 0)))
    dmonoc
  }
  def kappaHelper(): INDArray ={
    val k33 = "2.00E-8".toDouble
    val k11 = "2.00E-8".toDouble
    var kappa = Nd4j.create(Array(Array(k11,0,0),
      Array(0,k11,0),
      Array(0,0,k33)))
    kappa
  }
  def epsferroHelper():INDArray = {
    //sets up the epsferro33 array
    val epz0 = "2E-3".toDouble
    val mat_id = MME.matRep(Nd4j.eye(3),1,nbdom)
    mtimesx.mul(this.dir100,this.dir100,'T').mul(3).sub(mat_id).mul(epz0/2)
  }
  def data_write(): Unit ={
    var pw = new PrintWriter("DataFiles/mono/Cmonoc3333.txt")
    for(grain<-0 until this.Cmonoc3333.shape()(0)){
      for(dir<-0 until this.Cmonoc3333.shape()(1)){
        for(i<-0 until this.Cmonoc3333.shape()(2)){
          for(j<-0 until this.Cmonoc3333.shape()(3)){
            for(k<-0 until this.Cmonoc3333.shape()(4)){
              for(l<-0 until this.Cmonoc3333.shape()(5)){
                pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " k: " + k + " l: " + l  + " data: " + this.Cmonoc3333.getDouble(grain, dir, i, j,k,l) + "\n")
              }
            }
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/mono/kappa33.txt")
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
    pw = new PrintWriter("DataFiles/mono/dmonoc333.txt")
    for(grain<-0 until this.dmonoc333.shape()(0)){
      for(dir<-0 until this.dmonoc333.shape()(1)){
        for(i<-0 until this.dmonoc333.shape()(2)){
          for(j<-0 until this.dmonoc333.shape()(3)){
            for (k<-0 until this.dmonoc333.shape()(4)){
              pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " k: " + k +  " data: " + this.dmonoc333.getDouble(grain, dir, i, j,k) + "\n")
            }
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/mono/epsferro33.txt")
    for(grain<-0 until this.epsferro33.shape()(0)){
      for(dir<-0 until this.epsferro33.shape()(1)){
        for(i<-0 until this.epsferro33.shape()(2)){
          for(j<-0 until this.epsferro33.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.epsferro33.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/mono/dir100.txt")
    for(grain<-0 until this.dir100.shape()(0)){
      for(dir<-0 until this.dir100.shape()(1)){
        for(i<-0 until this.dir100.shape()(2)){
          for(j<-0 until this.dir100.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.dir100.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/mono/dir110.txt")
    for(grain<-0 until this.dir110.shape()(0)){
      for(dir<-0 until this.dir110.shape()(1)){
        for(i<-0 until this.dir110.shape()(2)){
          for(j<-0 until this.dir110.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.dir110.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    pw = new PrintWriter("DataFiles/mono/dir111.txt")
    for(grain<-0 until this.dir111.shape()(0)){
      for(dir<-0 until this.dir111.shape()(1)){
        for(i<-0 until this.dir111.shape()(2)){
          for(j<-0 until this.dir111.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + this.dir111.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
    val mat_id = MME.matRep(Nd4j.eye(3),1,nbdom)
    var temp = mtimesx.mul(this.dir100,this.dir100,'T').mul(3).sub(mat_id)
    pw = new PrintWriter("DataFiles/mono/dir100Tmul.txt")
    for(grain<-0 until temp.shape()(0)){
      for(dir<-0 until temp.shape()(1)){
        for(i<-0 until temp.shape()(2)){
          for(j<-0 until temp.shape()(3)){
            pw.write("grain: " + grain + " dir: " + dir + " i: " + i + " j: " + j + " data: " + temp.getDouble(grain, dir, i, j) + "\n")
          }
        }
      }
    }
    pw.close()
  }
}