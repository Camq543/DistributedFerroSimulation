import org.nd4j.linalg.api.ndarray.INDArray
import org.nd4j.linalg.factory.Nd4j
import java.io._



import vegas._
import vegas.render.WindowRenderer._

import math.sqrt
import math.cos
import math.sin

object MME{

  def main(args:Array[String]){
    //Set up Monocrystal
    var monoCrystal = new monoc()
//    monoCrystal.data_write()
    //Set up Polycrystal
    var polyCrystal = new polyc(monoCrystal)
//    polyCrystal.data_write()
    //plot directionality (unnecessary, only for checking)
    //first pass on polycrystal to determine mechanical and electric tensors
    var n = 51
//    var n = 1
    var E_value = Nd4j.linspace(0,2e6,n)
    var T_value = Nd4j.create(Array(-100.0,-50.0,0.0,50.0,100.0)).mul(1e6)
//    var E_value = Nd4j.create(Array(1e6))
//    var T_value = Nd4j.create(Array(10e7))
    var m = T_value.shape()(1)

    var polyc_D = Nd4j.zeros(n,7,3,1)
    var polyc_S = Nd4j.zeros(n,7,6,1)

    var Tmacro = Nd4j.zeros(6,1)
    var Emacro = Nd4j.zeros(3,1)
    var Tloc = Nd4j.zeros(polyCrystal.nbg,1,6,1)
    var Eloc = Nd4j.zeros(polyCrystal.nbg,1,3,1)
    //val pw = new PrintWriter(new File("iterCount.txt"))
    for(k <- 0 until m){
      Tmacro.putScalar(Array(2,0),T_value.getDouble(k))
      Tloc = matRep(Tmacro,polyCrystal.nbg,1)
      Eloc = Nd4j.zeros(polyCrystal.nbg,1,3,1)
      for(i<- 0 until n){
        Emacro.putScalar(Array(2,0),E_value.getDouble(i))
        println("i",i,"k",k)
        var retList:List[INDArray] = polyfunc.polyfunc(polyCrystal,Tmacro,Emacro,Tloc,Eloc)
        Tloc = retList(2)
        Eloc = retList(3)
        polyc_S.getRow(i).putRow(k,retList(0))
        polyc_D.getRow(i).putRow(k,retList(1))
//        try{
//          pw.write("i: " + i + " k: " + k + " iterations: " + retList(4).getDouble(0) + "\n")
//          pw.flush()
//        }
//        catch{
//          case iob: IndexOutOfBoundsException => println("ERROR, no iteration count")
//        }
      }
    }
    //pw.close()
    var plotData = List[Map[String,Double]]()
    var xval = 1.0
    var yval = 1.0
    plotData = List[Map[String,Double]]()
    for(i<-0 until n){
      for(k<-0 until m){
        xval = E_value.getDouble(i)
        yval = polyc_D.getRow(i).getRow(k).getRow(2).getDouble(0)
        //println("Xvalue",xval,"Yvalue",yval)
        plotData = plotData :+ Map("a" -> xval,"b" -> yval, "c" -> k.asInstanceOf[Double])
      }
    }
    var plot1 = Vegas("poly_D scatter.",width = 600.0,height = 600.0).
      withData(plotData).
      encodeX("a", Quantitative).
      encodeY("b", Quantitative).
      encodeColor("c",Nominal).
      mark(Line)
    plot1.show

    plotData = List[Map[String,Double]]()
    for(i<-0 until n){
      for(k<-0 until m){
        xval = E_value.getDouble(i)
        yval = polyc_S.getRow(i).getRow(k).getRow(2).getDouble(0) - polyc_S.getRow(0).getRow(k).getRow(2).getDouble(0)
        //println("Xvalue",xval,"Yvalue",yval)
        plotData = plotData :+ Map("a" -> xval,"b" -> yval, "c" -> k.asInstanceOf[Double])
      }
    }
    var plot2 = Vegas("Poly_S scatter",width=800.0,height=800.0).
      withData(plotData).
      encodeX("a", Quantitative).
      encodeY("b", Quantitative).
      encodeColor("c",Nominal).
      mark(Line)
    plot2.show

//    val nbT = 16
//    val nbE = 16
//    val alph = 1e-3
//    T_value = Nd4j.linspace(0,-150,nbT)
//    E_value = Nd4j.linspace(0,5,nbE)
//    var LApp = Nd4j.zeros(9,9,nbT,nbE)
//    var Tinit = Nd4j.zeros(polyCrystal.nbg,1,6,1)
//    var Einit = Nd4j.zeros(polyCrystal.nbg,1,3,1)
//    //hang onto tloc and eloc to remain as zero matrices
//    Tloc = Nd4j.zeros(polyCrystal.nbg,1,6,1)
//    Eloc = Nd4j.zeros(polyCrystal.nbg,1,3,1)
//    var mecVar = Nd4j.zeros(6,1)
//    var elecVar = Nd4j.zeros(3,1)
//    var TStat = Nd4j.zeros(6,1)
//    var EStat = Nd4j.zeros(3,1)
//    //need a bunch of arrays for storage
//    var S1:INDArray = null
//    var S2:INDArray = null
//    var D1:INDArray = null
//    var D2:INDArray = null
//    var T:INDArray = null
//    var E:INDArray = null
//    var delta = 0.0
//    var essai = 0.0
//
//    for(i<-0 until nbT){
//      TStat.putScalar(2,0,T_value.getDouble(i))
//      for(k<-0 until nbE){
//        EStat.putScalar(2,0,E_value.getDouble(k))
//        println("i",i,"k",k)
//        var retList:List[INDArray] = polyfunc.polyfunc(polyCrystal,TStat,EStat,Tloc,Eloc)
//        Tinit = retList(2)
//        Einit = retList(3)
//        for(j<- 0 until 9){
//          T = TStat
//          E = EStat
//          if(j<6){
//            mecVar = Nd4j.zeros(6,1)
//            if(T_value.getDouble(i) == 0){
//              delta = 100.0
//            }
//            else{
//              delta = -alph * T_value.getDouble(i)
//            }
//            mecVar.putScalar(j,delta)
//            T = TStat.add(mecVar)
//            retList = polyfunc.polyfunc(polyCrystal,T,E,Tinit.add(matRep(mecVar,polyCrystal.nbg,1)),Einit)
//            S1 = retList(0).getRow(0)
//            D1 = retList(1).getRow(0)
//            T = TStat.sub(mecVar)
//            retList = polyfunc.polyfunc(polyCrystal,T,E,Tinit.sub(matRep(mecVar,polyCrystal.nbg,1)),Einit)
//            S2 = retList(0).getRow(0)
//            D2 = retList(1).getRow(0)
//          }
//          else{
//            elecVar = Nd4j.zeros(3,1)
//            if(E_value.getDouble(k)==0){
//              delta = 100.0
//            }
//            else{
//              delta = alph * E_value.getDouble(k)
//            }
//            elecVar.putScalar(j-6,delta)
//            E = EStat.add(elecVar)
//            retList = polyfunc.polyfunc(polyCrystal,T,E,Tinit,Einit.add(matRep(elecVar,polyCrystal.nbg,1)))
//            S1 = retList(0).getRow(0)
//            D1 = retList(1).getRow(0)
//            E = EStat.sub(elecVar)
//            retList = polyfunc.polyfunc(polyCrystal,T,E,Tinit,Einit.sub(matRep(elecVar,polyCrystal.nbg,1)))
//            S2 = retList(0).getRow(0)
//            D2 = retList(1).getRow(0)
//          }
//          var temp = Nd4j.concat(0,S1.sub(S2).div(2*delta),D1.sub(D2).div(2*delta))
//          for(l<-0 until 9){
//            LApp.putScalar(l,j,i,k,temp.getDouble(l))
//          }
//        }
//      }
//    }
//    //below is plotting business
//    val eps11 = LApp.getRow(6).getRow(6)
//    val eps22 = LApp.getRow(7).getRow(7)
//    val eps33 = LApp.getRow(8).getRow(8)
//
//    val d31 = LApp.getRow(8).getRow(0)
//    val d32 = LApp.getRow(8).getRow(1)
//    val d33 = LApp.getRow(8).getRow(2)
//    val d24 = LApp.getRow(7).getRow(3)
//    val d15 = LApp.getRow(6).getRow(4)
//
//    val s11 = LApp.getRow(0).getRow(0)
//    val s22 = LApp.getRow(1).getRow(1)
//    val s33 = LApp.getRow(2).getRow(2)
//    val s12 = LApp.getRow(0).getRow(1)
//    val s13 = LApp.getRow(0).getRow(2)
//    val s23 = LApp.getRow(1).getRow(2)
//    val s44 = LApp.getRow(3).getRow(3)
//    val s55 = LApp.getRow(4).getRow(4)
//    val s66 = LApp.getRow(5).getRow(5)
//    plotData = List[Map[String,Double]]()
//    for(i<-0 until nbT){
//      xval = T_value.getDouble(i) * -1
//      //println("Xvalue",xval,"Yvalue",yval)
//      plotData = plotData :+ Map("a" -> xval,"b" -> eps11.getDouble(i,7)/eps11.getDouble(0,7), "c" -> 1.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> eps22.getDouble(i,7)/eps22.getDouble(0,7), "c" -> 2.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> eps33.getDouble(i,7)/eps33.getDouble(0,7), "c" -> 3.0)
//    }
//    var plot3 = Vegas("Eps graph.",width = 600.0,height = 600.0).
//      withData(plotData).
//      encodeX("a", Quantitative).
//      encodeY("b", Quantitative).
//      encodeColor("c",Nominal).
//      mark(Line)
//    plot3.show
//
//    plotData = List[Map[String,Double]]()
//    for(i<-0 until nbT){
//      xval = T_value.getDouble(i) * -1
//      //println("Xvalue",xval,"Yvalue",yval)
//      plotData = plotData :+ Map("a" -> xval,"b" -> d31.getDouble(i,7)/d31.getDouble(0,7), "c" -> 1.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> d32.getDouble(i,7)/d32.getDouble(0,7), "c" -> 2.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> d33.getDouble(i,7)/d33.getDouble(0,7), "c" -> 3.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> d24.getDouble(i,7)/d24.getDouble(0,7), "c" -> 4.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> d15.getDouble(i,7)/d15.getDouble(0,7), "c" -> 5.0)
//    }
//    var plot4 = Vegas("Eps graph.",width = 600.0,height = 600.0).
//      withData(plotData).
//      encodeX("a", Quantitative).
//      encodeY("b", Quantitative).
//      encodeColor("c",Nominal).
//      mark(Line)
//    plot4.show
//
//    plotData = List[Map[String,Double]]()
//    for(i<-0 until nbT){
//      xval = T_value.getDouble(i) * -1
//      //println("Xvalue",xval,"Yvalue",yval)
//      plotData = plotData :+ Map("a" -> xval,"b" -> s11.getDouble(i,7)/s11.getDouble(0,7), "c" -> 1.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> s22.getDouble(i,7)/s22.getDouble(0,7), "c" -> 2.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> s33.getDouble(i,7)/s33.getDouble(0,7), "c" -> 3.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> s12.getDouble(i,7)/s12.getDouble(0,7), "c" -> 4.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> s13.getDouble(i,7)/s13.getDouble(0,7), "c" -> 5.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> s23.getDouble(i,7)/s23.getDouble(0,7), "c" -> 6.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> s44.getDouble(i,7)/s44.getDouble(0,7), "c" -> 7.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> s55.getDouble(i,7)/s55.getDouble(0,7), "c" -> 8.0)
//      plotData = plotData :+ Map("a" -> xval,"b" -> s66.getDouble(i,7)/s66.getDouble(0,7), "c" -> 9.0)
//    }
//    var plot5 = Vegas("Eps graph.",width = 600.0,height = 600.0).
//      withData(plotData).
//      encodeX("a", Quantitative).
//      encodeY("b", Quantitative).
//      encodeColor("c",Nominal).
//      mark(Line)
//    plot5.show
  }
  //Below are functions to change matrix dimensions
  //They use Mandel-Voigt notation to transform symmetric tensors
  def voigt_33_61(S33:INDArray): INDArray ={
    //Takes a matrix of dimensions N x M x (3 x 3 x 3 x 3)
    //where N is the number of grains and M is the number of orientations.
    //Returns a matrix of dimensions N x M x (6 x 1)
    val transform_mat = Array(Array(0,1,2,1,2,0),
      Array(0,1,2,2,0,1),
      Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    val arr_size = S33.shape()
    var S61 = Nd4j.zeros(arr_size(0),arr_size(1),6,1)


    for(n<- 0 until arr_size(0)){
      for(m <-0 until arr_size(1)){
        for(i <- 0 until 6){
          val toPut = transform_mat(2)(i).asInstanceOf[Double]*S33.getDouble(n,m,transform_mat(0)(i).asInstanceOf[Int],transform_mat(1)(i).asInstanceOf[Int])
          S61.putScalar(Array(n,m,i,0),toPut)
        }
      }
    }
    S61
  }
  def voigt_333_36(d333:INDArray): INDArray ={
    //Takes a matrix of dimensions N x M x (3 x 3 x 3)
    //where N is the number of grains and M is the number of orientations.
    //Returns a matrix of dimensions N x M x (3 x 6)
    val transform_mat = Array(Array(0,1,2,1,2,0),
      Array(0,1,2,2,0,1),
      Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    val arr_size = d333.shape()
    var d36 = Nd4j.zeros(arr_size(0),arr_size(1),3,6)

    for(n<- 0 until (arr_size(0))){
      for(m <-0 until (arr_size(1))){
        for(i <- 0 until 3){
          for(j<-0 until 6) {
            val toPut = transform_mat(2)(j).asInstanceOf[Double] * d333.getDouble(n,m,i,transform_mat(0)(j).asInstanceOf[Int], transform_mat(1)(j).asInstanceOf[Int])
            d36.putScalar(Array(n,m,i,j), toPut)
          }
        }
      }
    }
    d36
  }
  def voigt_3333_66(C3333:INDArray): INDArray ={
    //Takes a matrix of dimensions N x M x (3 x 3 x 3 x 3)
    //where N is the number of grains and M is the number of orientations.
    //Returns a matrix of dimensions N x M x (6 x 6)
    val transform_mat = Array(Array(0,1,2,1,2,0),
      Array(0,1,2,2,0,1),
      Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    val arr_size = C3333.shape()
    var C66 = Nd4j.zeros(arr_size(0),arr_size(1),6,6)

    for(n<- 0 until arr_size(0)){
      for(m <-0 until arr_size(1)){
        for(i<-0 until 6){
          for(j<-0 until 6) {
            var toPut = transform_mat(2)(i).asInstanceOf[Double] *transform_mat(2)(j).asInstanceOf[Double] * C3333.getDouble(n,m,transform_mat(0)(i).asInstanceOf[Int], transform_mat(1)(i).asInstanceOf[Int], transform_mat(0)(j).asInstanceOf[Int],transform_mat(1)(j).asInstanceOf[Int])
            C66.putScalar(Array(n,m,i,j), toPut)
          }
        }
      }
    }
    C66
  }
  def voigt_36_333(d36:INDArray): INDArray ={
    //takes a matrix of dimensions 3 x 6
    //and returns a matrix of dimensions 3 x 3 x 3
    val transform_mat = Array(Array(0,1,2,1,2,0),
      Array(0,1,2,2,0,1),
      Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    var d333 = Nd4j.zeros(3,3,3)

    for(i<- 0 until 3){
      for(j <-0 until 6) {
        val toPut = (1 / (transform_mat(2)(j).asInstanceOf[Double]) * d36.getDouble(i, j))
        d333.putScalar(Array(i, transform_mat(0)(j).asInstanceOf[Int], transform_mat(1)(j).asInstanceOf[Int]), toPut)
        d333.putScalar(Array(i, transform_mat(1)(j).asInstanceOf[Int], transform_mat(0)(j).asInstanceOf[Int]), toPut)
      }
    }
    d333
  }
  def voigt_66_3333(C66:INDArray): INDArray ={
    //takes a matrix of dimensions 6 x 6
    //and returns a matrix of dimensions 3 x 3 x 3 x 3
    val transform_mat = Array(Array(0,1,2,1,2,0),
      Array(0,1,2,2,0,1),
      Array(1,1,1,sqrt(2),sqrt(2),sqrt(2)))
    var C3333 = Nd4j.zeros(3,3,3,3)

    for(t<- 0 until 6){
      for(u <-0 until 6) {
        val toPut = 1 / (transform_mat(2)(t).asInstanceOf[Double]*transform_mat(2)(u).asInstanceOf[Double])*C66.getDouble(t, u)
        C3333.putScalar(Array(transform_mat(0)(t).asInstanceOf[Int], transform_mat(1)(t).asInstanceOf[Int], transform_mat(0)(u).asInstanceOf[Int],transform_mat(1)(u).asInstanceOf[Int]), toPut)
        C3333.putScalar(Array(transform_mat(1)(t).asInstanceOf[Int], transform_mat(0)(t).asInstanceOf[Int], transform_mat(0)(u).asInstanceOf[Int],transform_mat(1)(u).asInstanceOf[Int]), toPut)
        C3333.putScalar(Array(transform_mat(0)(t).asInstanceOf[Int], transform_mat(1)(t).asInstanceOf[Int], transform_mat(1)(u).asInstanceOf[Int],transform_mat(0)(u).asInstanceOf[Int]), toPut)
        C3333.putScalar(Array(transform_mat(1)(t).asInstanceOf[Int], transform_mat(0)(t).asInstanceOf[Int], transform_mat(1)(u).asInstanceOf[Int],transform_mat(0)(u).asInstanceOf[Int]), toPut)
      }
    }
    C3333
  }
  //Below are rotations of matrices
  def rotate_ordre3(d333:INDArray,R:INDArray):INDArray ={
    //Takes a matrix of size N x M x (3 x 3 x 3)
    //where N is the number of grains and M is the number of orientations.
    //Rotates the matrix according to rotation matrix R and returns a new matrix
    val arr_size = d333.shape()
    var postRot = Nd4j.zeros(arr_size(0),arr_size(1),3,3,3)

    for(grain<-0 until arr_size(0)){
      for(dir<-0 until arr_size(1)){
        for(i<-0 until 3){
          for(j<-0 until 3){
            for(k<-0 until 3){
              for(m<-0 until 3){
                for(n<-0 until 3){
                  for(o<- 0 until 3){
                    var toPut = postRot.getDouble(grain,dir,i,j,k)+R.getDouble(grain,dir,i,m)*R.getDouble(grain,dir,j,n)*R.getDouble(grain,dir,k,o)*d333.getDouble(grain,dir,m,n,o)
                    postRot.putScalar(Array(grain,dir,i,j,k),toPut)
                  }
                }
              }
            }
          }
        }
      }
    }
    postRot
  }
  def rotate_ordre4(C3333:INDArray,R:INDArray):INDArray ={
    //Takes a matrix of size N x M x (3 x 3 x 3 x 3)
    //where N is the number of grains and M is the number of orientations.
    //Rotates the matrix according to rotation matrix R and returns a new matrix
    val arr_size = C3333.shape()
    var postRot = Nd4j.zeros(arr_size(0),arr_size(1),3,3,3,3)

    for(grain<-0 until arr_size(0)){
      for(dir<-0 until arr_size(1)){
        for(i<-0 until 3){
          for(j<-0 until 3){
            for(k<-0 until 3){
              for(l<-0 until 3) {
                for (m <-0 until 3) {
                  for (n <-0 until 3) {
                    for (o <-0 until 3) {
                      for (p <- 0 until 3) {
                        var toPut = postRot.getDouble(grain, dir, i, j, k, l) + (R.getDouble(grain, dir, i, m) * R.getDouble(grain, dir, j, n) * R.getDouble(grain, dir, k, o)*R.getDouble(grain, dir, l, p) * C3333.getDouble(grain, dir, m, n, o, p))
                        postRot.putScalar(Array(grain, dir, i, j, k, l), toPut)
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    postRot
  }
  def matrot(angles:INDArray): INDArray ={
    //converts a 3 x 1 orientation vector into a 3 x 3 matrix
    //UNSURE HOW DIMENSIONALITY OF ANGLES ARRAY WILL LINE UP IN DONN_POLYC BE CAREFUL
    var A_rot = Nd4j.zeros(3,3)
//    val pw = new FileWriter("DataFiles/matRotVals.txt",true)

    val c1 = cos(angles.getDouble(0))
    val s1 = sin(angles.getDouble(0))
    val c2 = cos(angles.getDouble(1))
    val s2 = sin(angles.getDouble(1))
    val c3 = cos(angles.getDouble(2))
    val s3 = sin(angles.getDouble(2))

//    pw.append(angles.toString() + "\n")
//    pw.append("c1: " + c1 + " s1: " + s1 + " c2: " + c2 + " s2: " + s2 + " c3: " + c3 + " s3: " + s3 + "\n")

    A_rot.putScalar(Array(0,0),c1*c3-s1*c2*s3)
    A_rot.putScalar(Array(0,1),-c1*s3-c2*c3*s1)
    A_rot.putScalar(Array(0,2),s1*s2)
    A_rot.putScalar(Array(1,0),c3*s1+c1*c2*s3)
    A_rot.putScalar(Array(1,1),c1*c2*c3-s1*s3)
    A_rot.putScalar(Array(1,2),-c1*s2)
    A_rot.putScalar(Array(2,0),s2*s3)
    A_rot.putScalar(Array(2,1),c3*s2)
    A_rot.putScalar(Array(2,2),c2)
//    pw.close()
    
    A_rot
  }
  //below are the eshelby functions
  def Eshelby33_muphylin(C_inf:INDArray,form_inclu:Int):INDArray = {
    var Ndemag33:INDArray = Nd4j.zeros(3,3)
    if(form_inclu == 1){
      Ndemag33 = Nd4j.eye(3)
      Ndemag33.muli(1.0/3)
    }else{
      Ndemag33.putScalar(0,0,1.0/2)
      Ndemag33.putScalar(1,1,1.0/2)
    }
    Ndemag33
  }
  def Eshelby66_muphylin(C_inf:INDArray,icalc:Int,semi_axis:Array[Array[Int]]):INDArray = {
    var SEsh66:INDArray = Nd4j.zeros(3)
    if (icalc != 1){
      println("ERROR IN ESHELBY66: ICALC INPUT WRONG")
    }else{
      var nu0:Double = C_inf.getDouble(0,1)/(C_inf.getDouble(0,0) + C_inf.getDouble(0,1))
      val ones = Nd4j.ones(3,3)
      val zeros = Nd4j.zeros(3,3)
      val eye = Nd4j.eye(3)
      var JJ = Nd4j.concat(0,Nd4j.concat(1,ones,zeros),Nd4j.concat(1,zeros,zeros)).mul(1.0/3)
      var KK = Nd4j.concat(0,Nd4j.concat(1,ones.mul(-1).add(eye.mul(3)),zeros),Nd4j.concat(1,zeros,eye.mul(3))).mul(1.0/3)
      SEsh66 = JJ.mul((1.0+nu0)/(3.0*(1-nu0))).add(KK.mul((2.0/15)*(4.0-5*nu0)/(1.0-nu0)))
    }
    SEsh66
  }
  //modified repmat function from matlab
  //given a matrix M and a desired shape, copies M to fill the desired shape
  def matRep(arr:INDArray, rows:Int, cols:Int):INDArray ={
    var inSize = arr.shape()
    inSize = cols +: inSize
    inSize = rows +: inSize
//    println("Rep size")
//    for(element <- inSize){
//      print(element)
//      print(", ")
//    }
//    println()
    var toReturn = Nd4j.zeros(inSize,'c')
    for(i<-0 until rows){
      for(j<-0 until cols){
        toReturn.getRow(i).putRow(j,arr)
      }
    }
//    print("Rep Return ")
//    print_shape(toReturn)
    toReturn
  }

  def print_shape(arr:INDArray):Unit ={
    print("Arr Size: ")
    for(element <- arr.shape()){
      print(element + " ")
    }
    println()
  }
}