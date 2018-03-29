import java.io.PrintWriter

import MME.matRep
import org.nd4j.linalg.api.ndarray.INDArray
import org.nd4j.linalg.factory.Nd4j
import org.nd4j.linalg.ops.transforms.Transforms.exp

import scala.math.log10

//Finally, the meat of the simulation, the polycrystal function
object polyfunc {
  val write = true
  def polyfunc(polyCrystal:polyc,Tmacro:INDArray, Emacro:INDArray, Tlocinit:INDArray, Elocinit:INDArray):List[INDArray] ={
    val nbg = polyCrystal.nbg

    var Smacro:INDArray = Nd4j.zeros(6,1)
    var Dmacro = Nd4j.zeros(3,1)
    var Tloc = Nd4j.zeros(nbg,1,6,1)
    var Eloc = Nd4j.zeros(nbg,1,3,1)
    var T61 = Tlocinit
    var E = Elocinit
    var S61:INDArray = null
    var D:INDArray = null
    var retHold:List[INDArray] = null

    var Tmem:INDArray = null
    var Emem:INDArray = null

    if(Tmacro.norm2Number().asInstanceOf[Double] + Emacro.norm2Number().asInstanceOf[Double]==0){
      print("Tmacro + Emacro norm is 0")
      Tloc = matRep(Tmacro,polyCrystal.nbg,1)
      Eloc = matRep(Emacro,polyCrystal.nbg,1)
      return List(Smacro, Dmacro, Tloc, Eloc, Nd4j.create(Array(0.0)))
    }
    var convergence_mec = Double.MaxValue
    var convergence_elec = Double.MaxValue
    var convergence = Double.MaxValue
    var convMecMem = Double.MaxValue
    var convElecMem = Double.MaxValue
    var diffMec:INDArray = null
    var diffElec:INDArray = null
    val alpha_num = 20
    var alpha_list = logspace(log10(0.9),-2,alpha_num)
    for (i <- 0 until alpha_num){
      alpha_list.putScalar(i, 1-alpha_list.getDouble(i))
    }
    var mec_index = 0
    var elec_index = 0
    var alpha_mec = alpha_list.getDouble(mec_index)
    var alpha_elec = alpha_list.getDouble(elec_index)
    var count = 0

    while(convergence > 1e-3){
      count += 1
      if(write){
        data_write("T61-iter-"+count,T61)
        data_write("E-iter-"+count,E)
      }
      retHold = LDC_monoc_vectorize(polyCrystal,T61,E)
      S61 = retHold(0)
      D = retHold(1)
      Smacro  = S61.mean(0)
      Dmacro = D.mean(0)

      Tmem = T61
      Emem = E
      T61 = matRep(Tmacro, nbg, 1).add(mtimesx.mul(polyCrystal.cStarMeca,matRep(Smacro,polyCrystal.nbg,1).sub(S61),'N'))
      E = matRep(Emacro, nbg, 1).add(mtimesx.mul(polyCrystal.cStarElec,matRep(Dmacro,polyCrystal.nbg,1).sub(D),'N'))

      T61 = Tmem.mul(alpha_mec).add(T61.mul(1-alpha_mec))
      E = Emem.mul(alpha_elec).add(E.mul(1-alpha_elec))

      diffMec = T61.sub(Tmem)
      diffElec = E.sub(Emem)

      convMecMem = convergence_mec
      convElecMem = convergence_elec

      convergence_mec = diffMec.norm2Number().asInstanceOf[Double] /(1 + T61.norm2Number().asInstanceOf[Double])
      convergence_elec = diffElec.norm2Number().asInstanceOf[Double] /(1 + E.norm2Number().asInstanceOf[Double])

      if(convergence_mec > convMecMem){
        //println("Diverging")
        mec_index += 1
        if(mec_index < alpha_num){
          alpha_mec = alpha_list.getDouble(mec_index)
        }
        else{
          println("ERROR: no convergence for mechanics")
          println(Tmacro)
          println(Emacro)
          count = 200
        }
      }
      if(convergence_elec > convElecMem){
        //println("Diverging")
        elec_index += 1
        if(elec_index < alpha_num){
          alpha_elec = alpha_list.getDouble(elec_index)
        }
        else{
          println("ERROR: no convergence for electricity")
          println(Tmacro)
          println(Emacro)
          count = 200
        }
      }
      convergence = Math.max(convergence_mec,convergence_elec)
      if(count>= 200){
        println("Convergence ERROR, max count reached")
        println("Tmacro",Tmacro)
        println("Emacro",Emacro)
        convergence = 0
      }
    }

    Tloc = T61
    Eloc = E
    Smacro  = S61.mean(0)
    Dmacro = D.mean(0)

    println("convergence")

    List(Smacro, Dmacro, Tloc, Eloc,Nd4j.create(Array(count.asInstanceOf[Double])))
  }
  //logspace function in matlab
  //given inputs a, b, and c, creates a logarithmically spaced vector
  //with c elements from 10^a - 10^b
  def logspace(a:Double, b:Double, c:Int): INDArray ={
    var toReturn = Nd4j.linspace(a, b, c)
    //var exponent = 1.0
    for(i<-0 until c){
      toReturn.putScalar(i,math.pow(10,toReturn.getDouble(i)))
    }
    toReturn
  }
  def LDC_monoc_vectorize(polyCrystal:polyc, T61:INDArray, E:INDArray): List[INDArray] ={
    var dpz = mtimesx.mul(polyCrystal.d36,T61,'N')
    var epspz = mtimesx.mul(polyCrystal.d36,E, 'F')

    var W = mtimesx.mul(E,polyCrystal.Pol, 'F').mul(-1).sub(mtimesx.mul(T61,polyCrystal.epsferro61, 'F')).sub(mtimesx.mul(E,dpz,'F').mul(2))

    var temp = exp(W.mul(-1).mul(polyCrystal.tetraAS)).sum(1)
    var denom = Nd4j.zeros(temp.shape()(0),6,temp.shape()(1), temp.shape()(2))
    for(i<-0 until denom.shape()(0)){
      for(j<-0 until denom.shape()(1)){
        denom.getRow(i).putRow(j,temp.getRow(i))
      }
    }
    var fvol = exp(W.mul(-1).mul(polyCrystal.tetraAS)).div(denom)
    var Dtemp = mtimesx.mul(mtimesx.mul(polyCrystal.kappa33,E,'N').add(dpz).add(polyCrystal.Pol),fvol,'N').sum(1)
    var D = Nd4j.zeros(Dtemp.shape()(0),1,Dtemp.shape()(1),Dtemp.shape()(2))
    for(i<-0 until Dtemp.shape()(0)){
      D.getRow(i).putRow(0,Dtemp.getRow(i))
    }
    var S61temp = mtimesx.mul(mtimesx.mul(polyCrystal.invC66,T61,'N').add(epspz).add(polyCrystal.epsferro61),fvol,'N').sum(1)
    var S61 = Nd4j.zeros(S61temp.shape()(0),1,S61temp.shape()(1),S61temp.shape()(2))
    for(i<-0 until S61temp.shape()(0)){
      S61.getRow(i).putRow(0,S61temp.getRow(i))
    }
    var Wtot = mtimesx.mul(fvol,W,'N').sum(1).sum(0)

    List(S61,D)
  }
  def data_write(fileName:String,toWrite:INDArray): Unit ={
    var pw = new PrintWriter("DataFiles/sim/" + fileName + ".csv")
    pw.write("grain,dir,i,j,data\n")
    pw.flush()
    for(grain<-0 until toWrite.shape()(0)){
      for(dir<-0 until toWrite.shape()(1)){
        for(i<-0 until toWrite.shape()(2)){
          for(j<-0 until toWrite.shape()(3)){
            pw.write(grain + "," + dir + "," + i + "," + j + "," + toWrite.getDouble(grain,dir,i,j) + "\n")
            pw.flush()
          }
        }
      }
    }
    pw.close()
  }
}
