/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;

import java.io.IOException;

/**
 * FILE TO TEST THE CONSTRUCTION OF A BINARY TREE WITH P_1 GEOMETRY
 *
 * @author siau
 */
public class TestP1 {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        String meshDir = new java.io.File(".").getCanonicalPath();
        /*
         meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
         //        ImportFlux mesh = new ImportFlux(meshDir + "/sphere/SPHERE_1884.DEC");
         ImportFlux mesh = new ImportFlux(meshDir + "/PLAQUE2000.DEC");
         //        ImportFlux mesh = new ImportFlux("VariateurDEC/ATV71_105963.DEC");
         ElementSet ES = mesh.getRegions(0).getElementSet();
         /*/
        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
//         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/plaque/plaque_4548.msh");   
//         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/cube.msh");   
        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + "/sphere/SPHERE_4893.msh");
        ElementSetHomogene ES = mesh1.createHomogeneESet();
        //*/
        int d = 3;
        int n = ES.getNbNoeud();

        System.out.println("nbNoeuds= " + ES.getNbNoeud());
        System.out.println("nbElmts= " + ES.getNbElement());

//        System.out.print("Computing the assembled matrix...");
//        NodalDeg1 Ces = new NodalDeg1((ElementSetHomogene) ES);
//        GalerkinIntegralFormulation g = new GalerkinIntegralFormulation(Ces, new MultG(), new SimpleTriangleLinearChargePotential());
//        g.matAssembly(Ces, 7, 7);
//        System.out.println("Done !");
//        
//        System.out.println("Checking the error...");
//        ColumnVector v = new ColumnVector(n);
//        for (int i = 0; i < n; i++) {
//            v.setElement(i, Math.random() * 100);
//        }
//        System.out.println("Doing the PMV");
//        ColumnVector exV = new ColumnVector(n);
//        exV.mul(((StorageFull) g.getStore()).getMatrix(), v);
//        double exN = exV.norm();
        int epsPow = 5;
//for( int epsPow = 2; epsPow<7;epsPow++){
        NodalDeg1 C = new NodalDeg1(ES);
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, C, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7, 7,
                Math.pow(10, -epsPow), 50, 30, epsPow - 1);
        //public Hmatrix(ElementSet ES, HmatrixCompatibleHCADeg1 HC, double epsilon, int k, int ordre, int leafSize, double eta)
        StorageHmatrix H = new StorageHmatrix(f);
        //*/

        /*
         Matrix ex = ((StorageFull)g.getStore()).getMatrix().copy();
        
         ex.permuteColumns(H.getIdx2Dof(), new int[n]);
         ex.permuteRows(H.getIdx2Dof(), new int[n]);
        
        
         H.CheckError(ex, Math.pow(10, -epsPow));
         //*/
//        ColumnVector vh = H.prodIndex(v);
//        vh.sub(exV);
//        System.out.println("eps= " + Math.pow(10, -epsPow)+" \t Error rel= " + vh.norm() / exN);
//        //*/
//        
////        System.out.println("Norme Ex=" + exN);
//        System.out.println("Done !");
//}
    }

}
/*
PLAQUE_2000.DEC
eps= 0.001 	 Error rel= 7.592002432093733E-5
eps= 1.0E-4 	 Error rel= 1.8059727983466773E-5
eps= 1.0E-5 	 Error rel= 8.341981070743761E-7
eps= 1.0E-6 	 Error rel= 1.315054964352821E-7
* 
* 
* 
CUBE_3349.msh
eps= 0.01 	 Error rel= 0.004533009263961797
eps= 0.001 	 Error rel= 3.6706934313548626E-4
eps= 1.0E-4 	 Error rel= 7.322887346056327E-5
eps= 1.0E-5 	 Error rel= 3.234623658664096E-6
eps= 1.0E-6 	 Error rel= 3.0300733731066785E-7
* 
* 
* 
SPHERE_1501.msh
eps= 0.01 	 Error rel= 0.0016015681453576194
eps= 0.001 	 Error rel= 1.1296976281934484E-4
eps= 1.0E-4 	 Error rel= 1.1987219655358486E-4
eps= 1.0E-5	 Error rel= 2.142611811717322E-6
eps= 1.0E-6 	 Error rel= 6.979047311338114E-8
*/
