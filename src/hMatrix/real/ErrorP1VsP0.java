/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;

import java.io.IOException;

/**
 * @author siau
 */
public class ErrorP1VsP0 {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        System.out.println("ErrorP1VsP0.java");
        double deb, fin;

        String meshDir = new java.io.File(".").getCanonicalPath();
        int kmax = 50, nmin = 30;
        //*
        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
//        ImportFlux mesh = new ImportFlux(meshDir + "/sphere/SPHERE_1884.DEC");
        ImportFlux mesh = new ImportFlux(meshDir + "/PLAQUE2000.DEC");
//        ImportFlux mesh = new ImportFlux("VariateurDEC/ATV71_105963.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();
        /*/
         meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
         //         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/plaque/plaque_4548.msh");   
         //         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/cube.msh");   
         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/sphere_gmsh/sphere_6481.msh");                    
         ElementSetHomogene ES = mesh1.createHomogeneESet();
         //*/
        int d = 3;
        int nbNodes = ES.getNbNoeud();
        int nbElmts = ES.getNbElement();

        System.out.println("nbNoeuds= " + ES.getNbNoeud());
        System.out.println("nbElmts= " + ES.getNbElement());

        System.out.print("Computing the assembled matrix...");
        NodalDeg1 Ces = new NodalDeg1(ES);
        System.out.print("1");
        GalerkinIntegralFormulation g = new GalerkinIntegralFormulationFull(Ces, Ces, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3);
        System.out.print("2");
        g.assembly();
        System.out.println("Done !");
        ColumnVector v = new ColumnVector(nbNodes);
        v.setAllElements(1);
//        for (int i = 0; i < nbNodes; i++) {
//            v.setElement(i, Math.random() * 100);
//        }
        System.out.println("Doing the PMV");
        ColumnVector exV = new ColumnVector(nbNodes);
        exV.mul(((StorageFull) g.getStore()).getMatrix(), v);
        double exN = exV.norm();
        /*
         * HMATRIX P1
         */
        System.out.println("Construction of the H-matrix P1");
        int epsPow = 5;
        double eps = Math.pow(10, -epsPow);
        NodalDeg1 alpha = new NodalDeg1(ES);
        deb = System.nanoTime();
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(alpha, alpha, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7, 7, eps, kmax, nmin, epsPow - 1);
        f.assembly();
        StorageHmatrix H = (StorageHmatrix) f.getStore();
        fin = System.nanoTime();
        double tHD1 = (fin - deb) * 1e-9;
        double stor1 = H.getStorageSizePerDofs();

        ColumnVector vh = H.prod(v);
        vh.sub(exV);
        double errD1 = vh.norm() / exN;
        System.out.println("P1: eps= " + Math.pow(10, -epsPow) + " \t Error rel= " + errD1);

        /*
         * *********************************************************************
         */
        System.out.print("Computing the assembled matrix...");
        Cell Ce = new Cell(ES);
        GalerkinIntegralFormulation g0 = new GalerkinIntegralFormulationFull(Ce, Ce, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7);
        g0.assembly();
        System.out.println("Done !");
        ColumnVector v0 = new ColumnVector(nbElmts);
        v0.setAllElements(1);
//        for (int i = 0; i < nbNodes; i++) {
//            v.setElement(i, Math.random() * 100);
//        }
        System.out.println("Doing the PMV");
        ColumnVector exV0 = new ColumnVector(nbElmts);
        exV0.mul(((StorageFull) g0.getStore()).getMatrix(), v0);
        double exN0 = exV0.norm();

        /*
         * HMATRIX P0
         */
        System.out.println("Construction of the H-matrix P0");
        Cell C = new Cell(ES);
        deb = System.nanoTime();
        GalerkinIntegralFormulationHCA f1 = new GalerkinIntegralFormulationHCA(C, C, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7, 7, eps, kmax, nmin, epsPow - 1);
//        Hmatrix Hhca = new Hmatrix(X,f,1e-3,kmax,n,nmin,2);
        f1.assembly();

        StorageHmatrix H0 = (StorageHmatrix) f1.getStore();
        fin = System.nanoTime();
        double tHD0 = (fin - deb) * 1e-9;
        double stor0 = H0.getStorageSizePerDofs();

        ColumnVector vh0 = H0.prod(v0);
        vh0.sub(exV0);
        double errD0 = vh0.norm() / exN0;
        System.out.println("P0: eps= " + Math.pow(10, -epsPow) + " \t Error rel= " + errD0);

        System.err.println("Time :  P1= " + tHD1 + " \t P0= " + tHD0 + "\n"
                + "Storage: P1= " + stor1 + " ko/dof \t P0= " + stor0 + " ko/dof");
    }

}
/*
PLAQUE 4370 / 2273 
Time to assembly the Hmatrix:  P1= 51.661786955000004 	 P0= 7.951330961000001
Storage: P1= 2.936216390160183 ko/dof 	 P0= 5.8243474685354695 ko/dof

CUBE 6238 / 3349
Time :  P1= 120.90706513900001 	 P0= 54.859102262
Storage: P1= 5.671085985091375 ko/dof 	 P0= 10.036232065565887 ko/dof

* SPHERE 6198 / 3101
Time :  P1= 119.02303390600001 	 P0= 85.144846806
Storage: P1= 4.505742779929009 ko/dof 	 P0= 8.737587982010325 ko/dof
 */
