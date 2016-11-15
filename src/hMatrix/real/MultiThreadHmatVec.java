/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;

import java.io.IOException;

/**
 * @author siau
 */
public class MultiThreadHmatVec {

    /**
     * @param args
     */
    public static void main(String[] args) {
        double deb, fin, eps = 1e-3;
        String dir = null;
        try {
            dir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
        }
        dir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
////        // 8, 64, 104, 248, 572, 1884, 2104, 3538
        //ImportFlux mesh = new ImportFlux("D:/WORK/ForgeG2ELAB/mesh/sphere/SPHERE_3538.DEC");
        //*
        ImportFlux mesh = new ImportFlux(dir + "/sphere/SPHERE_1884.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();
        /*/
         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,dir+"/sphere_gmsh/sphere_24091.msh");
         ElementSetHomogene ES = mesh1.createHomogeneESet();
         //*/
        int d = 3;
        int n = ES.getNbElement();

        /*
         * DEFINIE LA FONCTION (LE NOYAU) QUE L'ON DESIRE
         */
        Cell C = new Cell(ES);
        int kmax = 50, nmin = 30;
        GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3,
                eps, kmax, nmin);
        f.assembly();
        TruncationControl tol = new TruncationControl("rel", eps);
//        TruncationControl tol = new TruncationControl(0,eps);
        tol.setTruncationMethod(2);
        tol.printErrors();
        //*
        StorageHmatrix Hm0 = (StorageHmatrix) f.getStore();
        /*/
         BlockClusterTree Bc0 = new BlockClusterTree(n,8,d);
         Hmatrix Hm0 = new Hmatrix(Bc0, f, n, eps, kmax);
         //*/
        ColumnVector v = new ColumnVector(n);
        for (int i = 0; i < n; i++) {
            v.setElement(i, Math.random());
        }

        deb = System.nanoTime();
        ColumnVector x1 = Hm0.prod(v);
        fin = System.nanoTime();
        double t1 = (fin - deb) * 1e-9;
        System.out.println("Time Hmat-Vect normal: " + t1);

        GestionnaireTaches.getGestionnaireTaches().getNbTachesParalleles();
        deb = System.nanoTime();
        ColumnVector x2 = Hm0.prodMT(v);
        fin = System.nanoTime();
        double t2 = (fin - deb) * 1e-9;
        System.out.println("Time Hmat-Vect MT'd: " + t2);

        System.out.println("Ratio: " + t1 / t2);

        x2.sub(x1);
        System.out.println("Erreur relative= " + x2.norm() / x1.norm());
        GestionnaireTaches.getGestionnaireTaches().stop();

    }
}
