/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Arithmetics;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author siau
 */
public class CheckComplexiteHmatVect {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        double deb, fin, eps = 1e-2;
        String meshDir = new java.io.File(".").getCanonicalPath();

        int I[] = new int[11];
        I[0] = 6481;
        I[1] = 12529;
        I[2] = 19317;
        I[3] = 24091;
        I[4] = 34449;
        I[5] = 47431;
        I[6] = 62713;
        I[7] = 75631;
        I[8] = 94397;
        I[9] = 132241;
        I[10] = 301791;
//                I[0] = 4881;

        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere_gmsh";

        PrintWriter fic = new PrintWriter(new BufferedWriter(new FileWriter("TpsPMV.out")));

        StorageHmatrix Hm0;
        ColumnVector v;
        Matrix X;

        // Can go from -2 up to 10, ask me for the meshes or make your own !
        for (int j = 0; j <= 9; j++) {
            deb = System.nanoTime();
            ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + "/sphere_" + I[j] + ".msh");
            ElementSetHomogene ES = mesh1.createHomogeneESet();
            fin = System.nanoTime();
            System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);

            int d = 3;
            int n = ES.getNbElement();
            System.out.println(" n = " + n + "\t d = " + d);
            X = new Matrix(n, d);
            //*
            for (int i = 0; i < n; i++) {
                X.setRow(i, ES.getElements(i).getCentroid());
            }

            /*
             * DEFINIE LA FONCTION (LE NOYAU) QUE L'ON DESIRE
             */
            Cell FS = new Cell(ES);
//            f1 f = new f1(X);
            int nmin = 30;
            //*
            int kmax = 50;
            deb = System.nanoTime();
            GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(FS, FS, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7, eps, kmax, nmin);
            f.assembly();
            Hm0 = (StorageHmatrix) f.getStore();
            fin = System.nanoTime();

            System.out.println("Time to assembly the Hmatrix = " + (fin - deb) / 1e9);
            // System.out.println("Sparsity constant= "+Hm0.getCsp()+" \t Depth= "+Hm0.getDepth());
            //System.out.println("Go to Prod ! ");

            PrintWriter fic1 = new PrintWriter(new BufferedWriter(new FileWriter("TpsPMV" + n + ".out")));

            v = new ColumnVector(n);
            for (int i = 0; i < n; i++) {
                v.setElement(i, Math.random() * 1000);
            }

            int nbMoy = 20;
            double tps = 0;
            for (int k = 0; k < nbMoy; k++) {
                //System.out.println("i="+k+"\t tps="+tps/k);

                deb = System.nanoTime();
                ColumnVector u = Hm0.prodWithoutRenum(v);
                fin = System.nanoTime();
                tps += (fin - deb) * 1e-9;
                fic1.println(k + "  ,  " + tps / (k + 1));
            }
            fic1.close();
            tps /= nbMoy;
            System.out.println("Time to compute the 1st PMV= " + tps);
            fic.print(n + "  ,  " + tps);

            deb = System.nanoTime();
            Hm0.Agglomerate(new TruncationControl("rel", eps));
            fin = System.nanoTime();
            System.out.println("Time to coarsen= " + (fin - deb) / 1e9);

            tps = 0;
            for (int k = 0; k < nbMoy; k++) {
                //System.out.println("i="+k+"\t tps="+tps/k);
                deb = System.nanoTime();
                ColumnVector u = Hm0.prodWithoutRenum(v);
                fin = System.nanoTime();
                tps += (fin - deb) * 1e-9;
                fic1.println(k + "  ,  " + tps / (k + 1));
            }
            tps /= nbMoy;
            System.out.println("Time to compute the 2nd PMV= " + tps);
            fic.println("  ,  " + tps);
            System.out.println("++++++++++++++++++++++++++++++++++++++++++");
        }
        fic.close();
    }
}
