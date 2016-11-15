/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Complexite;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.Element;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.IOException;

/**
 * COMPUTE THE ASYMPTOTICAL BEHAVIOUR OF THE HMATRIX HCA DEGREE 0 AND 1.
 *
 * @author siau
 */
public class CheckComplexiteHCAs {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();

        /* SPHERES
         meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere_gmsh/sphere_"";
         int I[] = new int[11];
         I[0] = 6481;I[1] = 12529;I[2] = 19317;I[3] = 24091;I[4] = 34449;I[5] = 47431;
         I[6] = 62713;I[7] = 75631;I[8] = 94397;I[9] = 132241;I[10] = 301791;
         /*/ // PLAQUES
        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/plaque/plaque_";
        int I[] = new int[5];
        I[0] = 4548;
        I[1] = 13532;
        I[2] = 24228;
        I[3] = 53300;
        I[4] = 95006;
        // I[5] = 47431;I[6] = 62713;I[7] = 75631;I[8] = 94397;I[9] = 132241;I[10] = 301791;
        //*/

        double timeD1[] = new double[I.length];
        double timeD0[] = new double[I.length];

        double storD1[] = new double[I.length];
        double storD0[] = new double[I.length];

        int nEl[] = new int[I.length];
        int nNd[] = new int[I.length];

        ColumnVector v;

        for (int j = 0; j < I.length; j++) {
            deb = System.nanoTime();
            ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + I[j] + ".msh");
            ElementSetHomogene ES = mesh1.createHomogeneESet();
            fin = System.nanoTime();
            System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);

            nEl[j] = ES.getNbElement();
            nNd[j] = ES.getNbNoeud();

            System.out.println("NbNoeuds = " + ES.getNbNoeud() + " \t NbElmts = " + ES.getNbElement());


            /*
             * DEFINIE LA FONCTION (LE NOYAU) QUE L'ON DESIRE
             */
            NodalDeg1 alpha = new NodalDeg1(ES);
            deb = System.nanoTime();
            GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(alpha, alpha, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3, 3,
                    1e-5, 50, 30, 4);
            f.assembly();
            StorageHmatrix H = (StorageHmatrix) f.getStore();
            fin = System.nanoTime();
            timeD1[j] = (fin - deb) / 1e9;
            storD1[j] = H.getStorageSizePerDofs();
            System.out.println("Time to assembly the Hmatrix = " + timeD1[j]);
            System.err.println(ES.getNbNoeud() + " , " + timeD1[j]);

            /*
             * 
             * 
             * 
             */
            System.gc();

            Cell C = new Cell(ES);
            deb = System.nanoTime();
            GalerkinIntegralFormulationHCA f1 = new GalerkinIntegralFormulationHCA(C, C, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3, 3,
                    1e-5, 50, 30, 4);
            f1.assembly();
            Matrix X = new Matrix(ES.getNbElement(), 3);
            Element E;
            for (int i = 0; i < X.getRowCount(); i++) {
                // Load the next element
                E = ES.getElements(i);
                // Get the gravity center and Copy it in X !
                X.setRow(i, E.getCentroid());
            }
            StorageHmatrix H0 = (StorageHmatrix) f1.getStore();
            fin = System.nanoTime();
            timeD0[j] = (fin - deb) / 1e9;
            storD0[j] = H0.getStorageSizePerDofs();
            System.err.println(ES.getNbElement() + " , " + timeD0[j]);

            System.gc();

        }

        System.out.println("nEl , nNd , timeD0 , timeD1 , storeD0 , storeD1");
        for (int i = 0; i < I.length; i++) {
            System.out.format("%d , %d , %1.3e , %1.3e , %1.3e , %1.3e \n", nEl[i], nNd[i], timeD0[i], timeD1[i], storD0[i], storD1[i]);
        }

    }
}

/* SPHERE
3101 , 125.388991568
6071 , 325.066171874
9417 , 597.175728155
11774 , 713.296088615
16905 , 1044.623446559
23336 , 1565.625701391
30923 , 1892.373409358
37340 , 1831.316555837
46669 , 2269.93305047
65489 , 3212.937038621
* 
PLAQUE
nEl , nNd , timeD0 , timeD1 , storeD0 , storeD1
4370 , 2273 , 4,543e+00 , 1,554e+01 , 5,824e+00 , 2,936e+00 
13228 , 6765 , 1,385e+01 , 5,167e+01 , 7,838e+00 , 4,125e+00 
23822 , 12113 , 2,838e+01 , 9,678e+01 , 8,916e+00 , 4,757e+00 
52696 , 26649 , 8,132e+01 , 2,358e+02 , 1,048e+01 , 5,494e+00 
94200 , 47502 , 2,643e+02 , 4,602e+02 , 1,158e+01 , 6,149e+00 
*/
