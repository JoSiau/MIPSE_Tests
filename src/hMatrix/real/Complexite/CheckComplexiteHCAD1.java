/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package hMatrix.real.Complexite;


import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;

import java.io.IOException;

/**
 * Check the complexity of the HCA Degree 1
 *
 * @author jsiau
 */
public class CheckComplexiteHCAD1 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/Sphere_gmsh/sphere_";

        int nFile[] = new int[5];
        nFile[0] = 6481;
        nFile[1] = 12529;
        nFile[2] = 34449;
        nFile[3] = 62713;
        nFile[4] = 132241;
//        nFile[5] = 301791;

        int nF = nFile.length;

        double tHmat[] = new double[nF];
        double nNodes[] = new double[nF];
        double nElmts[] = new double[nF];

        for (int k = 0; k < 3; k++) {
            // Time counters and time-variables-templates
            double beg = 0, end = 0;

            ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + nFile[k] + ".msh");
            ElementSetHomogene ES = mesh1.createHomogeneESet();

            int d = 3;
            int n = ES.getNbNoeud();
            nNodes[k] = n;
            nElmts[k] = ES.getNbElement();

            System.out.println("nbNoeuds= " + ES.getNbNoeud() + "\t nbElmts= " + ES.getNbElement());
            /*
             * DEFINIE LA FONCTION (LE NOYAU) QUE L'ON DESIRE
             */
            NodalDeg1 alpha = new NodalDeg1(ES);
            deb = System.nanoTime();
            GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(alpha, alpha, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3, 3,
                    1e-3, 50, 30, 4);
            System.out.println("Calculs de quadratures = " + (System.nanoTime() - deb) * 1e-9);
            f.assembly();
            fin = System.nanoTime();
            tHmat[k] = (fin - deb) / 1e9;
            System.err.println("Time to assembly the Hmatrix = " + tHmat[k]);
        }
        System.out.println("\n\n");
        for (int i = 0; i < nF; i++) {
            System.out.println(nElmts[i] + " , " + nNodes[i] + " , " + tHmat[i]);
        }
    }
}
/*
6198.0 , 3101.0 , 26.355142485
12138.0 , 6071.0 , 49.755008912
33806.0 , 16905.0 , 158.620624602
61842.0 , 30923.0 , 320.032250364
130974.0 , 65489.0 , 745.55527317
*/