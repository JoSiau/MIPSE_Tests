/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package hMatrix.real.Complexite;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.CurlSEdgeDeg1;
import g2elab.mipse.meshCore.functionSpace.GradNodalDeg1;
import g2elab.mipse.meshCore.functionSpace.Hgrad;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.CrossDGmulAlpha;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;

import java.io.IOException;

/**
 * Verifie la complexite de lassemblage HCA sur la formulation magnetodynamque.
 *
 * @author jsiau
 */
public class HCAvect {

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

        int nF = nFile.length;

        double tHmat[] = new double[nF];
        double nNodes[] = new double[nF];
        double nElmts[] = new double[nF];

        int nbGauss = 3;
        for (int k = 0; k < nF; k++) {
            double dd = System.nanoTime();
            System.gc();
            System.out.println("Time gc= " + (System.nanoTime() - dd) * 1e-9);

            ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + nFile[k] + ".msh");
            ElementSetHomogene mesh = mesh1.createHomogeneESet();

            int d = 3;
            int n = mesh.getNbNoeud();
            nNodes[k] = n;
            nElmts[k] = mesh.getNbElement();

            System.out.println("nbNoeuds= " + mesh.getNbNoeud() + "\t nbElmts= " + mesh.getNbElement());
            /*
             * CREATION DES ESPACES FONCTIONELS
             */
            NodalDeg1 alpha = new NodalDeg1(mesh);
            Hgrad alphaN = alpha.createProjOnNormal(1);
            GradNodalDeg1 gradalpha = new GradNodalDeg1(mesh);
            CurlSEdgeDeg1 rotSalpha = new CurlSEdgeDeg1(alpha);
            /*
             * DEFINIE LA FONCTION (LE NOYAU) QUE L'ON DESIRE
             */
            deb = System.nanoTime();
            GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(alphaN, rotSalpha, new CrossDGmulAlpha(1), new SelfElementFixedGauss(3, new Cancel()), nbGauss, nbGauss,
                    1e-5, 50, 30, 4);
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
