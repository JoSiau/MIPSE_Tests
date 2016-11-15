/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author siau
 */
public class TestAgglo {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        ///////////////////////////////////////////////////////////////////////// 
        ////////////////////////// HMATRIX CONSTRUCTION /////////////////////////
        ///////////////////////////////////////////////////////////////////////// 
        String meshDir = new java.io.File(".").getCanonicalPath();
        /*
         meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
         ImportFlux mesh = new ImportFlux(meshDir+"/sphere/SPHERE_3538.DEC");       
         ElementSet ES  =  mesh.getRegions(0).getElementSet();      
         /*/
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/Sphere_gmsh/SPHERE_6481.msh");   
//        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + "/plaque/plaque_13532.msh");
        ElementSetHomogene ES = mesh1.createHomogeneESet();
        //*/

        ///////////////////////////////////////////////////////////////////////// 
        ///////////////////////// INITIALISATION DONNEES ////////////////////////
        ///////////////////////////////////////////////////////////////////////// 
        double deb, fin;
        int d = 3;
        int n = ES.getNbElement();
        int nbGauss = 7;
        double epsAssem = 1e-5;
        double epsAgglo[] = new double[14];
        epsAgglo[0] = 3e-3;
        epsAgglo[1] = 1e-2;
        epsAgglo[2] = 3e-2;
        epsAgglo[3] = 8e-2;
        epsAgglo[4] = 8.5e-2;
        epsAgglo[5] = 9e-2;
        epsAgglo[6] = 9.5e-2;
        epsAgglo[7] = 1e-1;
        epsAgglo[8] = 1.2e-1;
        epsAgglo[9] = 1.4e-1;
        epsAgglo[10] = 1.6e-1;
        epsAgglo[11] = 1.8e-1;
        epsAgglo[12] = 2e-1;
        epsAgglo[13] = 3e-1;
        ///////////////////////////////////////////////////////////////////////// 
        ////////////////////////// HMATRIX CONSTRUCTION /////////////////////////
        ///////////////////////////////////////////////////////////////////////// 
        Cell C = new Cell(ES);
        int kmax = 50, nmin = 30;
        GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss,
                epsAssem, kmax, nmin);
        f.assembly();
        // Hmatrix(Matrix P, HmatrixCompatible HC, double epsilon, int k, int dof, int leafSize)
        StorageHmatrix H0 = (StorageHmatrix) f.getStore();

        H0.Agglomerate(new TruncationControl("rel", epsAssem));

        int nb = 5;
        ColumnVector b[] = new ColumnVector[nb];
        double normb[] = new double[nb];
        for (int i = 0; i < nb; i++) {
            b[i] = new ColumnVector(n);
            for (int j = 0; j < n; j++) {
                b[i].setElement(j, Math.random() * Math.random() * 1000);
            }
            normb[i] = b[i].norm();
        }

        GestionnaireTaches.getGestionnaireTaches().setNbTaches(7);

        Cell Ces = new Cell((ElementSetHomogene) ES);
        GalerkinIntegralFormulation IFCell = new GalerkinIntegralFormulationFull(Ces, Ces, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss);
        IFCell.assembly();
        Matrix mEx = ((StorageFull) IFCell.getStore()).getMatrix();

        ColumnVector ex[] = new ColumnVector[nb];
        double normEx[] = new double[nb];
        for (int i = 0; i < nb; i++) {
            ex[i] = new ColumnVector(n);
            ex[i].mul(mEx, b[i]);

            normEx[i] = ex[i].norm();
        }

        int Csp[] = new int[epsAgglo.length];
        double errRel[] = new double[epsAgglo.length];
        double errAbs[] = new double[epsAgglo.length];
        ColumnVector x;
        for (int i = 0; i < epsAgglo.length; i++) {

            TruncationControl tol = new TruncationControl("rel", epsAgglo[i]);

            ///////////////////////////////////////////////////////////////////////// 
            ////////////////////////////   AGGLOMERATION   //////////////////////////
            ///////////////////////////////////////////////////////////////////////// 
            System.err.println("Go to Coarsening ! epsAgglo=" + epsAgglo[i]);

            StorageHmatrix Htmp = H0.copy(true);

            deb = System.nanoTime();
            Htmp.Agglomerate(tol);
            fin = System.nanoTime();
            System.out.println("Coarsening's over: " + (fin - deb) * 1e-9);
            Csp[i] = Htmp.getCsp();
            System.out.println("Csp= " + Csp[i]);

            double erRel = 0, erAbs = 0;
            for (int j = 0; j < nb; j++) {
                x = Htmp.prod(b[j]);
                x.sub(ex[j]);
                erAbs += x.norm();
                erRel += x.norm() / normEx[j];
            }
            errAbs[i] = erAbs / (double) nb;
            errRel[i] = erRel / (double) nb;
        }

        PrintWriter fic = new PrintWriter(new BufferedWriter(new FileWriter("Agglo.out")));
        for (int i = 0; i < epsAgglo.length; i++) {
            System.out.format("i= %2d  \t  eps= %1.3e  \t  errAbs= %1.3e  \t  errRel= %1.3e  \t  Csp=%d \n", i, epsAgglo[i], errAbs[i], errRel[i], Csp[i]);
            fic.println(i + " " + epsAgglo[i] + " " + errRel[i] + " " + Csp[i]);
        }
        fic.close();
    }

}
