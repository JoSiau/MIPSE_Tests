/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Assemblage;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;

import java.io.IOException;

/**
 * Check the errors of the Hmatrices ACA and HCAs using a matrix-vector product.
 * The results are computed with a fixed mesh and a variable error.
 *
 * @author jsiau
 */
public class CheckErrors {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/plaque";

        deb = System.nanoTime();//plaque_13532
        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + "/plaque_576.msh");
        ElementSetHomogene ES = mesh1.createHomogeneESet();
        fin = System.nanoTime();
        System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);


        Cell Cel = new Cell(ES);
        GalerkinIntegralFormulationFull IF = new GalerkinIntegralFormulationFull(Cel, Cel, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3);
        IF.assembly();
        double x[] = new double[ES.getNbElement()];
        for (int i = 0; i < x.length; i++) {
            x[i] = 1;
        }
        double res[] = new double[ES.getNbElement()];
        res = IF.getStore().product(x, res);
        ColumnVector vP0 = (new ColumnVector(res)).copy(); // La copie peut etre inutile

        NodalDeg1 alpha = new NodalDeg1(ES);
        IF = new GalerkinIntegralFormulationFull(alpha, alpha, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7);
        IF.assembly();
        x = new double[ES.getNbNoeud()];
        for (int i = 0; i < x.length; i++) {
            x[i] = 1;
        }
        res = new double[ES.getNbNoeud()];
        res = IF.getStore().product(x, res);
        ColumnVector vP1 = new ColumnVector(res);

        IF = null;

        double epsT[] = new double[5];
        for (int i = 0; i < epsT.length; i++) {
            epsT[i] = Math.pow(10, -(i + 1));
        }

        double tps[][] = new double[epsT.length][3];
        double err[][] = new double[epsT.length][3];
        double stor[][] = new double[epsT.length][3];

        ColumnVector vX1 = new ColumnVector(ES.getNbNoeud());
        vX1.setAllElements(1);
        ColumnVector vX0 = new ColumnVector(ES.getNbElement());
        vX0.setAllElements(1);
        double vP0norm = vP0.norm();
        double vP1norm = vP1.norm();
        int kmax = 50, nmin = 30;
        for (int i = 0; i < epsT.length; i++) {
            /*
            ACA DEG 0
            */
            deb = System.nanoTime();
            GalerkinIntegralFormulationACA faca = new GalerkinIntegralFormulationACA(Cel, Cel, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3,
                    epsT[i], kmax, nmin);
            faca.assembly();
            tps[i][0] = (System.nanoTime() - deb) * 1e-9;
            StorageHmatrix Haca = (StorageHmatrix) faca.getStore();
            stor[i][0] = Haca.getStorageSizePerDofs();

            ColumnVector tmp = Haca.prod(vX0);
            tmp.sub(vP0);
            err[i][0] = tmp.norm() / vP0norm;
            /*
            HCA DEG 0
            */
            deb = System.nanoTime();
            GalerkinIntegralFormulationHCA fhca = new GalerkinIntegralFormulationHCA(Cel, Cel, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3, 3,
                    epsT[i], kmax, nmin, i + 1);
            fhca.assembly();
            tps[i][1] = (System.nanoTime() - deb) * 1e-9;
            StorageHmatrix Hhca = (StorageHmatrix) fhca.getStore();
            stor[i][1] = Hhca.getStorageSizePerDofs();

            tmp = Hhca.prod(vX0);
            tmp.sub(vP0);
            err[i][1] = tmp.norm() / vP0norm;
            /*
            HCA DEG 1
            */
            deb = System.nanoTime();
            GalerkinIntegralFormulationHCA fhca1 = new GalerkinIntegralFormulationHCA(alpha, alpha, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7, 7,
                    epsT[i], kmax, nmin, i + 1);
            fhca1.assembly();
            tps[i][2] = (System.nanoTime() - deb) * 1e-9;
            StorageHmatrix Hhca1 = (StorageHmatrix) fhca1.getStore();
            stor[i][2] = Hhca1.getStorageSizePerDofs();

            tmp = Hhca1.prod(vX1);
            tmp.sub(vP1);
            err[i][2] = tmp.norm() / vP1norm;
        }

        System.out.println("\n\n");
        System.out.println("------------------------------------------------------------------------------------------------------------------------------------------");
        System.out.println(" eps \t | \t Temps d'assemblage (sec.) \t | \t\t Erreur relative obtenue \t | \t Stockage (ko/dof)");
        System.out.println("---------|------- ACA ---- HCA0 ---- HCA1 -------|--------- ACA ------- HCA0 ------- HCA1 -------|------- ACA ---- HCA0 ---- HCA1 -------|");
        for (int i = 0; i < epsT.length; i++) {
            System.out.format(" %.1e | \t %1.3f  ;  %1.3f  ; %1.3f \t | \t %1.3e  ;  %1.3e  ; %1.3e \t | \t  %1.3f  ;  %1.3f  ; %1.3f \t | \n",
                    epsT[i], tps[i][0], tps[i][1], tps[i][2], err[i][0], err[i][1], err[i][2], stor[i][0], stor[i][1], stor[i][2]);
        }
        GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
/*
13228 elmts / 6765 noeuds
------------------------------------------------------------------------------------------------------------------------------------------
 eps 	 | 	 Temps d'assemblage (sec.) 	 | 		 Erreur relative obtenue 	 | 	 Stockage (ko/dof)
---------|------- ACA ---- HCA0 ---- HCA1 -------|--------- ACA ------- HCA0 ------- HCA1 -------|------- ACA ---- HCA0 ---- HCA1 -------|
 1,0e-01 | 	 10,465  ;  6,470  ; 105,741 	 | 	 3,765e-04  ;  3,392e-04  ; 3,580e-04 	 | 	  5,321  ;  4,982  ; 3,274 	 | 
 1,0e-02 | 	 13,369  ;  7,286  ; 105,123 	 | 	 5,726e-05  ;  8,967e-05  ; 7,099e-05 	 | 	  6,769  ;  7,182  ; 4,066 	 | 
 1,0e-03 | 	 18,054  ;  8,411  ; 107,334 	 | 	 8,316e-06  ;  2,748e-05  ; 2,494e-05 	 | 	  8,715  ;  9,345  ; 4,905 	 | 
 1,0e-04 | 	 22,755  ;  10,278  ; 109,817 	 | 	 3,521e-07  ;  9,788e-07  ; 1,553e-06 	 | 	  10,758  ;  12,312  ; 5,938 	 | 
 1,0e-05 | 	 27,946  ;  13,055  ; 113,425 	 | 	 1,939e-08  ;  6,873e-08  ; 7,414e-08 	 | 	  12,910  ;  15,545  ; 7,094 	 | 
*/