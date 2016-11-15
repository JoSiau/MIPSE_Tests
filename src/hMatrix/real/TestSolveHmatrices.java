/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
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
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import g2elab.mipse.numericalTools.preconditioner.PrecondIdentityReal;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;

import java.io.IOException;

/**
 * @author jsiau
 */
public class TestSolveHmatrices {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/plaque";

        int[] nElmt = new int[6];
        nElmt[0] = 1568;
        nElmt[1] = 6152;
        nElmt[2] = 18928;
        nElmt[3] = 55476;
        nElmt[4] = 98312;
        nElmt[5] = 185696;

        double tpsAss[][] = new double[nElmt.length][3];
        double stor[][] = new double[nElmt.length][3];
        double tpsSolv[][] = new double[nElmt.length][3];
        int nIt[][] = new int[nElmt.length][3];

        int kmax = 50, nmin = 30, ordre = 3;
        double eps = 1e-4;

        PrecondIdentityReal precond = new PrecondIdentityReal();
        double configRes[] = {100, 1e-6, 1, -50};
        for (int i = 5; i < nElmt.length; i++) {
            deb = System.nanoTime();//plaque_13532
            ImportFlux mesh1 = new ImportFlux("D:/Meshs/Cube/CUBE_" + nElmt[i] + ".DEC");
            ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
            fin = System.nanoTime();
            System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);
            int nElt = ES.getNbElement(), nNod = ES.getNbNoeud();

            Cell Cel = new Cell(ES);
            NodalDeg1 alpha = new NodalDeg1(ES);
            /*
             ACA DEG 0
             */
            deb = System.nanoTime();
            GalerkinIntegralFormulationACA faca = new GalerkinIntegralFormulationACA(Cel, Cel, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7,
                    eps, kmax, nmin);
            faca.assembly();
            tpsAss[i][0] = (System.nanoTime() - deb) * 1e-9;
            StorageHmatrix Haca = (StorageHmatrix) faca.getStore();
            stor[i][0] = Haca.getStorageSizePerDofs();
            // RESOLUTION
            deb = System.nanoTime();
            FGMResReal solver = new FGMResReal(Haca, precond);

            solver.setInfoResolution(configRes);
            double xFMM[] = new double[nElt];
            double secondMembre[] = new double[nElt];
            for (int j = 0; j < nElt; j++) {
                secondMembre[i] = 1;
            }
            xFMM = solver.solve(xFMM, secondMembre);
            fin = System.nanoTime();
            tpsSolv[i][0] = (fin - deb) / 1e9;
            System.err.println("Time to solve the FMM = " + (fin - deb) / 1e9);
            nIt[i][0] = solver.nbIter;

            /*
             HCA DEG 0
             */
            deb = System.nanoTime();
            GalerkinIntegralFormulationHCA fhca = new GalerkinIntegralFormulationHCA(Cel, Cel, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7, 3,
                    eps, kmax, nmin, ordre);
            fhca.assembly();
            tpsAss[i][1] = (System.nanoTime() - deb) * 1e-9;
            StorageHmatrix Hhca = (StorageHmatrix) fhca.getStore();
            stor[i][1] = Hhca.getStorageSizePerDofs();
            // RESOLUTION
            deb = System.nanoTime();
            solver = new FGMResReal(Hhca, precond);

            solver.setInfoResolution(configRes);
            xFMM = solver.solve(xFMM, secondMembre);
            fin = System.nanoTime();
            tpsSolv[i][1] = (fin - deb) / 1e9;
            System.err.println("Time to solve the FMM = " + (fin - deb) / 1e9);
            nIt[i][1] = solver.nbIter;

            /*
             HCA DEG 1
             */
            deb = System.nanoTime();
            GalerkinIntegralFormulationHCA fhca1 = new GalerkinIntegralFormulationHCA(alpha, alpha, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7, 7,
                    eps, kmax, nmin, ordre);
            fhca1.assembly();
            tpsAss[i][2] = (System.nanoTime() - deb) * 1e-9;
            StorageHmatrix Hhca1 = (StorageHmatrix) fhca1.getStore();
            stor[i][2] = Hhca1.getStorageSizePerDofs();
            // RESOLUTION
            deb = System.nanoTime();
            solver = new FGMResReal(Hhca1, precond);
            solver.setInfoResolution(configRes);
            xFMM = new double[nNod];
            secondMembre = new double[nNod];
            for (int j = 0; j < nNod; j++) {
                secondMembre[i] = 1;
            }
            xFMM = solver.solve(xFMM, secondMembre);
            fin = System.nanoTime();
            tpsSolv[i][2] = (fin - deb) / 1e9;
            System.err.println("Time to solve the FMM = " + (fin - deb) / 1e9);
            nIt[i][2] = solver.nbIter;
        }

        System.out.println("\n\n");
        System.out.println("------------------------------------------------------------------------------------------------------------------------------------------");
        System.out.println(" nElmt \t | \t Temps d'assemblage (sec.) \t | \t\t Temps de resolution \t | \t Nombre iteration");
        System.out.println("---------|------- ACA ---- HCA0 ---- HCA1 -------|--------- ACA ------- HCA0 ------- HCA1 -------|------- ACA ---- HCA0 ---- HCA1 -------|");
        for (int i = 0; i < nElmt.length; i++) {
            System.out.format(" %d | \t %1.3e  ;  %1.3e  ; %1.3e \t | \t %1.3e  ;  %1.3e  ; %1.3e \t | \t  %d  ;  %d  ; %d \t | \n",
                    nElmt[i], tpsAss[i][0], tpsAss[i][1], tpsAss[i][2], tpsSolv[i][0], tpsSolv[i][1], tpsSolv[i][2], nIt[i][0], nIt[i][1], nIt[i][2]);
        }
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
/*
 ------------------------------------------------------------------------------------------------------------------------------------------
 nElmt 	 | 	 Temps d'assemblage (sec.)       	 | 		 Temps de resolution             | 	 Nombre iteration
 ---------|--------- ACA ------- HCA0 ------- HCA1 -------|--------- ACA ------- HCA0 ------- HCA1 -------|------ ACA -- HCA0 -- HCA1 ----|
 1568    | 	 2,533e+00  ;  2,051e+00  ; 1,919e+01 	 | 	 2,604e-01  ;  1,582e-01  ; 1,071e-01 	 | 	  37  ;  16  ; 48 	 | 
 6152    | 	 1,266e+01  ;  7,429e+00  ; 6,513e+01 	 | 	 1,804e+00  ;  8,958e-01  ; 1,456e+00 	 | 	  50  ;  17  ; 62 	 | 
 18928   | 	 5,235e+01  ;  2,637e+01  ; 2,013e+02 	 | 	 8,777e+00  ;  4,486e+00  ; 6,700e+00 	 | 	  59  ;  22  ; 77 	 | 
 55476   | 	 2,122e+02  ;  1,068e+02  ; 5,732e+02 	 | 	 4,908e+01  ;  2,395e+01  ; 3,556e+01 	 | 	  75  ;  30  ; 96 	 | 
 98312   | 	 4,401e+02  ;  2,051e+02  ; 1,144e+03 	 | 	 1,021e+02  ;  3,065e+01  ; 7,823e+01 	 | 	  82  ;  19  ; 108 	 | 
 185696 | 	 1,013e+03  ;  4,336e+02  ; 2,033e+03 	 | 	 2,504e+02  ;  7,417e+01  ; 1,838e+02 	 | 	  94  ;  23  ; 123 	 | 
 */
