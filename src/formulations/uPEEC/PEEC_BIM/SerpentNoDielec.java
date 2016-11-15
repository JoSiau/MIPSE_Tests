/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.PEEC_BIM;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class SerpentNoDielec {

    /**
     * Serpentin sans bornes.
     *
     * @param args
     */
    public static void main(String[] args) {

//        System.out.println("Entrez une frequence:");
//        Scanner sc = new Scanner(System.in);
//        double f[] = {1e6};
        double f[] = getF();

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Zf_Serpent.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Surf/S3_NB_FS.DEC");
        /*
         *******************************
         ***  Import du fichier Flux ***
         *******************************
         Nombre de regions importees : 3
         Region 0, Nom : CU, type : 2, 1400 elements
         Region 1, Nom : FR4, type : 2, 358 elements
         Region 2, Nom : NULL, type : 1, 716 elements
         *******************************
         ***       Fin Import        ***
         *******************************
         */
        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((SurfaceRegion) mesh.getRegion(0), 1 / 1.68e-8, 35e-6, null, null);

        solP.setPtsDeGaussInductifs(64, 25, 25);
        solP.setPtsDeGaussCapacitifs(9, 9);
//        solP.setPtsDeGaussInductifs(25, 16, 16);
//        solP.setPtsDeGaussCapacitifs(9, 9);
        solP.setAnalyticalIntegrationCapa(true);
        solP.setCorrectionVoisin(false);

//        solP.setEquiPot(new int[]{512, 1212});
        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int indNoeudBord = fd.getNbElement();
        int nbBranches = solP.getNbLignes();
        circuitPur.addSourceISimple(nbBranches, 404, 809, "Source", 1.0, 0.0);
        circuitPur.addSourceUSimple(1000001, 297, 702, "cct", 0, 0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

//        solP.setExportTopologie("D:");
//        solP.integration();
//        solP.checkMatricesSingularities();
//        solP.getZ(10);
        solP.setHmatrixCompression(false);

        double z[] = new double[f.length];
        for (int iF = 0; iF < f.length; iF++) {

            double[][] ib = solP.resolutionDirectePureConductor(f[iF]);
//        double[][] ib = solP.resolutionDirecte(f);

            System.err.println("Memory to store the matrices = " + (solP.getMemoryUsed() / 1024) + " ko");

//        double[][] ib = solP.resolutionIterativePureConductor(f);
            double I[] = new double[]{ib[0][2 * nbBranches], ib[0][2 * nbBranches + 1]};
            double U[] = new double[]{ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]};
            System.out.println("I= " + Arrays.toString(I));
            System.out.println("U= " + Arrays.toString(U));
            z[iF] = Math.hypot(U[0], U[1]) / Math.hypot(I[0], I[1]);
            System.out.println("|Z|= " + z[iF]);

            Matrix res = new Matrix(2, solP.getNbLignes());
            for (int i = 0; i < res.getColumnCount(); i++) {
                res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i] / 35e-6);
                res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1] / 35e-6);
            }

            if (iF % 10 == 0) {
                PEEC_FULLY_SURF.exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/SERPENT/f_" + f[iF], fd, res, solP.getDielecCell(), solP.getQ(), true);
            }
        }

        for (int i = 0; i < f.length; i++) {
            System.out.println(f[i] + " " + z[i]);

        }

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    public static double[] getF() {
        return new double[]{1000000, 1291549.665, 1668100.537, 2154434.69, 2782559.402,
            3593813.664, 4641588.834, 5994842.503, 7742636.827, 1.00E+07, 1.02E+07, 1.05E+07,
            1.07E+07, 1.10E+07, 1.12E+07, 1.15E+07, 1.18E+07, 1.20E+07, 1.23E+07, 1.26E+07,
            1.29E+07, 1.32E+07, 1.35E+07, 1.38E+07, 1.42E+07, 1.45E+07, 1.48E+07, 1.52E+07,
            1.56E+07, 1.59E+07, 1.63E+07, 1.67E+07, 1.71E+07, 1.75E+07, 1.79E+07, 1.83E+07,
            1.87E+07, 1.92E+07, 1.96E+07, 2.01E+07, 2.06E+07, 2.10E+07, 2.15E+07, 2.21E+07,
            2.26E+07, 2.31E+07, 2.36E+07, 2.42E+07, 2.48E+07, 2.54E+07, 2.60E+07, 2.66E+07,
            2.72E+07, 2.78E+07, 2.85E+07, 2.92E+07, 2.98E+07, 3.05E+07, 3.13E+07, 3.20E+07,
            3.27E+07, 3.35E+07, 3.43E+07, 3.51E+07, 3.59E+07, 3.68E+07, 3.76E+07, 3.85E+07,
            3.94E+07, 4.04E+07, 4.13E+07, 4.23E+07, 4.33E+07, 4.43E+07, 4.53E+07, 4.64E+07,
            4.75E+07, 4.86E+07, 4.98E+07, 5.09E+07, 5.21E+07, 5.34E+07, 5.46E+07, 5.59E+07,
            5.72E+07, 5.86E+07, 5.99E+07, 6.14E+07, 6.28E+07, 6.43E+07, 6.58E+07, 6.73E+07,
            6.89E+07, 7.05E+07, 7.22E+07, 7.39E+07, 7.56E+07, 7.74E+07, 7.92E+07, 8.11E+07,
            8.30E+07, 8.50E+07, 8.70E+07, 8.90E+07, 9.11E+07, 9.33E+07, 9.55E+07, 9.77E+07, 1.00E+08,
            2e8, 3e8, 4e8, 5e8, 6e8, 7e8, 8e8, 9e8, 1e9
        };
    }
}
