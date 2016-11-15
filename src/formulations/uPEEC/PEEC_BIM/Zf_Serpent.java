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
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;

/**
 *
 * @author jsiau
 */
public class Zf_Serpent {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Zf_Serpent.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        int casem;
        /*
         ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Surf/S4_FS.DEC");
         casem = 0;
         /*/
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Surf/S6_NB_FS.DEC");
        casem = 1;
        //*/
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
        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((SurfaceRegion) mesh.getRegion(0), 1 / 1.68e-8, 35e-6, null, null,// (LineRegion) mesh.getRegion(2), (LineRegion) mesh.getRegion(3),
                (SurfaceRegion) mesh.getRegion(1), 4.7 * eps0,
                (SurfaceRegion) mesh.getRegion(0)
        );

//        solP.setPtsDeGaussInductifs(25, 16, 16);
//        solP.setPtsDeGaussCapacitifs(25, 64);
        solP.setPtsDeGaussInductifs(225, 64, 64);
////        solP.setPtsDeGaussCapacitifs(4, 4);
////        solP.setAnalyticalIntegrationCapa(true);
//        solP.setPtsDeGaussInductifs(9, 4, 4);
        solP.setPtsDeGaussCapacitifs(16, 16);
        solP.setAnalyticalIntegrationCapa(true);
//        solP.setPtsDeGaussInductifs(64, 25, 25);
//        solP.setPtsDeGaussCapacitifs(25, 64);
//        solP.setAnalyticalIntegrationCapa(true);
//        solP.setCorrectionVoisin(true);

//        solP.setEquiPot(new int[]{512, 1212});
        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int indNoeudBord = fd.getNbElement();
        int nbBranches = solP.getNbLignes();
//        circuitPur.addSourceUSimple(nbBranches, indNoeudBord, indNoeudBord + 1, "Source", 240.0, 0.0);
        if (casem == 0) {
            circuitPur.addSourceUSimple(nbBranches, 699, 1399, "Source", 1.0, 0.0);
        } else {
            circuitPur.addSourceUSimple(nbBranches, 1625, 3251, "Source", 1.0, 0.0);
        }
//        circuitPur.addSourceUSimple(nbBranches+1, 512, 1212, "Source", 0.0, 0.0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

        double f[];
        double fmin, fmax;
        int nF;
        if (args.length != 0) {
            double emin = Double.valueOf(args[0]);
            double emax = Double.valueOf(args[1]);
            nF = Integer.valueOf(args[2]);
            f = Zf_Serpent.compteFrequenciesLog(emin, emax, nF);
            fmin = Math.pow(10, emin);
            fmax = Math.pow(10, emax);

        } else {
            fmin = 1e6;
            fmax = 1e9;
            nF = 100;
            f = compteFrequenciesLog(6, 9, nF);
        }
        System.out.println("f.length= " + f.length);
        double Zmod[] = new double[f.length];

        boolean bull = true;
        for (int i = 0; i < f.length; i++) {
            System.out.println("Frequency: " + (i + 1) + " / " + f.length);
            double[][] ib = solP.resolutionDirecte(f[i]);

            double I[] = new double[]{ib[0][2 * nbBranches], ib[0][2 * nbBranches + 1]};
            double U[] = new double[]{ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]};
            System.out.println("I= " + Arrays.toString(I));
            System.out.println("U= " + Arrays.toString(U));
            Zmod[i] = Math.hypot(U[0], U[1]) / Math.hypot(I[0], I[1]);
            System.out.println("|Z|= " + Zmod[i]);

            Matrix res = new Matrix(2, solP.getNbLignes());
            for (int j = 0; j < res.getColumnCount(); j++) {
                res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(j), ib[0][2 * j] / 35e-6);
                res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(j), ib[0][2 * j + 1] / 35e-6);
            }
            if (bull) {
                PEEC_FULLY_SURF.exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/SERPENT/f_" + f[i], fd, res, solP.getDielecCell(), solP.getQ(), bull);
                bull = false;
            }
        }

        for (int i = 0; i < f.length; i++) {
            System.out.println(f[i] + " " + Zmod[i]);
        }

        try {
            Ecriture file = new Ecriture("D:/Z_" + fmin + "_to_" + fmax + ".out");
            for (int i = 0; i < f.length; i++) {
                file.ecrire(f[i] + " " + Zmod[i] + "\n");
            }
            file.close();
        } catch (IOException ex) {
            Logger.getLogger(Zf_Serpent.class.getName()).log(Level.SEVERE, null, ex);
        }

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    protected static double[] compteFrequenciesLog(double logfmin, double logfmax, int nbF) {
        double f[] = new double[nbF];
        double steplog = (logfmax - logfmin) / (f.length - 1);
        for (int i = 0; i < f.length; i++) {
            f[i] = Math.pow(10, logfmin + i * steplog);
        }
        return f;
    }

}
