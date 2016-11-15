/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.impedanceCurves;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.Tetraedre;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.LineRegion;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF.eps0;

/**
 *
 * @author jsiau
 */
public class Zf_PEEC_RLMPC_SURF_LOOPANTENA {

    public static void main(String[] args) {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Zf_PEEC_RLMPC_SURF_LOOPANTENA.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";

        // roCU = 1.68e-8 , roMetal = 2.836e-8
        double epsilon_fr4 = 4.7;
        double epsilon_air = 1;
        // Loop Antena
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/LA_WD_SURF_2CD.DEC");//SNAKE_DIELEC_SURF_4  LA_WD_SURF2
        double sigma = 1 / 1.72E-8;
        double epsilon = 5 * eps0;
        double ep = 2e-5;

        int ind[] = new int[]{1, 4, 0, 2, 3};

        PEEC_RLMPC_SURF solP = new PEEC_RLMPC_SURF((SurfaceRegion) mesh.getRegion(ind[0]), (LineRegion) mesh.getRegion(ind[1]), sigma, ep, (VolumeRegion) mesh.getRegion(ind[2]), epsilon, (LineRegion) mesh.getRegion(ind[3]), (LineRegion) mesh.getRegion(ind[4]),
                3, false, false);
        if (mesh.getRegion(0).getElementSet().getElements(0) instanceof Tetraedre) {
            solP.setPtsDeGaussInductifsVol(15, 4, 4);// Tetraedre
            solP.setPtsDeGaussInductifsSurf(13, 3, 3);// Triangle
            solP.setPtsDeGaussCapacitifs(3, 3);
        } else {
            solP.setPtsDeGaussInductifsVol(27, 8, 8);// Hexaedre
            solP.setPtsDeGaussInductifsSurf(64, 9, 9);// Quandrangle
            //        solP.setPtsDeGaussInductifsSurf(64, 144);// Quandrangle bugg
            solP.setPtsDeGaussCapacitifs(4, 4);
        }

//        Lecture lecf = null;
//        try {
//            lecf = new Lecture("d:/f.save");
//        } catch (IOException ex) {
//            Logger.getLogger(PEEC_RLMPC_SNAKE_Zf.class.getName()).log(Level.SEVERE, null, ex);
//        }
//        int nbF = lecf.getNbLignes();
//        System.out.println("Nombre de frequences a analyser= " + nbF);
//        double ff[] = new double[nbF];
//        for (int i = 0; i < nbF; i++) {
//            ff[i] = lecf.lireLigneDouble();
//        }
//        // 167
//        double f[] = new double[84];
//        for (int i = 0; i < f.length; i++) {
//            f[i] = ff[i + 2 * 167 + 83];
//        }
//        double fmin = f[0];
//        double fmax = f[f.length - 1];
        double fmin = 1e9;
        double fmax = 2.5e9;
        double f[] = compteFrequenciesLin(fmin, fmax, 30);
        try {
            Ecriture file = new Ecriture("D:/f_" + fmin + "_to_" + fmax + ".out");
            file.ecrire(Arrays.toString(f));
            file.close();
        } catch (IOException ex) {
            Logger.getLogger(Zf_PEEC_RLMPC_SURF_LOOPANTENA.class.getName()).log(Level.SEVERE, null, ex);
        }
//        
        System.out.println("f = [ " + fmin + " , " + fmax + " ];");
        System.out.println("f= " + Arrays.toString(f));
//        solP.setCompression(false);
//        solP.integration(f);
//        Matrix L[][] = solP.getL();
//        System.out.println("L00=\n" + L[0][0]);
//        Basic2D M = new Basic2D(solP.getZbFull(new double[solP.getNbLignes()][2 * solP.getNbLignes()], 0, solP.getNbLignes(), 0, solP.getNbLignes(), f));
//        System.out.println("M=\n" + M);
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        int indNodeInf = solP.getIndNodeInfinite();
        int B1[] = solP.getIndexNodeBorne1();
        int B2[] = solP.getIndexNodeBorne2();

        System.out.print("B1=[");
        for (int i = 0; i < B1.length; i++) {
            System.out.print(" " + B1[i]);
        }
        System.out.println("]");
        System.out.print("B2=[");
        for (int i = 0; i < B2.length; i++) {
            System.out.print(" " + B2[i]);
        }
        System.out.println("]");

        int indB1 = indNodeInf + B1.length + B2.length + 1;
        int indB2 = indB1 + 1;
        System.out.println("indCommunB1= " + indB1);

        System.out.println("indCommunB2= " + indB2);

        for (int i = 0; i < B1.length; i++) {
            circuitPur.addSourceUSimple(solP.getNbLignes() + i, B1[i], indB1, "Commun Borne1", 0, 0);
        }

        for (int i = 0; i < B2.length; i++) {
            circuitPur.addSourceUSimple(solP.getNbLignes() + B1.length + i, B2[i], indB2, "Commun Borne2", 0, 0);
        }

        int nbBranches = solP.getNbLignes() + B1.length + B2.length;
        // Le nombre de branche correspond aussi a lindice de lalimentation !
        circuitPur.addSourceISimple(nbBranches, indB1, indB2, "Source I", 1.0, 0.0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

        FaceDeg1 FSd = solP.getFD1_Dielectric();
        FaceDeg1 FSc = solP.getFD1_Conductor();

        System.out.println("circuitPur.NbLigne= " + circuitPur.getNbLignes());
        System.out.println("solP.NbLigne= " + solP.getNbLignes());
        System.err.println("NbBranches= " + nbBranches);

        solP.setCompression(false);// Enleve la compression
        double ib[][] = solP.resolutionDirecte(f);

        System.out.println(" Frequency = " + Arrays.toString(f) + ";");
        System.out.println("I=");
        for (int i = 0; i < ib.length / 2; i++) {
            System.out.println(ib[i][2 * nbBranches] + " + j* " + ib[i][2 * nbBranches + 1]);
        }
        System.out.println("U=");
        for (int i = 0; i < ib.length / 2; i++) {
            System.out.println(ib[i + (ib.length / 2)][2 * nbBranches] + " + j* " + ib[i + (ib.length / 2)][2 * nbBranches + 1]);
        }
        try {
            Ecriture file = new Ecriture("D:/Z_" + fmin + "_to_" + fmax + ".out");
            for (int i = 0; i < ib.length / 2; i++) {
                file.ecrire(f[i] + " " + Math.hypot(ib[i + (ib.length / 2)][2 * nbBranches], ib[i + (ib.length / 2)][2 * nbBranches + 1]) + " " + Math.atan(ib[i + (ib.length / 2)][2 * nbBranches + 1] / ib[i + (ib.length / 2)][2 * nbBranches]) + "\n");
            }
            file.close();
        } catch (IOException ex) {
            Logger.getLogger(Zf_PEEC_RLMPC_SURF_LOOPANTENA.class.getName()).log(Level.SEVERE, null, ex);
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

    protected static double[] compteFrequenciesLin(double fmin, double fmax, int nbF) {
        double f[] = new double[nbF];
        double step = (fmax - fmin) / (f.length - 1);
        for (int i = 0; i < f.length; i++) {
            f[i] = fmin + i * step;
        }
        return f;
    }

}
