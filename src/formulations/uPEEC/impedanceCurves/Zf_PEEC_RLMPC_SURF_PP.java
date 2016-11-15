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
public class Zf_PEEC_RLMPC_SURF_PP {

    public static void main(String[] args) {

        if (args == null) {
            GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);
        }

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Zf_PEEC_RLMPC_SURF_PP.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";

        // fr4
        double epsilon_fr4 = 4.7;

        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/2PP_SURF_NB.DEC");

        double sigma = 1 / 1.68E-8;
        double epsilon = epsilon_fr4 * eps0;
        double ep = 35e-6;

//        int ind[] = new int[]{1, 4, 0, 2, 3};
        int ind[] = new int[]{1, 2, 0};

        PEEC_RLMPC_SURF solP = new PEEC_RLMPC_SURF((SurfaceRegion) mesh.getRegion(ind[0]), (LineRegion) mesh.getRegion(ind[1]), sigma, ep, (VolumeRegion) mesh.getRegion(ind[2]), epsilon, null, null,
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
//                cct = new int[]{470, 471}; // Pour le maillage 2PP_SURF_NB
////                circuitPur.addSourceISimple(nbBranches, 889, 888, "Source I", 1.0, 0.0);
//                circuitPur.addSourceISimple(nbBranches, 889, 1000000000, "Source I", 1.0, 0.0);
//                circuitPur.addInductanceMutuelle(nbBranches+1, 1000000000, 888, "Source I", new double[]{1e-8}, new int[]{nbBranches+1});
//        
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
        double f[];
        double fmin = 1e6, fmax = 1e9;
        int nF = 40;
        if (args != null) {
            fmin = Double.valueOf(args[0]);
            fmax = Double.valueOf(args[1]);
            nF = Integer.valueOf(args[2]);
        }
        f = compteFrequenciesLin(fmin, fmax, nF);
        try {
            Ecriture file = new Ecriture("D:/f_" + fmin + "_to_" + fmax + ".out");
            file.ecrire(Arrays.toString(f));
            file.close();
        } catch (IOException ex) {
            Logger.getLogger(Zf_PEEC_RLMPC_SURF_PP.class.getName()).log(Level.SEVERE, null, ex);
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
        solP.setCct(470, 471);

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int nbBranches = solP.getNbLignes();
////                circuitPur.addSourceISimple(nbBranches, 889, 888, "Source I", 1.0, 0.0);
//        circuitPur.addSourceISimple(nbBranches, 889, 1000000000, "Source I", 1.0, 0.0);
////        circuitPur.addInductanceMutuelle(nbBranches + 1, 1000000000, 888, "Source I", new double[]{1e-8}, new int[]{nbBranches + 1});        
////        circuitPur.addCapacite(nbBranches+2, 888, 889, "C12", 262e-12);
//        circuitPur.addResistance(nbBranches + 1, 1000000000, 1000000000 + 1, "R", 1e-3);
//        circuitPur.addInductanceMutuelle(nbBranches + 2, 1000000000 + 1, 888, "Source I", new double[]{10e-9}, new int[]{nbBranches + 2});
//        circuitPur.addCapacite(nbBranches + 3, 888, 889, "C12", 262e-12);
        //*
        circuitPur.addSourceISimple(nbBranches, 889, 888, "Source I", 1.0, 0.0);
        /*/
         int ndInter = 1000000000;
         circuitPur.addSourceISimple(nbBranches, 889, ndInter, "Source I", 1.0, 0.0);
         //        circuitPur.addResistance(nbBranches + 1, ndInter++, ndInter, "R", 100e-3);
         //        circuitPur.addInductanceMutuelle(nbBranches + 2, ndInter, 888, "Source I", new double[]{100e-9}, new int[]{nbBranches + 2});
         circuitPur.addCapacite(nbBranches + 3, 888, 889, "C12", 262e-12);//;262e-12);
         //*/

//        int indNodeInf = solP.getIndNodeInfinite();
//        int B1[] = solP.getIndexNodeBorne1();
//        int B2[] = solP.getIndexNodeBorne2();
//
//        System.out.print("B1=[");
//        for (int i = 0; i < B1.length; i++) {
//            System.out.print(" " + B1[i]);
//        }
//        System.out.println("]");
//        System.out.print("B2=[");
//        for (int i = 0; i < B2.length; i++) {
//            System.out.print(" " + B2[i]);
//        }
//        System.out.println("]");
//
//        int indB1 = indNodeInf + B1.length + B2.length + 1;
//        int indB2 = indB1 + 1;
//        System.out.println("indCommunB1= " + indB1);
//
//        System.out.println("indCommunB2= " + indB2);
//
//        for (int i = 0; i < B1.length; i++) {
//            circuitPur.addSourceUSimple(solP.getNbLignes() + i, B1[i], indB1, "Commun Borne1", 0, 0);
//        }
//
//        for (int i = 0; i < B2.length; i++) {
//            circuitPur.addSourceUSimple(solP.getNbLignes() + B1.length + i, B2[i], indB2, "Commun Borne2", 0, 0);
//        }
//
//        int nbBranches = solP.getNbLignes() + B1.length + B2.length;
//        // Le nombre de branche correspond aussi a lindice de lalimentation !
//        circuitPur.addSourceISimple(nbBranches, indB1, indB2, "Source I", 1.0, 0.0);
        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);
        System.out.println("Circuit Pur= ");
        System.out.println(circuitPur.toString());

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
        System.out.println("f \t |Z| \t Phase");
        for (int i = 0; i < ib.length / 2; i++) {
            System.out.println(f[i] + " " + Math.hypot(ib[i + (ib.length / 2)][2 * nbBranches], ib[i + (ib.length / 2)][2 * nbBranches + 1]) + " " + Math.atan(ib[i + (ib.length / 2)][2 * nbBranches + 1] / ib[i + (ib.length / 2)][2 * nbBranches]));
        }
        try {
            Ecriture file = new Ecriture("D:/Z_" + fmin + "_to_" + fmax + ".out");
            for (int i = 0; i < ib.length / 2; i++) {
                file.ecrire(f[i] + " " + Math.hypot(ib[i + (ib.length / 2)][2 * nbBranches], ib[i + (ib.length / 2)][2 * nbBranches + 1]) + " " + Math.atan(ib[i + (ib.length / 2)][2 * nbBranches + 1] / ib[i + (ib.length / 2)][2 * nbBranches]) + "\n");
            }
            file.close();
        } catch (IOException ex) {
            Logger.getLogger(Zf_PEEC_RLMPC_SURF_PP.class.getName()).log(Level.SEVERE, null, ex);
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
