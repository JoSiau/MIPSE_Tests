/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.impedanceCurves;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static formulations.uPEEC.impedanceCurves.Zf_PEEC_RLMPC_SURF_SNAKE_CID.compteFrequenciesLin;
import static formulations.uPEEC.impedanceCurves.Zf_PEEC_RLMPC_SURF_SNAKE_CID.compteFrequenciesLog;
import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume.eps0;

/**
 *
 * @author jsiau
 */
public class Zf_RLMPCvol_CondInDielec {

    public static void main(String[] args) {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_Volume.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Volumic/SNAKE_IDB.DEC");//T1B.DEC");
        mesh.getRegions();
        VolumeRegion Cond = (VolumeRegion) mesh.getRegion(1);
        VolumeRegion Dielec = (VolumeRegion) mesh.getRegion(0);

        SurfaceRegion b1 = (SurfaceRegion) mesh.getRegion(2);
        SurfaceRegion b2 = (SurfaceRegion) mesh.getRegion(3);

        double sigma = 1 / 1.68e-8;
        double epsilon_fr4 = 4.7;
        double epsilon_air = 1;
        double epsilon = 4.7 * eps0;
        double ep = 35e-6;
        /*
         *******************************
         Nombre de regions importees : 5
         Region 0, Nom : COND, type : 3, 49 elements
         Region 1, Nom : DIELEC, type : 3, 1134 elements
         Region 2, Nom : B1, type : 2, 7 elements
         Region 3, Nom : B2, type : 2, 7 elements
         Region 4, Nom : NULL, type : 2, 112 elements
         *******************************
         */

        PEEC_RLMPC_Volume solP = new PEEC_RLMPC_Volume(Cond, sigma,
                new VolumeRegion[]{Dielec}, new double[]{epsilon},
                b1, b2, null, null);

        solP.setPtsDeGaussInductifsVol(512, 125, 125);
        solP.setPtsDeGaussCapacitifs(25, 25);
        solP.setFullAnalyticalP(true);

//        solP.setCheckNaN(true);
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        int indNodeInf = solP.getIndNodeInfinite();
        int B1[] = solP.getIndexNodeBorne1();
        int B2[] = solP.getIndexNodeBorne2();
        int nbBranches;

        int indB1 = indNodeInf + B1.length + B2.length + 1;
        int indB2 = indB1 + 1;
        for (int i = 0; i < B1.length; i++) {
            circuitPur.addSourceUSimple(solP.getNbLignes() + i, B1[i], indB1, "Commun Borne1", 0, 0);
        }
        for (int i = 0; i < B2.length; i++) {
            circuitPur.addSourceUSimple(solP.getNbLignes() + B1.length + i, B2[i], indB2, "Commun Borne2", 0, 0);
        }
        nbBranches = solP.getNbLignes() + B1.length + B2.length;
        // Le nombre de branche correspond aussi a lindice de lalimentation !
        circuitPur.addSourceISimple(nbBranches, indB1, indB2, "Source I", 1.0, 0.0);

        circuitPur.finSaisie();

        System.out.println("circuitPur=\n" + circuitPur.toStringOrdonance());

        solP.setCircuitElectrique(circuitPur);

        FaceDeg1 FSd = solP.getFD1_Dielectric();
        FaceDeg1 FSc = solP.getFD1_Conductor();

        double fmin = 1e6;
        double fmax = 1e9;
        int nF = 100;
        double f[];
        if (args.length != 0) {
            fmin = Double.valueOf(args[0]);
            fmax = Double.valueOf(args[1]);
            nF = Integer.valueOf(args[2]);
            f = compteFrequenciesLin(fmin, fmax, nF);
        } else {
            f = compteFrequenciesLog(6, 9, 100);
        }
        try {
            Ecriture file = new Ecriture("D:/f_" + fmin + "_to_" + fmax + ".out");
            file.ecrire(Arrays.toString(f));
            file.close();
        } catch (IOException ex) {
            Logger.getLogger(Zf_RLMPCvol_CondInDielec.class.getName()).log(Level.SEVERE, null, ex);
        }
//        
        System.out.println("f = [ " + fmin + " , " + fmax + " ];");
        System.out.println("f= " + Arrays.toString(f));
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // Resolution
        /*
         boolean HmatComp = false;
         solP.setCompression(HmatComp);// Enleve la com pression
         solP.setParamIterativeSolver(1, 100, -1e-10, 1, -1);
         solP.setParamPreconditionner(2, 500, 0, new double[]{1, 1e-4, -1});
         double ib[][] = solP.resolutionIterative(f);
         /*/
        boolean HmatComp = false;
        solP.setCompression(HmatComp);// Enleve la compression
        double ib[][] = solP.resolutionDirecte(f);
        //*/
        int nDof = FSd.getActiveDofCount() + FSc.getActiveDofCount();

        System.out.println("nDof= " + nDof);
        Matrix res = new Matrix(2, nDof);
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }

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
            Logger.getLogger(Zf_RLMPCvol_CondInDielec.class.getName()).log(Level.SEVERE, null, ex);
        }

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
