/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.impedanceCurves;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMesh;
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

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;
import static formulations.uPEEC.impedanceCurves.Zf_PEEC_RLMPC_SURF_SERPENT.compteFrequenciesLog;
import g2elab.mipse.mipseCore.matrixCompression.Compression;

/**
 *
 * @author jsiau
 */
public class Zf_BobZacharie {

    public static void main(String[] args) {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Zf_BobZacharie.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportGmshMesh mesh = new ImportGmshMesh(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/bob2.msh");
        mesh.meshSummary();
        //
        // Gather the meshes
        //
        // Conductor
        VolumeRegion conductor = (VolumeRegion) mesh.getRegion(6);
        // Put the  2 dieledctric regions in an array
        VolumeRegion dielectrics[] = new VolumeRegion[]{(VolumeRegion) mesh.getRegion(7), (VolumeRegion) mesh.getRegion(8)};
        // Interface between the conductor and the dielectrics
        SurfaceRegion interfaceCondDielec = (SurfaceRegion) mesh.getRegion(1);
        // Get the dielectrics border
        SurfaceRegion bobAir = (SurfaceRegion) mesh.getRegion(2);
        SurfaceRegion formerAir = (SurfaceRegion) mesh.getRegion(3);
        // And compute their union
        SurfaceRegion dielecAir = (SurfaceRegion) bobAir.generateUnionWith(formerAir);
        dielecAir.exportGmsh("d:/dielecAir.msh");
        SurfaceRegion fluxPos = (SurfaceRegion) mesh.getRegion(4);
        SurfaceRegion fluxNeg = (SurfaceRegion) mesh.getRegion(5);

        boolean computeCapa = true;

        PEEC_RLMPC_Volume solP = new PEEC_RLMPC_Volume(conductor, 5.8502e7,
                dielectrics, new double[]{2.8 * eps0, 3.0 * eps0},
                fluxPos, fluxNeg, interfaceCondDielec, dielecAir,
                computeCapa);

        solP.setPtsDeGaussInductifsVol(64, 27, 27);
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

        int nF = 5;
        double f[];
        /*
         if (args.length != 0) {
         fmin = Double.valueOf(args[0]);
         fmax = Double.valueOf(args[1]);
         nF = Integer.valueOf(args[2]);
         f = compteFrequenciesLin(fmin, fmax, nF);
         } else {
         f = compteFrequenciesLog(6, 7, nF);
         }
         /*/
//        nF = 12;
//        f = new double[nF];
//        System.arraycopy(compteFrequenciesLog(6, 8, 12), 0, f, 0, 12);
////////////////////////////////////////////////////////////////////////////////
//        nF = 13;
//        f = new double[nF];
//        f[0] = 40;
//        System.arraycopy(compteFrequenciesLog(2, 7, 12), 0, f, 1, 12);
////////////////////////////////////////////////////////////////////////////////
        nF = 45;
        f = new double[nF];
        f[0] = 40;
        f[1] = 1e3;
        f[2] = 1e4;
        f[3] = 1e5;
        System.arraycopy(compteFrequenciesLog(6, 8, 40), 0, f, 4, 40);
        f[44] = 1.1e8;

        double ff[] = new double[10];
        System.arraycopy(f, 0, ff, 0, ff.length);
        f = ff;
        //*/
        nF = f.length;
        double fmin = f[0];
        double fmax = f[nF - 1];
        try {
            Ecriture file = new Ecriture("D:/f_" + fmin + "_to_" + fmax + ".out");
            file.ecrire(Arrays.toString(f));
            file.close();
        } catch (IOException ex) {
            Logger.getLogger(Zf_BobZacharie.class.getName()).log(Level.SEVERE, null, ex);
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
        boolean HmatComp = true;
        solP.setCompression(1e-4, 100, 50, 3);// Active la compression
        solP.setParamIterativeSolver(1, 100, -1e-8, 1, -1);
        solP.setParamPreconditionner(2, 1000, 0, new double[]{1, 1e-4, -1});
        double ib[][] = solP.resolutionIterative(f);
        /*/
         boolean HmatComp = false;
         solP.setCompression(Compression.No);// Enleve la compression
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
            Logger.getLogger(Zf_BobZacharie.class.getName()).log(Level.SEVERE, null, ex);
        }

        GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
