/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.impedanceCurves;

import g2elab.mipse.circuit.solverCircuitComplex.SolveurMaillesIndependantes;
import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.SurfaceImpedance.PEEC_RLMPC_SURFACE_IMPEDANCE;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMesh;
import g2elab.mipse.meshCore.region.LineRegion;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;
import static formulations.uPEEC.impedanceCurves.Zf_PEEC_RLMPC_SURF_SERPENT.compteFrequenciesLog;

/**
 *
 * @author jsiau
 */
public class SaveMatrixSolutions {

    public static void main(String args[]) {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(4);

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_Volume.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportGmshMesh mesh = new ImportGmshMesh(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/bob16.msh");
        mesh.meshSummary();
        //
        // Gather the meshes
        //
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
        LineRegion fluxPos = (LineRegion) ((SurfaceRegion) mesh.getRegion(4)).generateBorder();
        LineRegion fluxNeg = (LineRegion) ((SurfaceRegion) mesh.getRegion(5)).generateBorder();

        PEEC_RLMPC_SURFACE_IMPEDANCE solP = new PEEC_RLMPC_SURFACE_IMPEDANCE(interfaceCondDielec, null, 5.8502e7,
                dielectrics, new double[]{3.2 * eps0, 3.0 * eps0},
                fluxPos, fluxNeg, 3, false, true);
        /*
         solP.setPtsDeGaussInductifsVol(27, 8, 8);
         solP.setPtsDeGaussInductifsSurf(25, 16, 16);
         solP.setPtsDeGaussCapacitifs(4, 4);
         /*/
        solP.setPtsDeGaussInductifsVol(225, 100, 100);
        solP.setPtsDeGaussInductifsSurf(45, 36, 36);
        solP.setPtsDeGaussCapacitifs(45, 45);
        solP.setFullAnalyticalP(true);
        //*/
        //
        //
        //
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int indNodeInf = solP.getIndNodeInfinite();
        int B1[] = solP.getIndexNodeBorne1();
        System.out.println("B1= " + Arrays.toString(B1));
        int B2[] = solP.getIndexNodeBorne2();
        System.out.println("B2= " + Arrays.toString(B2));
        int nbBranches;
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
        nbBranches = solP.getNbLignes() + B1.length + B2.length;
        // Le nombre de branche correspond aussi a lindice de lalimentation !
        circuitPur.addSourceISimple(nbBranches, indB1, indB2, "Source I", 1.0, 0.0);
        circuitPur.finSaisie();

        solP.setCircuitElectrique(circuitPur);
        solP.setCompression(false);
//        solP.setExportTopologie("D:");
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        double fmin = 40;
        double fmax = 1.1e8;
        int nF = 45;
        double f[];
        f = new double[nF];
        f[0] = 40;
        f[1] = 1e3;
        f[2] = 1e4;
        f[3] = 1e5;
        System.arraycopy(compteFrequenciesLog(6, 8, 40), 0, f, 4, 40);
        f[44] = 1.1e8;

//        double ff[] = new double[24];
//        System.arraycopy(f, 10, ff, 0, ff.length);
//        f = ff;
        nF = f.length;
        fmin = f[0];
        fmax = f[nF - 1];
        System.out.println("f = " + Arrays.toString(f));
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        solP.integration(f[0]);

        SolveurMaillesIndependantes solver = solP.getSolveurCircuit();
        for (int i = 0; i < nF; i++) {

            double b[] = solver.constructionSecondMembre(f[i]);
            saveVectorComplex(b, "D:/tmp/b_f" + f[i]);

            double z[][] = solver.getMatriceAResoudreFull(f[i]);
            saveMatrixComplex(z, "D:/tmp/z_f" + f[i]);

        }
    }


    protected static void saveMatrixComplex(double z[][], String file) {
        try {
            Ecriture ec = new Ecriture(file);

            for (int i = 0; i < z.length; i++) {
                ec.ecrire(z[i], ' ');
                ec.aLaLigne();
            }
            ec.close();
        } catch (IOException ex) {
            Logger.getLogger(SaveMatrixSolutions.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    protected static void saveVectorComplex(double b[], String file) {
        try {
            Ecriture ec = new Ecriture(file);

            ec.ecrire(b, ' ');

            ec.close();
        } catch (IOException ex) {
            Logger.getLogger(SaveMatrixSolutions.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

}
