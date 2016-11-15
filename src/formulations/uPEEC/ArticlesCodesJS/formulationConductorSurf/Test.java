/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.ArticlesCodesJS.formulationConductorSurf;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.LineRegion;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF.eps0;

/**
 *
 * @author jsiau
 */
public class Test {

    public static void main(String[] args) {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_SURF.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";

        double sigma = 1 / 1.68E-8;
        double epsilon_fr4 = 4.7;
        double epsilon_air = 1;
        double epsilon = epsilon_fr4 * eps0;
        double ep = 35e-6;
        //
        // Conducteur, Null, Dielectric, Borne Pos, Borne Neg.
        int[] ind = new int[]{1, 2, 0, -1, -3};
        int cct[] = null;

        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/NEWSERPENT_NB.DEC");
        SurfaceRegion conductor = (SurfaceRegion) mesh.getRegion(ind[0]);
        LineRegion bordNull = (LineRegion) mesh.getRegion(ind[1]);
        VolumeRegion Dielec = (VolumeRegion) mesh.getRegion(ind[2]);
        LineRegion borne1 = null, borne2 = null;
        if (ind[3] >= 0) {
            borne1 = (LineRegion) mesh.getRegion(ind[3]);
            borne2 = (LineRegion) mesh.getRegion(ind[4]);
        }
        //
        double f = 1e6;
        int formulation = 3;
        boolean negligeC = false;
        boolean negligeL = false;
        //
        PEEC_RLMPC_SURF solP = new PEEC_RLMPC_SURF(conductor, bordNull, sigma, ep, Dielec, epsilon, borne1, borne2,
                formulation, negligeC, negligeL);
        //
        solP.setPtsDeGaussInductifsVol(27, 27);// Hexaedre
//        solP.setPtsDeGaussInductifsSurf(64, 25, 25);// Quandrangle
//        solP.setPtsDeGaussInductifsSurf(25, 16, 16);// Quandrangle
             solP.setPtsDeGaussInductifsSurf(225, 64, 64);
        solP.setPtsDeGaussCapacitifs(64, 64);
        solP.setFullAnalyticalP(true);
        //
        //
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        int nbBranches = solP.getNbLignes();
        circuitPur.addSourceISimple(nbBranches, 1401, 1402, "Source I", 1.0, 0.0);
        cct = new int[]{1181, 1182};

        circuitPur.finSaisie();

        System.out.println("circuitPur=\n" + circuitPur.toStringOrdonance());

        solP.setCircuitElectrique(circuitPur);
        if (cct != null) {
            solP.setCct(cct[0], cct[1]);
        }

        FaceDeg1 FSd = solP.getFD1_Dielectric();
        FaceDeg1 FSc = solP.getFD1_Conductor();

        System.out.println("circuitPur.NbLigne= " + circuitPur.getNbLignes());
        System.out.println("solP.NbLigne= " + solP.getNbLignes());
        System.err.println("NbBranches= " + nbBranches);
        solP.setExportTopologie("D:");

        // Resolution
        solP.setCompression(false);// Enleve la compression
        double ib[][] = solP.resolutionIterative(f);
//        double ib[][] = solP.resolutionDirecte(f);
        System.err.println("Memory to store the matrices = " + solP.getMemoryUsed() / 1024 + " kO");
        //*/
        int nDof;
        if (formulation == 1 && negligeC) {
            nDof = FSc.getActiveDofCount();
        } else {
            nDof = FSd.getActiveDofCount() + FSc.getActiveDofCount();
        }
        System.out.println("nDof= " + nDof);
        Matrix res = new Matrix(2, nDof);
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }

        //* Scale seulement les DOF du conducteur
        res.submatrix(0, FSd.getActiveDofCount(), 2, FSc.getActiveDofCount()).scale(1 / ep);
        /*/ // Scale tous les dof
         res.scale(1/ep);
         //*/

        System.out.println("");
        System.out.println("");
        System.out.println("I= " + ib[0][2 * nbBranches] + " + j* " + ib[0][2 * nbBranches + 1]);
        System.out.println("U= " + ib[1][2 * nbBranches] + " + j* " + ib[1][2 * nbBranches + 1]);
        double Zmod = Math.hypot(ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]) / Math.hypot(ib[0][2 * nbBranches], ib[0][2 * nbBranches + 1]);
        System.out.println("|Z|= " + Zmod);
        double omega = 2 * Math.PI * f;
        System.out.println("C12= " + 1 / Zmod / omega);
        System.out.println("");
        System.out.println("");

        //
        //
        // EXPORT ALL THE QUANTITIES
        //
        //
        RealFaceQuantity Jreal_d = new RealFaceQuantity(FSd, res.row(0).subvector(0, FSd.getActiveDofCount()));
        RealFaceQuantity Jimag_d = new RealFaceQuantity(FSd, res.row(1).subvector(0, FSd.getActiveDofCount()));

        String path_JrealD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Snake/" + nDof + "f" + f + "_d_RE.msh";
        String path_JimD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Snake/" + nDof + "f" + f + "_d_IM.msh";
        String path_JmodD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Snake/" + nDof + "f" + f + "_d_MOD.msh";

        ExportGmshHdiv exportJreal = new ExportGmshHdiv(FSd, path_JrealD);
        exportJreal.addQuantity(Jreal_d, "Jreal_d");
        ExportGmshHdiv exportJimag = new ExportGmshHdiv(FSd, path_JimD);
        exportJimag.addQuantity(Jimag_d, "Jimag_d");
        ComplexFaceQuantity J = new ComplexFaceQuantity(FSd, res.submatrix(0, 0, 2, FSd.getActiveDofCount()));
        exportJreal = new ExportGmshHdiv(FSd, path_JmodD);
        exportJreal.addQuantityExportMod(J, "Jmod_d");

        //
        //
        //
        RealFaceQuantity Jreal_c = new RealFaceQuantity(FSc, res.row(0).subvector(FSd.getActiveDofCount(), FSc.getActiveDofCount()));
        RealFaceQuantity Jimag_c = new RealFaceQuantity(FSc, res.row(1).subvector(FSd.getActiveDofCount(), FSc.getActiveDofCount()));

        String path_JrealC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Snake/" + nDof + "f" + f + "_c_RE.msh";
        String path_JimC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Snake/" + nDof + "f" + f + "_c_IM.msh";
        String path_JmodC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Snake/" + nDof + "f" + f + "_c_MOD.msh";

        exportJreal = new ExportGmshHdiv(FSc, path_JrealC);
        exportJreal.addQuantity(Jreal_c, "Jreal_c");
        exportJimag = new ExportGmshHdiv(FSc, path_JimC);
        exportJimag.addQuantity(Jimag_c, "Jimag_c");
        J = new ComplexFaceQuantity(FSc, res.submatrix(0, FSd.getActiveDofCount(), 2, FSc.getActiveDofCount()));
        exportJreal = new ExportGmshHdiv(FSc, path_JmodC);
        exportJreal.addQuantityExportMod(J, "Jmod_c");

        String path_Merge = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Snake/" + nDof + "f" + f + "_MERGE.msh";
        try {
            Ecriture merge = new Ecriture(path_Merge);
            merge.ecrire("Merge '" + path_JrealC + "';\n");
            merge.ecrire("Merge '" + path_JimC + "';\n");
            merge.ecrire("Merge '" + path_JmodC + "';\n");

            merge.ecrire("Merge '" + path_JrealD + "';\n");
            merge.ecrire("Merge '" + path_JimD + "';\n");
            merge.ecrire("Merge '" + path_JmodD + "';\n");
            merge.close();

            File out_file = new File(path_Merge);
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file);
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_SURF.class.getName()).log(Level.SEVERE, null, ex);
        }
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
