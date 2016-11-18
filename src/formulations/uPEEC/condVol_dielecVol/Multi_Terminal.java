/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.condVol_dielecVol;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume.eps0;
import g2elab.mipse.mipseCore.matrixCompression.Compression;

/**
 *
 * @author jsiau
 */
public class Multi_Terminal {

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        System.out.println("\nChoisir la frequence: ");
        double f = sc.nextDouble();

        GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Multi_Terminal.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Volumic/MULTI_TERM.DEC");//T1B.DEC");
        mesh.getRegions();
        VolumeRegion Cond = (VolumeRegion) mesh.getRegion(1);
        VolumeRegion Dielec = (VolumeRegion) mesh.getRegion(0);

        SurfaceRegion b1p = (SurfaceRegion) mesh.getRegion(2);
        SurfaceRegion b1m = (SurfaceRegion) mesh.getRegion(4);
        SurfaceRegion b2p = (SurfaceRegion) mesh.getRegion(3);
        SurfaceRegion b2m = (SurfaceRegion) mesh.getRegion(5);

        SurfaceRegion border = (SurfaceRegion) mesh.getRegion(6);
        SurfaceRegion inter = (SurfaceRegion) mesh.getRegion(7);

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
                new SurfaceRegion[]{b1p, b1m, b2p, b2m}, inter, border, true);

        solP.setPtsDeGaussInductifsVol(27, 8, 8);
        solP.setPtsDeGaussCapacitifs(4, 4);
//            solP.setFullAnalyticalP(true);

        solP.setCheckNaN(true);
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        int indNodeInf = solP.getIndNodeInfinite();
        int B[] = solP.getIndexTerminals();
        int nbBranches;
        
        circuitPur.addSourceISimple(solP.getNbLignes(), B[0], B[1], "Source I1", 1.0, 0.0);
        circuitPur.addSourceISimple(solP.getNbLignes(), B[2], B[3], "Source I2", 1.0, 0.0);
        nbBranches = solP.getNbLignes();
        
        circuitPur.finSaisie();

        System.out.println("circuitPur=\n" + circuitPur.toStringOrdonance());

        solP.setCircuitElectrique(circuitPur);

        FaceDeg1 FSd = solP.getFD1_Dielectric();
        FaceDeg1 FSc = solP.getFD1_Conductor();

        System.out.println("circuitPur.NbLigne= " + circuitPur.getNbLignes());
        System.out.println("solP.NbLigne= " + solP.getNbLignes());
        System.err.println("NbBranches= " + nbBranches);
        solP.setExportTopologie("D:");
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // Resolution
        /*
         boolean HmatComp = false;
         solP.setCompression(HmatComp?Compression.HCA:Compression.No);// Enleve la com pression
         solP.setParamIterativeSolver(1, 100, -1e-10, 1, -1);
         solP.setParamPreconditionner(2, 500, 0, new double[]{1, 1e-4, -1});
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

        try {
            // Save the datas.
            Ecriture saveRes = new Ecriture("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/CondInDielec/f" + f + "_RES.out");
            for (int i = 0; i < 2; i++) {
                double[] tmp = new double[nDof];
                res.row(i).get(tmp);
                saveRes.ecrire(tmp, ',');
                saveRes.aLaLigne();
            }
            saveRes.close();
        } catch (IOException ex) {
            Logger.getLogger(Multi_Terminal.class.getName()).log(Level.SEVERE, null, ex);
        }

        System.out.println("");
        System.out.println("");
        System.out.println("I= " + ib[0][2 * nbBranches] + " + j* " + ib[0][2 * nbBranches + 1]);
        System.out.println("U= " + ib[1][2 * nbBranches] + " + j* " + ib[1][2 * nbBranches + 1]);
        double Zmod = Math.hypot(ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]);
        System.out.println("|Z|= " + Zmod);
        double omega = 2 * Math.PI * f;
        System.out.println("C12= " + 1 / Zmod / omega);
        System.out.println("");
        System.out.println("");

//        System.out.println("P:\n"+solP.getP());
        //
        //
        // EXPORT ALL THE QUANTITIES
        //
        //
        RealFaceQuantity Jreal_d = new RealFaceQuantity(FSd, res.row(0).subvector(0, FSd.getActiveDofCount()));
        RealFaceQuantity Jimag_d = new RealFaceQuantity(FSd, res.row(1).subvector(0, FSd.getActiveDofCount()));

        String path_JrealD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/CondInDielec/" + nDof + "f" + f + "_d_RE.msh";
        String path_JimD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/CondInDielec/" + nDof + "f" + f + "_d_IM.msh";
        String path_JmodD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/CondInDielec/" + nDof + "f" + f + "_d_MOD.msh";

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

        String path_JrealC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/CondInDielec/" + nDof + "f" + f + "_c_RE.msh";
        String path_JimC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/CondInDielec/" + nDof + "f" + f + "_c_IM.msh";
        String path_JmodC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/CondInDielec/" + nDof + "f" + f + "_c_MOD.msh";

        exportJreal = new ExportGmshHdiv(FSc, path_JrealC);
        exportJreal.addQuantity(Jreal_c, "Jreal_c");
        exportJimag = new ExportGmshHdiv(FSc, path_JimC);
        exportJimag.addQuantity(Jimag_c, "Jimag_c");
        J = new ComplexFaceQuantity(FSc, res.submatrix(0, FSd.getActiveDofCount(), 2, FSc.getActiveDofCount()));
        exportJreal = new ExportGmshHdiv(FSc, path_JmodC);
        exportJreal.addQuantityExportMod(J, "Jmod_c");

        String path_Merge = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/CondInDielec/" + nDof + "f" + f + "_MERGE.msh";
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
            Logger.getLogger(Multi_Terminal.class.getName()).log(Level.SEVERE, null, ex);
        }

        double b[] = new double[2 * solP.getNbLignes()];
        for (int i = 0; i < b.length; i++) {
            b[i] = Math.random() * 100;
        }
        double pf[] = solP.produit(b, new double[2 * solP.getNbLignes()], f);

        if (HmatComp) {
            solP.setCompression(Compression.No);
            solP.integrationMono();
        }
        Basic2D Mf = new Basic2D(solP.getZbFull(new double[solP.getNbLignes()][2 * solP.getNbLignes()], 0, solP.getNbLignes() - 1, 0, solP.getNbLignes() - 1, f));

        double rf[] = Mf.product(b, new double[2 * solP.getNbLignes()]);

        ColumnVector vf = new ColumnVector(rf);
        ColumnVector vp = new ColumnVector(pf);

        vp.sub(vf);
        System.out.println("ErrorRel pmv's= " + vp.norm() / vf.norm());
        System.out.println("ErrorAbs pmv's= " + vp.norm());

        GestionnaireTaches.getGestionnaireTaches().stop();

    }

}
