/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.surfaceImpedance;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.SurfaceImpedance.PEEC_RLMPC_SURFACE_IMPEDANCE;
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
import java.util.Arrays;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;

/**
 *
 * @author jsiau
 */
public class TwoDielec {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(TwoDielec.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/2D.DEC");

        SurfaceRegion conductor = (SurfaceRegion) mesh.getRegion(2);
        LineRegion bordNull = (LineRegion) mesh.getRegion(3);
        VolumeRegion Dielec[] = new VolumeRegion[]{(VolumeRegion) mesh.getRegion(0), (VolumeRegion) mesh.getRegion(1)};
        LineRegion borne1 = (LineRegion) mesh.getRegion(4);
        LineRegion borne2 = (LineRegion) mesh.getRegion(5);

        double epsilon[] = new double[]{4.7 * eps0, 100 * eps0};

        PEEC_RLMPC_SURFACE_IMPEDANCE solP = new PEEC_RLMPC_SURFACE_IMPEDANCE(conductor, bordNull, 1 / 1.68e-8,
                Dielec, epsilon,
                borne1, borne2,
                3, false, true);
        
        
//        solP.setPtsDeGaussInductifsVol(125, 64, 64);
        solP.setPtsDeGaussInductifsVol(27, 8, 8);
        solP.setPtsDeGaussInductifsSurf(25, 16, 16);
//        solP.setPtsDeGaussInductifsSurf(64, 25, 25);
        solP.setPtsDeGaussCapacitifs(16, 16);
        solP.setFullAnalyticalP(true);

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

        FaceDeg1 FSd = solP.getFD1_Dielectric();
        FaceDeg1 FSc = solP.getFD1_Conductor();

//        solP.setExportTopologie("D:");
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        double f = 1e6;
        Scanner sc = new Scanner(System.in);
        boolean HmatComp;
        do {
            // Resolution
        /*
             boolean HmatComp = true;
             solP.setCompression(HmatComp);
             solP.setParamIterativeSolver(1, 10000, 1e-3, 1, -150);
             solP.setParamPreconditionner(2, 1500, 0, new double[]{1, 1e-4, 30000});
             double ib[][] = solP.resolutionIterative(f);
             /*/
            HmatComp = false;
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

            try {
                // Save the datas.
                Ecriture saveRes = new Ecriture("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/Surface_Impedance/f" + f + "_RES.out");
                for (int i = 0; i < 2; i++) {
                    double[] tmp = new double[nDof];
                    res.row(i).get(tmp);
                    saveRes.ecrire(tmp, ',');
                    saveRes.aLaLigne();
                }
                saveRes.close();
            } catch (IOException ex) {
                Logger.getLogger(TwoDielec.class.getName()).log(Level.SEVERE, null, ex);
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

            String path_JrealD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/Surface_Impedance/" + nDof + "f" + f + "_d_RE.msh";
            String path_JimD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/Surface_Impedance/" + nDof + "f" + f + "_d_IM.msh";
            String path_JmodD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/Surface_Impedance/" + nDof + "f" + f + "_d_MOD.msh";

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

            String path_JrealC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/Surface_Impedance/" + nDof + "f" + f + "_c_RE.msh";
            String path_JimC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/Surface_Impedance/" + nDof + "f" + f + "_c_IM.msh";
            String path_JmodC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/Surface_Impedance/" + nDof + "f" + f + "_c_MOD.msh";

            exportJreal = new ExportGmshHdiv(FSc, path_JrealC);
            exportJreal.addQuantity(Jreal_c, "Jreal_c");
            exportJimag = new ExportGmshHdiv(FSc, path_JimC);
            exportJimag.addQuantity(Jimag_c, "Jimag_c");
            J = new ComplexFaceQuantity(FSc, res.submatrix(0, FSd.getActiveDofCount(), 2, FSc.getActiveDofCount()));
            exportJreal = new ExportGmshHdiv(FSc, path_JmodC);
            exportJreal.addQuantityExportMod(J, "Jmod_c");

            String path_Merge = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/Surface_Impedance/" + nDof + "f" + f + "_MERGE.msh";
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
                Logger.getLogger(PEEC_RLMPC_Volume.class.getName()).log(Level.SEVERE, null, ex);
            }
            System.out.println("Entrez -1 pour quitter ou une autre frequence a resoudre: ");
            f = sc.nextDouble();
        } while (f >= 0);

        

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}