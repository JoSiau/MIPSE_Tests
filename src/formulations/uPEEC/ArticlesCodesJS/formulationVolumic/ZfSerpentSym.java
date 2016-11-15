/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.ArticlesCodesJS.formulationVolumic;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.elements.HexaedreDroit;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.File;
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
public class ZfSerpentSym {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        HexaedreDroit.setAdaptatif(true);
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(7);

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(ZfSerpentSym.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/ArticlesCodesJS/FormulationVolumique/SERP_SYM_VOL.DEC";
        ImportFlux mesh = new ImportFlux(meshDir);
        //
        // Gather the meshes
        //
        // Conductor
        VolumeRegion conductor = null;
        // Put the dieledctric regions in an array
        VolumeRegion dielectrics[] = null;
        if (meshDir.endsWith("SERP_SYM_VOL.DEC")) {
            // Conductor
            conductor = (VolumeRegion) mesh.getRegion(1);
            // Put the dieledctric regions in an array
            dielectrics = new VolumeRegion[]{(VolumeRegion) mesh.getRegion(0)};
        } else if (meshDir.endsWith("SERP_SYM_VOL_5.DEC")) {
            // Conductor
            conductor = (VolumeRegion) mesh.getRegion(0);
            // Put the dieledctric regions in an array
            dielectrics = new VolumeRegion[]{(VolumeRegion) mesh.getRegion(1)};
        }
        //
//        SurfaceRegion fluxPos = (SurfaceRegion) mesh.getRegion(2);
//        SurfaceRegion fluxNeg = (SurfaceRegion) mesh.getRegion(3);

        boolean computeCapa = true;

        double sigma = 1 / 1.72e-8;
        PEEC_RLMPC_Volume solP = new PEEC_RLMPC_Volume(conductor, sigma,
                dielectrics, new double[]{4.7 * eps0},
                null, null, null, null,
                // fluxPos, fluxNeg, null, null,
                computeCapa);

        //
        if (HexaedreDroit.isAdaptatif()) {
            if (meshDir.endsWith("SERP_SYM_VOL_5.DEC")) {
                solP.setPtsDeGaussInductifsVol(-1, 27, 27, -1, 32, 32, false);
                //
            } else if (meshDir.endsWith("SERP_SYM_VOL.DEC")) {
                solP.setPtsDeGaussInductifsVol(-1, 75, 75, -1, 32, 32, false);
                //
            } else {
                throw new InternalError();
            }
        } else {
            solP.setPtsDeGaussInductifsVol(27, 8, 8);
        }
//        solP.setPtsDeGaussInductifsVol(125, 64, 64);
//        solP.setPtsDeGaussCapacitifs(16, 16);
        solP.setPtsDeGaussCapacitifs(25, 25);
        solP.setFullAnalyticalP(true);
//        solP.setCheckNaN(true);

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
        /*
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
         /*/
        nbBranches = solP.getNbLignes();
        if (meshDir.endsWith("SERP_SYM_VOL.DEC")) {
            circuitPur.addSourceISimple(nbBranches, 3318, 4000, "Source", 1.0, 0);// 3
            // cct
//            circuitPur.addSourceUSimple(nbBranches + 1, 3647, 4327, "cct", 0.0, 0);
            //
            //
        } else if (meshDir.endsWith("SERP_SYM_VOL_5.DEC")) {
            circuitPur.addSourceISimple(nbBranches, 8802, 8762, "Source", 1.0, 0);// 5
            // cct
//            circuitPur.addSourceUSimple(nbBranches + 1, 5934, 5790, "cct", 0.0, 0);
        }
        //*/
        circuitPur.finSaisie();

        solP.setCircuitElectrique(circuitPur);

        FaceDeg1 FSd = solP.getFD1_Dielectric();
        FaceDeg1 FSc = solP.getFD1_Conductor();

        solP.setExportTopologie("D:");
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // Resolution       
        boolean openFile = true;
        /*
         double ftab[] = new double[100];
         ftab = compteFrequenciesLog(8, 9, ftab.length);
         /*/
        double[] ftab = getFrequencies();
        //*/
        double imp[] = new double[ftab.length];
        for (int iF = 0; iF < ftab.length; iF++) {
            solP.setVerbose(false);
            double f = ftab[iF];
            /*
             boolean HmatComp = true;
             solP.setCompression(HmatComp);// Enleve la com pression
             solP.setParamIterativeSolver(1, 10000, -1e-8, 1, -150);
             solP.setParamPreconditionner(1, 500, 0, new double[]{500, -3e-1, 500});
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

            try {
                // Save the datas.
                Ecriture saveRes = new Ecriture("D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/f" + f + "_RES.out");
                for (int i = 0; i < 2; i++) {
                    double[] tmp = new double[nDof];
                    res.row(i).get(tmp);
                    saveRes.ecrire(tmp, ',');
                    saveRes.aLaLigne();
                }
                saveRes.close();
            } catch (IOException ex) {
                Logger.getLogger(ZfSerpentSym.class.getName()).log(Level.SEVERE, null, ex);
            }

            System.out.println("");
            System.out.println("");
            System.out.println("I= " + ib[0][2 * nbBranches] + " + j* " + ib[0][2 * nbBranches + 1]);
            System.out.println("U= " + ib[1][2 * nbBranches] + " + j* " + ib[1][2 * nbBranches + 1]);
            imp[iF] = Math.hypot(ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]);
            System.out.println("|Z|= " + imp[iF]);
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

            String path_JrealD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/" + nDof + "f" + f + "_d_RE.msh";
            String path_JimD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/" + nDof + "f" + f + "_d_IM.msh";
            String path_JmodD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/" + nDof + "f" + f + "_d_MOD.msh";

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

            String path_JrealC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/" + nDof + "f" + f + "_c_RE.msh";
            String path_JimC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/" + nDof + "f" + f + "_c_IM.msh";
            String path_JmodC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/" + nDof + "f" + f + "_c_MOD.msh";

            exportJreal = new ExportGmshHdiv(FSc, path_JrealC);
            exportJreal.addQuantity(Jreal_c, "Jreal_c");
            exportJimag = new ExportGmshHdiv(FSc, path_JimC);
            exportJimag.addQuantity(Jimag_c, "Jimag_c");
            J = new ComplexFaceQuantity(FSc, res.submatrix(0, FSd.getActiveDofCount(), 2, FSc.getActiveDofCount()));
            exportJreal = new ExportGmshHdiv(FSc, path_JmodC);
            exportJreal.addQuantityExportMod(J, "Jmod_c");

            String path_Merge = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/" + nDof + "f" + f + "_MERGE.msh";
            try {
                Ecriture merge = new Ecriture(path_Merge);
                merge.ecrire("Merge '" + path_JrealC + "';\n");
                merge.ecrire("Merge '" + path_JimC + "';\n");
                merge.ecrire("Merge '" + path_JmodC + "';\n");

                merge.ecrire("Merge '" + path_JrealD + "';\n");
                merge.ecrire("Merge '" + path_JimD + "';\n");
                merge.ecrire("Merge '" + path_JmodD + "';\n");
                merge.close();

                if (openFile) {
                    File out_file = new File(path_Merge);
                    Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file);
                    openFile = false;
                }
            } catch (IOException ex) {
                Logger.getLogger(ZfSerpentSym.class.getName()).log(Level.SEVERE, null, ex);
            }
            for (int i = 0; i < iF + 1; i++) {
                System.out.println(ftab[i] + " ; " + imp[i]);
            }
        }
        try {
            Ecriture saveImp = new Ecriture("d:/z.out");
            for (int i = 0; i < ftab.length; i++) {
                saveImp.ecrire(ftab[i] + " ; " + imp[i] + "\n");
            }
            saveImp.close();
        } catch (IOException ex) {
            Logger.getLogger(ZfSerpentSym.class.getName()).log(Level.SEVERE, null, ex);
        }

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    public static double[] getFrequencies() {
        return new double[]{
            1.00E+06, 1.74E+06, 2.95E+06, 5.14E+06, 8.59E+06, 1.00E+07,
            1.01E+07, 1.03E+07, 1.04E+07, 1.06E+07, 1.08E+07, 1.09E+07, 1.11E+07, 1.12E+07,
            1.14E+07, 1.16E+07, 1.17E+07, 1.19E+07, 1.20E+07, 1.22E+07, 1.24E+07, 1.26E+07,
            1.29E+07, 1.31E+07, 1.33E+07, 1.35E+07, 1.37E+07, 1.39E+07, 1.41E+07, 1.43E+07,
            1.45E+07, 1.47E+07, 1.50E+07, 1.52E+07, 1.54E+07, 1.56E+07, 1.58E+07, 1.60E+07,
            1.63E+07, 1.65E+07, 1.68E+07, 1.71E+07, 1.74E+07, 1.76E+07, 1.79E+07, 1.82E+07,
            1.85E+07, 1.88E+07, 1.90E+07, 1.93E+07, 1.96E+07, 1.99E+07, 2.01E+07, 2.04E+07,
            2.07E+07, 2.10E+07, 2.13E+07, 2.17E+07, 2.20E+07, 2.24E+07, 2.28E+07, 2.31E+07,
            2.35E+07, 2.39E+07, 2.42E+07, 2.46E+07, 2.50E+07, 2.53E+07, 2.57E+07, 2.61E+07,
            2.64E+07, 2.68E+07, 2.72E+07, 2.75E+07, 2.79E+07, 2.83E+07, 2.88E+07, 2.93E+07,
            2.98E+07, 3.02E+07, 3.07E+07, 3.12E+07, 3.17E+07, 3.22E+07, 3.27E+07, 3.31E+07,
            3.36E+07, 3.41E+07, 3.46E+07, 3.51E+07, 3.56E+07, 3.60E+07, 3.65E+07, 3.71E+07,
            3.77E+07, 3.83E+07, 3.90E+07, 3.96E+07, 4.03E+07, 4.09E+07, 4.15E+07, 4.22E+07,
            4.28E+07, 4.35E+07, 4.41E+07, 4.47E+07, 4.54E+07, 4.60E+07, 4.66E+07, 4.73E+07,
            4.79E+07, 4.86E+07, 4.93E+07, 5.01E+07, 5.10E+07, 5.18E+07, 5.27E+07, 5.35E+07,
            5.43E+07, 5.52E+07, 5.60E+07, 5.69E+07, 5.77E+07, 5.85E+07, 5.94E+07, 6.02E+07,
            6.11E+07, 6.19E+07, 6.28E+07, 6.36E+07, 6.46E+07, 6.57E+07, 6.68E+07, 6.79E+07,
            6.90E+07, 7.01E+07, 7.12E+07, 7.23E+07, 7.34E+07, 7.46E+07, 7.57E+07, 7.68E+07,
            7.79E+07, 7.90E+07, 8.01E+07, 8.12E+07, 8.23E+07, 8.34E+07, 8.46E+07, 8.58E+07,
            8.73E+07, 8.88E+07, 9.02E+07, 9.17E+07, 9.32E+07, 9.46E+07, 9.61E+07, 9.76E+07,
            9.90E+07, 1.00E+08, 1.02E+08, 1.03E+08, 1.05E+08, 1.06E+08, 1.08E+08, 1.09E+08,
            1.11E+08, 1.12E+08, 1.14E+08, 1.16E+08, 1.18E+08, 1.20E+08, 1.22E+08, 1.24E+08,
            1.26E+08, 1.28E+08, 1.30E+08, 1.32E+08, 1.34E+08, 1.36E+08, 1.38E+08, 1.39E+08,
            1.41E+08, 1.43E+08, 1.45E+08, 1.47E+08, 1.49E+08, 1.52E+08, 1.55E+08, 1.57E+08,
            1.60E+08, 1.62E+08, 1.65E+08, 1.67E+08, 1.70E+08, 1.72E+08, 1.75E+08, 1.78E+08,
            1.80E+08, 1.83E+08, 1.85E+08, 1.88E+08, 1.90E+08, 1.93E+08, 1.96E+08, 1.99E+08,
            2.02E+08, 2.06E+08, 2.09E+08, 2.13E+08, 2.16E+08, 2.19E+08, 2.23E+08, 2.26E+08,
            2.29E+08, 2.33E+08, 2.36E+08, 2.40E+08, 2.43E+08, 2.46E+08, 2.50E+08, 2.53E+08,
            2.56E+08, 2.60E+08, 2.65E+08, 2.69E+08, 2.74E+08, 2.78E+08, 2.82E+08, 2.87E+08,
            2.91E+08, 2.96E+08, 3.00E+08, 3.05E+08, 3.09E+08, 3.14E+08, 3.18E+08, 3.22E+08,
            3.27E+08, 3.31E+08, 3.36E+08, 3.41E+08, 3.47E+08, 3.53E+08, 3.58E+08, 3.64E+08,
            3.70E+08, 3.76E+08, 3.82E+08, 3.88E+08, 3.94E+08, 4.00E+08, 4.05E+08, 4.11E+08,
            4.17E+08, 4.23E+08, 4.29E+08, 4.35E+08, 4.41E+08, 4.46E+08, 4.53E+08, 4.61E+08,
            4.69E+08, 4.76E+08, 4.84E+08, 4.92E+08, 5.00E+08, 5.07E+08, 5.15E+08, 5.23E+08,
            5.31E+08, 5.38E+08, 5.46E+08, 5.54E+08, 5.62E+08, 5.69E+08, 5.77E+08, 5.85E+08,
            5.94E+08, 6.04E+08, 6.14E+08, 6.24E+08, 6.34E+08, 6.45E+08, 6.55E+08, 6.65E+08,
            6.75E+08, 6.86E+08, 6.96E+08, 7.06E+08, 7.16E+08, 7.26E+08, 7.37E+08, 7.47E+08,
            7.57E+08, 7.67E+08, 7.77E+08, 7.89E+08, 8.03E+08, 8.16E+08, 8.30E+08, 8.43E+08,
            8.57E+08, 8.70E+08, 8.84E+08, 8.97E+08, 9.10E+08, 9.24E+08, 9.37E+08, 9.51E+08,
            9.64E+08, 9.78E+08, 9.91E+08, 1.00E+09};
    }

}
