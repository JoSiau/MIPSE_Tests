/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.condVol_dielecVol;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMesh;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.numericalTools.matrix.MatriceIncidence;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;
import g2elab.mipse.mipseCore.matrixCompression.Compression;

/**
 *
 * @author jsiau
 */
public class TestNewBobine {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_Volume.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportGmshMesh mesh = new ImportGmshMesh(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/bob3.msh");
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
        System.out.println("Nb Face Bob/Air= " + bobAir.getElementSet().getNbElement());
        //*
        SurfaceRegion formerAir = (SurfaceRegion) mesh.getRegion(3);
        System.out.println("Nb Face CoilFormer/Air= " + formerAir.getElementSet().getNbElement());
        /*/
         SurfaceRegion formerAir = (SurfaceRegion) mesh.getRegion(2).generateBorder();
         System.out.println("Nb Face CoilFormer/Air= " + formerAir.getElementSet().getNbElement());        
         //*/
        bobAir.exportGmsh("D:/tmp/bobAir.msh");
        formerAir.exportGmsh("D:/tmp/formerAir.msh");

        // And compute their union
        SurfaceRegion dielecAir = (SurfaceRegion) bobAir.generateUnionWith(formerAir);
        System.out.println("Nb Face Dielec/Air= " + dielecAir.getElementSet().getNbElement());
        dielecAir.exportGmsh("d:/dielecAir.msh");
        SurfaceRegion fluxPos = (SurfaceRegion) mesh.getRegion(4);
        SurfaceRegion fluxNeg = (SurfaceRegion) mesh.getRegion(5);

        double f = 8e6;//6454614.473;//5223345.074;
        boolean computeCapa = true;

        // Sigma = [ 5.6268e7 ; 6.0781e7 ] et  5.8502e7 la valeur nominale
        PEEC_RLMPC_Volume solP = new PEEC_RLMPC_Volume(conductor, 5.8502e7,
                dielectrics, new double[]{3.2 * eps0, 3.0 * eps0},
                fluxPos, fluxNeg, interfaceCondDielec, dielecAir,
                computeCapa);

        solP.setPtsDeGaussInductifsVol(27, 8, 8);
//        solP.setPtsDeGaussInductifsVol(125, 64, 64);
        solP.setPtsDeGaussCapacitifs(4, 4);
//        solP.setFullAnalyticalP(true);
        solP.setCheckNaN(true);

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

        MatriceIncidence mi = solP.getSolveurCircuit().getMI().transpose();
        mi.saveMTX( "Mi");

        FaceDeg1 FSd = solP.getFD1_Dielectric();
        FaceDeg1 FSc = solP.getFD1_Conductor();

        solP.setExportTopologie("D:");
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // Resolution
        //*
        boolean HmatComp = true;
        solP.setCompression(HmatComp?Compression.HCA:Compression.No);// Enleve la com pression
        solP.setParamIterativeSolver(1, 10000, -1e-8, 1, -150);
        solP.setParamPreconditionner(1, 500, 0, new double[]{500, -3e-1, 500});
        double ib[][] = solP.resolutionIterative(f);
        /*/
         boolean HmatComp = false;
         solP.setCompression(HmatComp?Compression.HCA:Compression.No);// Enleve la compression
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
            Ecriture saveRes = new Ecriture("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/f" + f + "_RES.out");
            for (int i = 0; i < 2; i++) {
                double[] tmp = new double[nDof];
                res.row(i).get(tmp);
                saveRes.ecrire(tmp, ',');
                saveRes.aLaLigne();
            }
            saveRes.close();
        } catch (IOException ex) {
            Logger.getLogger(LoopAntenna_WithDielec.class.getName()).log(Level.SEVERE, null, ex);
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

        String path_JrealD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/" + nDof + "f" + f + "_d_RE.msh";
        String path_JimD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/" + nDof + "f" + f + "_d_IM.msh";
        String path_JmodD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/" + nDof + "f" + f + "_d_MOD.msh";

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

        String path_JrealC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/" + nDof + "f" + f + "_c_RE.msh";
        String path_JimC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/" + nDof + "f" + f + "_c_IM.msh";
        String path_JmodC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/" + nDof + "f" + f + "_c_MOD.msh";

        exportJreal = new ExportGmshHdiv(FSc, path_JrealC);
        exportJreal.addQuantity(Jreal_c, "Jreal_c");
        exportJimag = new ExportGmshHdiv(FSc, path_JimC);
        exportJimag.addQuantity(Jimag_c, "Jimag_c");
        J = new ComplexFaceQuantity(FSc, res.submatrix(0, FSd.getActiveDofCount(), 2, FSc.getActiveDofCount()));
        exportJreal = new ExportGmshHdiv(FSc, path_JmodC);
        exportJreal.addQuantityExportMod(J, "Jmod_c");

        String path_Merge = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Bob_Zach/" + nDof + "f" + f + "_MERGE.msh";
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

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
