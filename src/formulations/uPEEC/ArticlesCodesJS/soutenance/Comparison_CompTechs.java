/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.ArticlesCodesJS.soutenance;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.elements.ElementFactory;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D;
import g2elab.mipse.numericalTools.vector.full.VectorFull;
import g2elab.mipse.numericalTools.vector.full.VectorFullComplex;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.files.Exec;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * COMPARE THE
 *
 * @author jsiau
 */
public class Comparison_CompTechs {

    public static ImportFlux mesh;
    public static String meshChoice = "SERP_1EP";
//    private static String meshChoice = "SERP_3EP";

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(1);
        double f = 1e6;
//        solve(Compression.HCA, f);
        testProds(f);
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    public static void testProds(double f) {
        PEEC_RLMPC_Volume solFMM = getSolver(Compression.FMM, f);
        PEEC_RLMPC_Volume solFULL = getSolver(Compression.No, f);
        //
        VectorFullComplex x = VectorFullComplex.generateRandom(solFMM.getNbLignes());
        VectorFullComplex bfmm = new VectorFullComplex(solFMM.produit(x.getValues(true), new double[2 * x.length()], f));
        VectorFullComplex bfull = new VectorFullComplex(solFULL.produit(x.getValues(true), new double[2 * x.length()], f));
        //
        VectorFull dif = bfmm.sub(bfull, null);
        System.out.println("Erreur prod = " + dif.norm2() / bfull.norm2());

        System.out.println("--Ecriture de matFF--");
        Basic2D.getFullMatrix(solFMM.getFinalMatrix(), solFMM.getNbLignes(), solFMM.getNbLignes()).saveASCII("d:/tmp/matFMM.mat", ',');
        System.out.println("--Ecriture de matFF terminee--");
        System.out.println("--Ecriture de matFULL--");
        Basic2D.getFullMatrix(solFULL.getFinalMatrix(), solFULL.getNbLignes(), solFULL.getNbLignes()).saveASCII("d:/tmp/matFULL.mat", ',');
        System.out.println("--Ecriture de matFULL terminee--");
    }

    /**
     * LOAD THE MESH
     */
    public static void loadMesh() {
        ElementFactory ory = new ElementFactory(ElementFactory.STEP, ElementFactory.STEP, ElementFactory.STEP);
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Comparison_CompTechs.class.getName()).log(Level.SEVERE, null, ex);
        }
        mesh = new ImportFlux(meshDir + "/src/formulations/uPEEC/ArticlesCodesJS/soutenance/" + meshChoice + ".DEC");
        ory.addImport(mesh);
    }

    /**
     *
     * @param compression
     * @param f
     * @return
     */
    public static PEEC_RLMPC_Volume getSolver(Compression compression, double f) {
        if (mesh == null) {
            loadMesh();
        }
        ///////////// Gather the meshes
        // Conductor
        VolumeRegion conductor = (VolumeRegion) mesh.getRegion(0);
        // Put the  2 dieledctric regions in an array
        VolumeRegion dielectrics[] = new VolumeRegion[]{(VolumeRegion) mesh.getRegion(1)};
        // Interface between the conductor and the dielectrics
        SurfaceRegion interfaceCondDielec = (SurfaceRegion) mesh.getRegion(2);
        // And compute their union
        SurfaceRegion dielecAir = (SurfaceRegion) mesh.getRegion(3);
        SurfaceRegion fluxPos = null;
        SurfaceRegion fluxNeg = null;
        //
        PEEC_RLMPC_Volume solP = new PEEC_RLMPC_Volume(conductor, 1 / 1.68e-8,
                dielectrics, new double[]{4.7 * eps0},
                fluxPos, fluxNeg, interfaceCondDielec, dielecAir);
        //
        solP.setPtsDeGaussInductifsVol(27, 8, 8);
        solP.setPtsDeGaussCapacitifs(4, 4);
//            solP.setFullAnalyticalP(true);
        //
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int nbBranches = solP.getNbLignes();
        if (mesh.equals("SERP_3EP")) {
            circuitPur.addSourceISimple(nbBranches, 14351, 14271, "Source I", 1.0, 0.0);
        } else {

        }
        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);
        solP.setExportTopologie("D:");
        solP.setCompression(compression);
        //
        solP.integration(f);
        //
        return solP;
    }

    /**
     *
     * @param compression
     * @param f
     * @return
     */
    public static double[][] solve(Compression compression, double f) {
        // CONSTRUCTION OF THE SOLVER
        PEEC_RLMPC_Volume solP = getSolver(compression, f);
        // SOLUTION
        double ib[][] = solP.resolutionIterative(f);
        //
        exportRes(solP, ib, f, compression);
        //
        return ib;
    }

    public static void exportRes(PEEC_RLMPC_Volume solP, double ib[][], double f, Compression compression) {
        String outPath = "D:/jsiau/_Backup_Sources/Resultats/EtudeCompTech/" + compression;
        FaceDeg1 FSd = solP.getFD1_Dielectric();
        FaceDeg1 FSc = solP.getFD1_Conductor();
        int nbBranches = solP.getNbLignes();
        int nDof = FSd.getActiveDofCount() + FSc.getActiveDofCount();

        System.out.println("nDof= " + nDof);
        Matrix res = new Matrix(2, nDof);
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }
        try {
            // Save the datas.
            Ecriture saveRes = new Ecriture(outPath + "/f" + f + "_RES.out");
            for (int i = 0; i < 2; i++) {
                double[] tmp = new double[nDof];
                res.row(i).get(tmp);
                saveRes.ecrire(tmp, ',');
                saveRes.aLaLigne();
            }
            saveRes.close();
        } catch (IOException ex) {
            Logger.getLogger(Comparison_CompTechs.class.getName()).log(Level.SEVERE, null, ex);
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

        String path_JrealD = outPath + "/" + nDof + "f" + f + "_d_RE.msh";
        String path_JimD = outPath + "/" + nDof + "f" + f + "_d_IM.msh";
        String path_JmodD = outPath + "/" + nDof + "f" + f + "_d_MOD.msh";

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

        String path_JrealC = outPath + "/" + nDof + "f" + f + "_c_RE.msh";
        String path_JimC = outPath + "/" + nDof + "f" + f + "_c_IM.msh";
        String path_JmodC = outPath + "/" + nDof + "f" + f + "_c_MOD.msh";

        exportJreal = new ExportGmshHdiv(FSc, path_JrealC);
        exportJreal.addQuantity(Jreal_c, "Jreal_c");
        exportJimag = new ExportGmshHdiv(FSc, path_JimC);
        exportJimag.addQuantity(Jimag_c, "Jimag_c");
        J = new ComplexFaceQuantity(FSc, res.submatrix(0, FSd.getActiveDofCount(), 2, FSc.getActiveDofCount()));
        exportJreal = new ExportGmshHdiv(FSc, path_JmodC);
        exportJreal.addQuantityExportMod(J, "Jmod_c");

        String path_Merge = outPath + "/" + nDof + "f" + f + "_MERGE.msh";
        try {
            Ecriture merge = new Ecriture(path_Merge);
            merge.ecrire("Merge '" + path_JrealC + "';\n");
            merge.ecrire("Merge '" + path_JimC + "';\n");
            merge.ecrire("Merge '" + path_JmodC + "';\n");

            merge.ecrire("Merge '" + path_JrealD + "';\n");
            merge.ecrire("Merge '" + path_JimD + "';\n");
            merge.ecrire("Merge '" + path_JmodD + "';\n");
            merge.close();
        } catch (IOException ex) {
            Logger.getLogger(Comparison_CompTechs.class.getName()).log(Level.SEVERE, null, ex);
        }
        Exec.openFile(path_Merge);
    }

}
