/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF;
import formulations.uPEEC.condVol_dielecVol.LoopAntenna_WithDielec;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.IO.paraview.ExportVtkHdiv;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.LineRegion;
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

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv.eps0;

/**
 *
 * @author jsiau
 */
public class CAPA_PLAN {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.out.println("0= Volumique \t 1= Surfacique");

        Scanner sc = new Scanner(System.in);
        int form = sc.nextInt();

        CAPA_PLAN c = new CAPA_PLAN(form);

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    public CAPA_PLAN(int formulation) {
        switch (formulation) {
            case 0:
                this.vol();
                break;
            case 1:
                this.surf();
                break;
            default:
                throw new IllegalArgumentException();
        }
    }

    double epsilon = 4.7* eps0;
    double f = 100;
    double sigma = 1 / 2.836E-8;

    private void checkC(double Zmod, double f, double epsilon) {
        double Cth = epsilon * 25e-4 / 1e-3;
        double cZ = 1 / (2 * Math.PI * f * Zmod);
        System.out.println("Capa Plan Theorique= " + Cth);
        System.out.println("Capa Plan Calculee= " + cZ);
        System.out.println("Erreur relative= " + Math.abs(Cth - cZ) / Math.abs(Cth));
    }

    private void surf() {
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_SURF.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/COND_SURF.DEC";

        ImportFlux mesh = new ImportFlux(meshDir);
        int ind[] = new int[]{1, 4, 0, 2, 3};
        int formulation = 3;

        double ep = 2e-5;

        PEEC_RLMPC_SURF solP = new PEEC_RLMPC_SURF((SurfaceRegion) mesh.getRegion(ind[0]), (LineRegion) mesh.getRegion(ind[1]), sigma, ep, (VolumeRegion) mesh.getRegion(ind[2]), epsilon, (LineRegion) mesh.getRegion(ind[3]), (LineRegion) mesh.getRegion(ind[4]),
                formulation, false, false);
        solP.setPtsDeGaussInductifsVol(27, 8, 8);// Hexaedre
        solP.setPtsDeGaussInductifsSurf(25, 4, 4);// Quandrangle
        //        solP.setPtsDeGaussInductifsSurf(64, 144);// Quandrangle bugg
        solP.setPtsDeGaussCapacitifs(4, 4);

        //
        //
        //
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int indNodeInf = solP.getIndNodeInfinite();
        int B1[] = solP.getIndexNodeBorne1();
        int B2[] = solP.getIndexNodeBorne2();

        System.out.print("B1=[");
        for (int i = 0; i < B1.length; i++) {
            System.out.print(" " + B1[i]);
        }
        System.out.println("]");
        System.out.print("B2=[");
        for (int i = 0; i < B2.length; i++) {
            System.out.print(" " + B2[i]);
        }
        System.out.println("]");

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

        int nbBranches = solP.getNbLignes() + B1.length + B2.length;
        // Le nombre de branche correspond aussi a lindice de lalimentation !
        circuitPur.addSourceISimple(nbBranches, indB1, indB2, "Source I", 1.0, 0.0);

        circuitPur.finSaisie();

        System.out.println("circuitPur=\n" + circuitPur.toStringOrdonance());

        solP.setCircuitElectrique(circuitPur);

        FaceDeg1 FSd = solP.getFD1_Dielectric();
        FaceDeg1 FSc = solP.getFD1_Conductor();

        System.out.println("circuitPur.NbLigne= " + circuitPur.getNbLignes());
        System.out.println("solP.NbLigne= " + solP.getNbLignes());
        System.err.println("NbBranches= " + nbBranches);
//        solP.setExportTopologie("D:");
        // Resolution
        /*
         solP.setCompression(false);// Enleve la compression
         double ib[][] = solP.resolutionIterative(f);
         /*/
        solP.setCompression(false);// Enleve la compression
        double ib[][] = solP.resolutionDirecte(f);
        //*/
        int nDof = FSd.getActiveDofCount() + FSc.getActiveDofCount();

        System.out.println("nDof= " + nDof);
        Matrix res = new Matrix(2, nDof);
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }
        System.out.println("res=\n" + res);
        //* Scale seulement les DOF du conducteur
        res.submatrix(0, FSd.getActiveDofCount(), 2, FSc.getActiveDofCount()).scale(1 / ep);
        /*/ // Scale tous les dof
         res.scale(1/ep);
         //*/

        try {
            // Save the datas.
            Ecriture saveRes = new Ecriture("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Condensator/" + nDof + "f" + f + "_RES.out");
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

        System.out.println("I= " + ib[0][2 * nbBranches] + " + j* " + ib[0][2 * nbBranches + 1]);
        System.out.println("U= " + ib[1][2 * nbBranches] + " + j* " + ib[1][2 * nbBranches + 1]);

        System.out.println("");
        this.checkC(Math.hypot(ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]), f, epsilon);
        System.out.println("");

//        System.out.println("P:\n"+solP.getP());
        //
        //
        // EXPORT ALL THE QUANTITIES
        //
        //
        RealFaceQuantity Jreal_d = new RealFaceQuantity(FSd, res.row(0).subvector(0, FSd.getActiveDofCount()));
        RealFaceQuantity Jimag_d = new RealFaceQuantity(FSd, res.row(1).subvector(0, FSd.getActiveDofCount()));

        String path_JrealD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Condensator/" + nDof + "f" + f + "_d_RE.msh";
        String path_JimD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Condensator/" + nDof + "f" + f + "_d_IM.msh";
        String path_JmodD = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Condensator/" + nDof + "f" + f + "_d_MOD.msh";

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

        String path_JrealC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Condensator/" + nDof + "f" + f + "_c_RE.msh";
        String path_JimC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Condensator/" + nDof + "f" + f + "_c_IM.msh";
        String path_JmodC = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Condensator/" + nDof + "f" + f + "_c_MOD.msh";

        exportJreal = new ExportGmshHdiv(FSc, path_JrealC);
        exportJreal.addQuantity(Jreal_c, "Jreal_c");
        exportJimag = new ExportGmshHdiv(FSc, path_JimC);
        exportJimag.addQuantity(Jimag_c, "Jimag_c");
        J = new ComplexFaceQuantity(FSc, res.submatrix(0, FSd.getActiveDofCount(), 2, FSc.getActiveDofCount()));
        exportJreal = new ExportGmshHdiv(FSc, path_JmodC);
        exportJreal.addQuantityExportMod(J, "Jmod_c");

        String path_Merge = "D:/jsiau/_Backup_Sources/Resultats/Dielectric/Conductor_surface/Condensator/" + nDof + "f" + f + "_MERGE.msh";
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

        ExportVtkHdiv exV = new ExportVtkHdiv(FSd, "D:/test.vtk");

    }

    private void vol() {
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_SURF.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/COND_VOL16.DEC";

        ImportFlux Imp = new ImportFlux(meshDir);
        int nDof[] = new int[Imp.getRegions().length];
        for (int i = 0; i < nDof.length; i++) {
            nDof[i] = Imp.getRegion(i).getElementSet().getNbElement();
        }
        /*
         CONSTRUCTION DU SUPPORT DES CONDUCTIVITES EQUIVALENTES
         */
        /*
         VolumeRegion rV[] = new VolumeRegion[1];
         rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice

         double Sigma[][] = new double[2][nDof[0]];// Size = [2][nbElmttotal]
         // Conducteur
         for (int i = 0; i < nDof[0]; i++) {
         Sigma[0][i] = sigma;// Conductivite du cuivre
         Sigma[1][i] = eps0;// Epsilon0
         }
         /*/
        // Dielectrique
        VolumeRegion rV[] = new VolumeRegion[2];
        rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice
        rV[1] = (VolumeRegion) Imp.getRegion(1);// Region dielectrique

        double Sigma[][] = new double[2][nDof[0] + nDof[1]];// Size = [2][nbElmttotal]
        // Conducteur
        for (int i = 0; i < nDof[0]; i++) {
            Sigma[0][i] = sigma;// Conductivite du cuivre
            Sigma[1][i] = eps0;// Epsilon0
        }
        for (int i = nDof[0]; i < nDof[0] + nDof[1]; i++) {
            Sigma[0][i] = 0.0;// 0 pour le dielec
            Sigma[1][i] = epsilon;// Espilon
        }
         //*/

        // Create the solver
        PEEC_RLMPC_ConducEquiv solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma, (SurfaceRegion) Imp.getRegion(3), (SurfaceRegion) Imp.getRegion(4), true);
        solP.setPtsDeGaussInductifs(8, 8, 27);
        solP.setPtsDeGaussCapacitifs(4, 4);

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        FaceDeg1 Omega = solP.getFunctionSpace();
        System.out.println("Activedofcount= " + Omega.getActiveDofCount());

        int nbBranches = solP.getNbLignes();
        int C1 = Omega.getElementSet().getNbElement() + 1;
        int C2 = C1 + 1;
        int brancheSource = nbBranches;
        circuitPur.addSourceISimple(brancheSource, C1, C2, "SourceI", 1, 0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

        // Resolution
        /*
         double ib[][] = solP.resolutionIterative(f);
         /*/
        double ib[][] = solP.resolutionDirecte(f);
        //*/
        Matrix res = new Matrix(2, Omega.getActiveDofCount());
        for (int i = 0; i < Omega.getActiveDofCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }
        double[] IsurSource = new double[2];
        IsurSource[0] = ib[0][2 * brancheSource];
        IsurSource[1] = ib[0][2 * brancheSource + 1];

        double[] UsurSource = new double[2];
        UsurSource[0] = ib[1][2 * brancheSource];
        UsurSource[1] = ib[1][2 * brancheSource + 1];

        System.err.println("UsurSource= " + ib[1][2 * brancheSource] + " + j. " + ib[1][2 * brancheSource + 1]);
        System.err.println("IsurSource= " + IsurSource[0] + " + j. " + IsurSource[1]);

        System.out.println("");
        this.checkC(Math.hypot(UsurSource[0], UsurSource[1]), f, epsilon);
        System.out.println("");

        // On reorganise les dof   
        RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
        RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));

        ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/PP/S_RE_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJreal.addQuantity(Jreal, "Jreal");

        ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/PP/S_IM_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJimag.addQuantity(Jimag, "Jimag");

        ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
        exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/PP/S_Jmod" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJreal.addQuantityExportMod(J, "Jmod");

        File out_file = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/PP/S_Jmod" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        File out_file1 = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/PP/S_RE_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        File out_file2 = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/PP/S_IM_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        try {
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);

        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_ConducEquiv.class.getName()).log(Level.SEVERE, null, ex);
        }

        double b[] = new double[2 * solP.getNbLignes()];
        for (int i = 0; i < b.length; i++) {
            b[i] = Math.random() * 100;
        }

        Basic2D Mf = new Basic2D(solP.getZbFull(new double[solP.getNbLignes()][2 * solP.getNbLignes()], 0, solP.getNbLignes() - 1, 0, solP.getNbLignes() - 1, f));

        double rf[] = Mf.product(b, new double[2 * solP.getNbLignes()]);
        double pf[] = solP.produit(b, new double[2 * solP.getNbLignes()], f);

        ColumnVector vf = new ColumnVector(rf);
        ColumnVector vp = new ColumnVector(pf);

        vp.sub(vf);
        System.out.println("Error pmv's= " + vp.norm() / vf.norm());

    }

}
