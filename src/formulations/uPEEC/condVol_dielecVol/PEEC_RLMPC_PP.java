/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.condVol_dielecVol;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv.eps0;

/**
 * EXAMPLE OF USE OF THE PEEC_RLMPC FOR 2 PARALLEL PLATES.
 *
 * @author jsiau
 */
public class PEEC_RLMPC_PP {

    public static void main(String[] args) {
        File fil = new File("");
        String path = fil.getAbsolutePath();

        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/PP_D16.DEC";
        /*
         Nombre de regions importees : 4
         Region 0, Nom : VOLUME, type : 3, 200 elements
         Region 1, Nom : DIELECTRIC, type : 3, 200 elements
         Region 2, Nom : BORNE_POS, type : 2, 10 elements
         Region 3, Nom : BORNE_NEG, type : 2, 10 elements
         */
        ImportFlux Imp = new ImportFlux(file);
        int nDof[] = new int[Imp.getRegions().length];
        for (int i = 0; i < nDof.length; i++) {
            nDof[i] = Imp.getRegion(i).getElementSet().getNbElement();
        }

        double f =1e7;
        /*
         CONSTRUCTION DU SUPPORT DES CONDUCTIVITES EQUIVALENTES
         */
        /*
         VolumeRegion rV[] = new VolumeRegion[1];
         rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice
         double Sigma[][] = new double[2][nDof[0]];// + nDof[1]];// Size = [2][nbElmttotal]
         // Conducteur
         double Re = 1 / 2.836E-8;
         double Im = eps0;
         for (int i = 0; i < nDof[0]; i++) {
         Sigma[0][i] = Re;// Conductivite du cuivre
         Sigma[1][i] = Im;// Epsilon0
         }
         /*/ // Dielectrique
        VolumeRegion rV[] = new VolumeRegion[2];
        rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice
        rV[1] = (VolumeRegion) Imp.getRegion(1);// Region dielectrique
        double Sigma[][] = new double[2][nDof[0] + nDof[1]];// Size = [2][nbElmttotal]
        // Conducteur
        double Re = 1 / 1.68E-8;
        double Im = eps0;
        for (int i = 0; i < nDof[0]; i++) {
            Sigma[0][i] = Re;// Conductivite du cuivre
            Sigma[1][i] = Im;// Epsilon0
        }
        double epsilon = 4.7 * eps0;
        for (int i = nDof[0]; i < nDof[0] + nDof[1]; i++) {
            Sigma[0][i] = 0;// 0 pour le dielec
            Sigma[1][i] = epsilon;// Espilon
        }
        //*/
        // Create the solver
        PEEC_RLMPC_ConducEquiv solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma, (SurfaceRegion) Imp.getRegion(2), (SurfaceRegion) Imp.getRegion(3), true);
        solP.setPtsDeGaussInductifs(512, 27, 27);
        solP.setPtsDeGaussCapacitifs(225, 225);

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
        double C = eps0 * 1e-4 / 1e-3;
        System.err.println("C = eps0 * S / e = " + C);

        double omega = 2 * Math.PI * f;
        double Cm = 1 / UsurSource[1] / omega;
        System.err.println("C1 = Im(I) / omega = " + Cm);// U = 1.0
        double c2 = 1 / Math.hypot(UsurSource[0], UsurSource[1]) / omega;
        System.err.println("C2 = |I| / omega = " + c2);// U = 1.0

        System.err.println("Erreur relative C1= " + Math.abs(C - Cm) / Math.abs(C));
        System.err.println("Erreur relative C2= " + Math.abs(C - c2) / Math.abs(C));
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

        double x[] = new double[2 * solP.getNbLignes()];
        for (int i = 0; i < x.length; i++) {
            x[i] = Math.random();
        }

        double zb[][] = new double[solP.getNbLignes()][2 * solP.getNbLignes()];
        solP.getZbFull(zb, 0, solP.getNbLignes() - 1, 0, solP.getNbLignes() - 1, f);
        Basic2D Z = new Basic2D(zb);
        double vf[] = Z.product(x, new double[2 * solP.getNbLignes()]);
        ColumnVector rvf = new ColumnVector(vf);

        double vp[] = solP.produit(x, new double[2 * solP.getNbLignes()], f);
        ColumnVector rvp = new ColumnVector(vp);

        rvp.sub(rvf);
        System.out.println("Erreur abs= " + rvp.norm());
        System.out.println("erreur relative= " + rvp.norm() / rvf.norm());

        GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
