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
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv.eps0;

/**
 *
 * @author jsiau
 */
public class LoopAntenna_WithDielec {

    public static void main(String[] args) {
        File fil = new File("");
        String path = fil.getAbsolutePath();
        /*
         Nombre de regions importees : 4
         Region 0, Nom : LOOP_ANTENNA, type : 3, 1305 elements
         Region 1, Nom : DIELECTRIC, type : 3, 2354 elements
         Region 2, Nom : POSITIVE, type : 2, 4 elements
         Region 3, Nom : NEGATIVE, type : 2, 4 elements
         */
//        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/LOOPANTENNA_SMALL.DEC";
        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/LOOPANTENNA_DIELEC_VOL.DEC";

        ImportFlux Imp = new ImportFlux(file);
        int nDof[] = new int[4];
        for (int i = 0; i < nDof.length; i++) {
            nDof[i] = Imp.getRegion(i).getElementSet().getNbElement();
        }
        VolumeRegion rV[] = new VolumeRegion[2];
        rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice
        rV[1] = (VolumeRegion) Imp.getRegion(1);// Region dielectrique

        double f = 10e9;
        /*
         CONSTRUCTION DU SUPPORT DES CONDUCTIVITES EQUIVALENTES
         */
        double Sigma[][] = new double[2][nDof[0] + nDof[1]];// Size = [2][nbElmttotal]
        // Conducteur
        double Re = 1/ 1.68e-8;
        double Im = eps0;
        for (int i = 0; i < nDof[0]; i++) {
            Sigma[0][i] = Re;// Conductivite du cuivre
            Sigma[1][i] = Im;// Epsilon0
        }
        // Dielectrique
        for (int i = nDof[0]; i < nDof[0] + nDof[1]; i++) {
            Sigma[0][i] = 0;// 0 pour le dielec
            Sigma[1][i] = 5 * eps0;// Espilon
        }
        // Create the solver
        PEEC_RLMPC_ConducEquiv solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma, (SurfaceRegion) Imp.getRegion(2), (SurfaceRegion) Imp.getRegion(3), true);
        solP.setPtsDeGaussInductifs(15, 4);
        solP.setPtsDeGaussCapacitifs(3, 3);

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        FaceDeg1 Omega = solP.getFunctionSpace();
        System.out.println("Activedofcount= " + Omega.getActiveDofCount());

        /*
         CREATION DES BORNES
         */
        int nbBranches = solP.getNbLignes();
        int C1 = Omega.getNbElement() + 1;
        int C2 = C1 + 1;

        circuitPur.addSourceISimple(nbBranches, C1, C2, "Source I", 1, 0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);
        
        solP.setMultiThread(true);

        // Resolution
//        double ib[][] = solP.resolutionIterative(f);
        double ib[][] = solP.resolutionDirecte(f);
        Matrix res = new Matrix(2, Omega.getActiveDofCount());
        for (int i = 0; i < Omega.getActiveDofCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }
        try {
            // Save the datas.
            Ecriture saveRes = new Ecriture("D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/"+ Omega.getActiveDofCount() + "f" + f + "_RES.out");
            for (int i = 0; i < 2; i++) {
                double[] tmp = new double[Omega.getActiveDofCount()];
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

        // On reorganise les dof   
        RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
        RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));

        ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/" + Omega.getActiveDofCount() + "f" + f + "_RE.msh");
        exportJreal.addQuantity(Jreal, "Jreal");
//        exportJreal.vizualiseFaceValue(Jreal);

        ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/" + Omega.getActiveDofCount() + "f" + f + "_IM.msh");
        exportJimag.addQuantity(Jimag, "Jimag");
//        exportJimag.vizualiseFaceValue(Jimag);

        ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
        exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/" + Omega.getActiveDofCount() + "f" + f + "_MOD.msh");
        exportJreal.addQuantityExportMod(J, "Jmod");

//        ExportVtkHdiv tmp = new ExportVtkHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/Jmod.vtk");
//        tmp.addQuantity(J, "J");

        File out_file = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/" + Omega.getActiveDofCount() + "f" + f + "_RE.msh");
        File out_file1 = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/" + Omega.getActiveDofCount() + "f" + f + "_IM.msh");
        File out_file2 = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/" + Omega.getActiveDofCount() + "f" + f + "_MOD.msh");
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
