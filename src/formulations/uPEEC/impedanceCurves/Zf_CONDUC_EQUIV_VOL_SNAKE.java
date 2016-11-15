/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.impedanceCurves;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import formulations.uPEEC.condVol_dielecVol.LoopAntenna_WithDielec;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv;
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
import java.util.logging.Level;
import java.util.logging.Logger;

import static formulations.uPEEC.impedanceCurves.Zf_PEEC_RLMPC_SURF_SERPENT.compteFrequenciesLog;
import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv.eps0;

/**
 *
 * @author clusterdev
 */
public class Zf_CONDUC_EQUIV_VOL_SNAKE {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // Enable the adaptive integration
        HexaedreDroit.setAdaptatif(true);
        //
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
        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/SERPENTDIELEC.DEC";

        ImportFlux Imp = new ImportFlux(file);
        int nDof[] = new int[4];
        for (int i = 0; i < nDof.length; i++) {
            nDof[i] = Imp.getRegion(i).getElementSet().getNbElement();
        }
        VolumeRegion rV[] = new VolumeRegion[2];
        rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice
        rV[1] = (VolumeRegion) Imp.getRegion(1);// Region dielectrique

        /*
         CONSTRUCTION DU SUPPORT DES CONDUCTIVITES EQUIVALENTES
         */
        double Sigma[][] = new double[2][nDof[0] + nDof[1]];// Size = [2][nbElmttotal]
        // Conducteur
        double sigma = 1 / 1.68e-8;
        double Re = sigma;
        double Im = eps0;
        for (int i = 0; i < nDof[0]; i++) {
            Sigma[0][i] = Re;// Conductivite du cuivre
            Sigma[1][i] = Im;// Epsilon0
        }
        // Dielectrique
        for (int i = nDof[0]; i < nDof[0] + nDof[1]; i++) {
            Sigma[0][i] = 0;// 0 pour le dielec
            Sigma[1][i] = 4.7 * eps0;// Espilon
        }
        // Create the solver
        PEEC_RLMPC_ConducEquiv solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma,
                (SurfaceRegion) Imp.getRegion(2), (SurfaceRegion) Imp.getRegion(3),
                false);
        //
        //
        if (HexaedreDroit.isAdaptatif()) {
            solP.setPtsDeGaussInductifs(9, 9, 16);
        } else {
            solP.setPtsDeGaussInductifs(8, 8, 27);
        }
        solP.setAnalyticalIntegrationCapa(true);
        solP.setPtsDeGaussCapacitifs(4, 4);

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

        double ftab[] = new double[50];
        ftab = compteFrequenciesLog(6, 9, ftab.length);
        
        
        double imp[] = new double[ftab.length];
        for (int ifr = 0; ifr < ftab.length; ifr++) {
            double f = ftab[ifr];
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
                Ecriture saveRes = new Ecriture("D:/tmp/" + Omega.getActiveDofCount() + "f" + f + "_RES.out");
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

            double ir = ib[0][2 * nbBranches], ii = ib[0][2 * nbBranches + 1], ur = ib[1][2 * nbBranches], ui = ib[1][2 * nbBranches + 1];
            System.out.println("I= " + ir + " + j* " + ii);
            System.out.println("U= " + ur + " + j* " + ui);
            imp[ifr] = Math.hypot(ur, ui) / Math.hypot(ir, ii);
            System.out.println("Z(" + f + ") = " + imp[ifr]);
            // On reorganise les dof   
            RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
            RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));

            String fold = "D:/tmp/postProSnakeVol/";
            ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, fold + "f" + f + "_RE.msh");
            exportJreal.addQuantity(Jreal, "Jreal");
//        exportJreal.vizualiseFaceValue(Jreal);

            ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, fold + "f" + f + "_IM.msh");
            exportJimag.addQuantity(Jimag, "Jimag");
//        exportJimag.vizualiseFaceValue(Jimag);

            ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
            exportJreal = new ExportGmshHdiv(Omega, fold + "f" + f + "_MOD.msh");
            exportJreal.addQuantityExportMod(J, "Jmod");

            for (int i = 0; i < ifr; i++) {
                System.out.println("Z(" + ftab[i] + ") = " + imp[i]);
            }
        }
        try {
            Ecriture sa = new Ecriture("D:/tmp/Z.txt");
            sa.ecrire(imp, ' ');
            sa.close();
        } catch (IOException ex) {
            Logger.getLogger(Zf_CONDUC_EQUIV_VOL_SNAKE.class.getName()).log(Level.SEVERE, null, ex);
        }

        GestionnaireTaches.getGestionnaireTaches().stop();

    }

}
