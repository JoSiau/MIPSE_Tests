/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.impedanceCurves;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.files.Lecture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv.eps0;

/**
 * USED TO COMPUTE THE Z(F) FOR THIS SNAKE
 *
 * @author jsiau
 */
public class PEEC_RLMPC_SNAKE_Zf1 {

    public static void main(String[] args) throws IOException {
        File fil = new File("");
        String path = fil.getAbsolutePath();

        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/SERPENT_3C.DEC";
        /*
         Nombre de regions importees : 6
         Region 0, Nom : SNAKE, type : 3, 490 elements
         Region 1, Nom : DIELEC, type : 3, 476 elements
         Region 2, Nom : NORD_OUEST, type : 2, 2 elements
         Region 3, Nom : SUD_OUEST, type : 2, 2 elements
         Region 4, Nom : NORD_EST, type : 2, 2 elements
         Region 5, Nom : SUD_EST, type : 2, 2 elements
         */
        ImportFlux Imp = new ImportFlux(file);
        int nDof[] = new int[Imp.getRegions().length];
        for (int i = 0; i < nDof.length; i++) {
            nDof[i] = Imp.getRegion(i).getElementSet().getNbElement();
        }
//        VolumeRegion rV[] = new VolumeRegion[2];
//        rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice
//        rV[1] = (VolumeRegion) Imp.getRegion(1);// Region dielectrique

        /*
         CONSTRUCTION DU SUPPORT DES CONDUCTIVITES EQUIVALENTES
         */
        double Sigma[][] = new double[2][nDof[0] + nDof[1]];// Size = [2][nbElmttotal]
        // Conducteur
        double Re = 5.814 * 1e8;
        double Im = 4.7*eps0;
        for (int i = 0; i < nDof[0]; i++) {
            Sigma[0][i] = Re;// Conductivite du cuivre
            Sigma[1][i] = Im;// Epsilon0
        }
        // Dielectrique
//        Re = 0;
//        Im = 4.7 * eps0;
//        for (int i = nDof[0]; i < nDof[0] + nDof[1]; i++) {
//            Sigma[0][i] = Re;// 0 pour le dielec
//            Sigma[1][i] = Im;// Espilon
//        }
        double f[];
        int nbF;
        //*
        Lecture fic = null;
        try {
            fic = new Lecture("d:/f.txt");
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_SNAKE_Zf1.class.getName()).log(Level.SEVERE, null, ex);
        }
        nbF = fic.getNbLignes();
        System.out.println("Nombre de frequences a analyser= " + nbF);
        f = new double[nbF];
        for (int i = 0; i < nbF; i++) {
            f[i] = fic.lireLigneDouble();
        }
        /*/
         nbF = 100;
         double f_dep = 1e6;
         double f_fin = 1e9;
         double pas = (f_fin - f_dep) / nbF;
         System.out.println("pas= "+pas);
         f = new double[nbF];
         for (int i = 0; i < nbF; i++) {
         f[i] = f_dep + i * pas;
         }
         //*/
        Ecriture save = new Ecriture("D:/Zf_snake_rlmpc_noye.txt");
        double ib[][];
        Matrix res;
        PEEC_RLMPC_ConducEquiv solP;
        VolumeRegion rV[] = new VolumeRegion[1];
        int nPas = 20;
        for (int j = 0; j < nbF; j += nPas) {
            System.out.println("f= " + f[j]);
            if (f[j] > 1e8) {
                nPas = 6;
            }
            Imp = new ImportFlux(file);

            rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice
//            rV[1] = (VolumeRegion) Imp.getRegion(1);// Region dielectrique
            // Create the solver
            solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma);
            solP.setPtsDeGaussInductifs(8, 1);
            solP.setPtsDeGaussCapacitifs(4, 4);

            BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
            FaceDeg1 Omega = solP.getFunctionSpace();
            System.out.println("Activedofcount= " + Omega.getActiveDofCount());


            int nbBranches = Omega.getActiveDofCount();
            int C1 = Omega.getElementSet().getNbElement() + 1;
            int C2 = C1 + 1;
            /*
             CREER LES SOURCES
             */
            //*
            int branchRef = nbBranches;
            circuitPur.addSourceISimple(branchRef, C1, C2, "Source I", 1, 0);
            /*/
             int branchRef = nbBranches;
             circuitPur.addSourceISimple(nbBranches, 761, 266, "Source I", 1, 0);
             //*/
            /*
             RELIE LES PLAQUES
             */
            /*
             int[] borne3 = getIndicesBornes(Omega, (SurfaceRegion) Imp.getRegion(2));
             int[] borne4 = getIndicesBornes(Omega, (SurfaceRegion) Imp.getRegion(3));
             for (int i = 0; i < borne3.length; i++) {
             circuitPur.addSourceUSimple(nbBranches + borne1.length + borne2.length + 1 + i, borne3[i], borne4[i], "cable " + i, 0, 0);
             }
             /*/
//            circuitPur.addSourceUSimple(nbBranches + 1 , 12, 3, "cct", 0, 0);
            //*/
            circuitPur.finSaisie();
            solP.setCircuitElectrique(circuitPur);
            solP.setVerbose(false);
            // Resolution
            /*
             double ib[][] = solP.resolutionIterative(f);
             /*/
            ib = solP.resolutionDirecte(f[j]);
            //*/
            res = new Matrix(2, Omega.getActiveDofCount());
            for (int i = 0; i < Omega.getActiveDofCount(); i++) {
                res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
                res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
            }

            System.err.println("Z(" + f[j] + ")= " + ib[1][2 * branchRef] + " + j*" + ib[1][2 * branchRef + 1]);
            System.err.println("|Z|= " + Math.hypot(ib[1][2 * branchRef], ib[1][2 * branchRef + 1]));
            System.err.println("I= " + ib[0][2 * branchRef] + " + j*" + ib[0][2 * branchRef + 1]);

            save.ecrire(f[j] + " " + Math.hypot(ib[1][2 * branchRef], ib[1][2 * branchRef + 1]) + "\n");

            // On reorganise les dof   
            RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
            RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));

            ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/S_RE_f" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
            exportJreal.addQuantity(Jreal, "Jreal");

            ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/S_IM_f" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
            exportJimag.addQuantity(Jimag, "Jimag");

            ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
            exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/S_Jmod" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
            exportJreal.addQuantityExportMod(J, "Jmod");

            if (j == 0) {
                File out_file = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/S_Jmod" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
                File out_file1 = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/S_RE_f" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
                File out_file2 = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/SNAKE/S_IM_f" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
                try {
                    Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);
                } catch (IOException ex) {
                    Logger.getLogger(PEEC_RLMPC_ConducEquiv.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
            if (f[j + 1] > 1e9) {
                save.close();
                GestionnaireTaches.getGestionnaireTaches().stop();
                return;
            }
        }
        save.close();
        GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
