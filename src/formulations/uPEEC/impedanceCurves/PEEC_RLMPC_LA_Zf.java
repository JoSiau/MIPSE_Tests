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
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * USED TO COMPUTE THE Z(F) FOR THIS SNAKE
 *
 * @author jsiau
 */
public class PEEC_RLMPC_LA_Zf {

    public static void main(String[] args) throws IOException {
        String defaultPath = null;
        if (args.length != 0) {
            defaultPath = args[0] + args[1];
            new File(defaultPath).mkdirs();
        }
        File fil = new File("");
        String path = fil.getAbsolutePath();

        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/LOOPANTENNA_DIELEC_VOL.DEC";
        /*
         Nombre de regions importees : 6
         Region 0, Nom : SNAKE, type : 3, 490 elements
         Region 1, Nom : DIELEC, type : 3, 476 elements
         Region 2, Nom : NORD_OUEST, type : 2, 2 elements
         Region 3, Nom : SUD_OUEST, type : 2, 2 elements
         Region 4, Nom : NORD_EST, type : 2, 2 elements
         Region 5, Nom : SUD_EST, type : 2, 2 elements
         */
        ImportFlux IF = new ImportFlux(file);
        int nDof[] = new int[IF.getRegions().length];
        for (int i = 0; i < nDof.length; i++) {
            nDof[i] = IF.getRegion(i).getElementSet().getNbElement();
        }
        VolumeRegion rV[] = new VolumeRegion[01];
        double Sigma[][] = new double[2][nDof[0]];// Size = [2][nbElmttotal]
        // Conducteur
        double Re = 1 / 2.836e-8;
        double Im = PEEC_RLMPC_ConducEquiv.eps0;
        for (int i = 0; i < nDof[0]; i++) {
            Sigma[0][i] = Re;// Conductivite du cuivre
            Sigma[1][i] = Im;// Epsilon0
        }

        double f[];
        int nbF;
        /*
         Lecture fic = null;
         try {
         fic = new Lecture("d:/f.txt");
         } catch (IOException ex) {
         Logger.getLogger(PEEC_RLMPC_LA_Zf.class.getName()).log(Level.SEVERE, null, ex);
         }
         nbF = fic.getNbLignes();
         System.out.println("Nombre de frequences a analyser= " + nbF);
         f = new double[nbF];
         for (int i = 0; i < nbF; i++) {
         f[i] = fic.lireLigneDouble();
         }
         /*/
        nbF = 100;
        double f_dep = 1e9;
        double f_fin = 1e10;
        double pas = (f_fin - f_dep) / (nbF - 1);
        System.out.println("pas= " + pas);
        f = new double[nbF];
        for (int i = 0; i < nbF; i++) {
            f[i] = f_dep + i * pas;
        }
        //*/
        Ecriture save = new Ecriture("D:/Zf_LA_rlmpc.txt");
        double ib[][];
        Matrix res;
        PEEC_RLMPC_ConducEquiv solP;
        double Z[] = new double[nbF];
//        VolumeRegion rV[] = new VolumeRegion[1];
        for (int j = 0; j < nbF; j += 1) {
            System.out.println(j + "/" + nbF + "f= " + f[j]);

            IF = new ImportFlux(file);

            rV[0] = (VolumeRegion) IF.getRegion(0);// Region conductrice
            // Create the solver
            solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma, (SurfaceRegion) IF.getRegion(2), (SurfaceRegion) IF.getRegion(3));
            solP.setPtsDeGaussInductifs(15, 4);
            solP.setPtsDeGaussCapacitifs(4, 4);

            BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
            FaceDeg1 Omega = solP.getFunctionSpace();
            System.out.println("Activedofcount= " + Omega.getActiveDofCount());

            int nbBranches = solP.getNbLignes();
            /*
             CREER LES SOURCES
             */
            int branchRef = nbBranches;
            circuitPur.addSourceISimple(nbBranches, Omega.getNbElement() + 1, Omega.getNbElement() + 2, "Source I", 1, 0);

            circuitPur.finSaisie();
            solP.setCircuitElectrique(circuitPur);
            solP.setVerbose(false);
            // Resolution
            /*
             double ib[][] = solP.resolutionIterative(f);
             /*/
            ib = solP.resolutionDirecte(f[j]);
            //*/
            for (int i = 0; i < j; i++) {
                System.out.println(f[i]+" "+Z[i]);
            }
            Z[j] = Math.hypot(ib[1][2 * branchRef], ib[1][2 * branchRef + 1]);
            System.err.println("|Z|= " + Z[j]);
            System.err.println("Z(" + f[j] + ")= " + ib[1][2 * branchRef] + " + j*" + ib[1][2 * branchRef + 1]);
            System.err.println("I= " + ib[0][2 * branchRef] + " + j*" + ib[0][2 * branchRef + 1]);

            save.ecrire(f[j] + " " + Z[j] + "\n");

            if (j == 0) {
                res = new Matrix(2, Omega.getActiveDofCount());
                for (int i = 0; i < Omega.getActiveDofCount(); i++) {
                    res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
                    res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
                }
                // On reorganise les dof   
                RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
                RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));

                ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/S_RE_f" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
                exportJreal.addQuantity(Jreal, "Jreal");

                ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/S_IM_f" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
                exportJimag.addQuantity(Jimag, "Jimag");

                ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
                exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/S_Jmod" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
                exportJreal.addQuantityExportMod(J, "Jmod");

                File out_file = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources") + "/LoopAntenna/S_Jmod" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
                File out_file1 = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources") + "/LoopAntenna/S_RE_f" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
                File out_file2 = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources") + "/LoopAntenna/S_IM_f" + f[j] + "_nDof" + Omega.getActiveDofCount() + ".msh");
                try {
                    Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);
                } catch (IOException ex) {
                    Logger.getLogger(PEEC_RLMPC_ConducEquiv.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        save.close();
        GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
