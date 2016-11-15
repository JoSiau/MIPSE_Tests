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
import g2elab.mipse.meshCore.region.VolumeRegion;
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
public class PEEC_RLMPC_WO {

    public static void main(String[] args) throws IOException {
        String defaultPath = null;
        if (args.length != 0) {
            defaultPath = args[0] + args[1];
            new File(defaultPath).mkdirs();
        }
        File fil = new File("");
        String path = fil.getAbsolutePath();
        /*
         Nombre de regions importees : 4
         Region 0, Nom : DIELEC, type : 3, 1035 elements
         Region 1, Nom : CONDUCTOR, type : 3, 3240 elements
         Region 2, Nom : BORNE1, type : 2, 76 elements
         Region 3, Nom : BORNE2, type : 2, 76 elements
         */
        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/WIRE_OPENED_PARFAIT_4C_TINY.DEC";

        ImportFlux Imp = new ImportFlux(file);
        int nDof[] = new int[Imp.getRegions().length];
        for (int i = 0; i < nDof.length; i++) {
            nDof[i] = Imp.getRegion(i).getElementSet().getNbElement();
        }
        VolumeRegion rV[] = new VolumeRegion[2];
        rV[0] = (VolumeRegion) Imp.getRegion(1);// Region conductrice
        rV[1] = (VolumeRegion) Imp.getRegion(0);// Region dielectrique

        double f = 1e6;
        /*
         CONSTRUCTION DU SUPPORT DES CONDUCTIVITES EQUIVALENTES
         */
        double Sigma[][] = new double[2][nDof[0] + nDof[1]];// Size = [2][nbElmttotal]
        // Conducteur
        double Re = 59.6 * 1e6;
        double Im = eps0;
        for (int i = 0; i < nDof[1]; i++) {
            Sigma[0][i] = Re;// Conductivite du cuivre
            Sigma[1][i] = Im;// Epsilon0
        }
        // Dielectrique
        Re = 0;
        Im = 2.25;
        for (int i = nDof[1]; i < nDof[0] + nDof[1]; i++) {
            Sigma[0][i] = Re;// 0 pour le dielec
            Sigma[1][i] = Im;// Espilon
        }

        // Create the solver
        PEEC_RLMPC_ConducEquiv solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma);
        solP.setPtsDeGaussInductifs(15, 15);
        solP.setPtsDeGaussCapacitifs(3, 3);

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        FaceDeg1 Omega = solP.getFunctionSpace();
        System.out.println("Activedofcount= " + Omega.getActiveDofCount());
        int nbBranches = Omega.getActiveDofCount();
        int C1 = Omega.getElementSet().getNbElement() + 1;
        int C2 = C1 + 1;

        circuitPur.addSourceUSimple(nbBranches, C1, C2, "Source I", 1, 0);

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

        // On reorganise les dof   
        RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
        RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));

        ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, (defaultPath != null ? defaultPath : "D:") + "/Resultats/Dielectric/WIRE_OPENED/WORE_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJreal.addQuantity(Jreal, "Jreal");

        ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, (defaultPath != null ? defaultPath : "D:") + "/Resultats/Dielectric/WIRE_OPENED/WOIM_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJimag.addQuantity(Jimag, "Jimag");

        ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
        exportJreal = new ExportGmshHdiv(Omega, (defaultPath != null ? defaultPath : "D:") + "/Resultats/Dielectric/WIRE_OPENED/WOJmod" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJreal.addQuantityExportMod(J, "Jmod");

        File out_file = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources/Resultats/Dielectric/WIRE_OPENED") + "/WOJmod" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        File out_file1 = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources/Resultats/Dielectric/WIRE_OPENED") + "/WORE_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        File out_file2 = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources/Resultats/Dielectric/WIRE_OPENED") + "/WOIM_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        try {
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);

        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_ConducEquiv.class.getName()).log(Level.SEVERE, null, ex);
        }
        //GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
