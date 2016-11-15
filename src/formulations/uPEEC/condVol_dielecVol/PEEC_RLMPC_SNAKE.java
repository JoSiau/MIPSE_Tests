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
 * EXAMPLE OF USE OF THE PEEC_RLMPC.
 *
 * @author jsiau
 */
public class PEEC_RLMPC_SNAKE {

    public static void main(String[] args) throws IOException {
        String defaultPath = null;
        if (args.length != 0) {
            defaultPath = args[0] + args[1];
            new File(defaultPath).mkdirs();
        }
        File fil = new File("");
        String path = fil.getAbsolutePath();

        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/SERPENTDIELEC_3C.DEC";
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
        VolumeRegion rV[] = new VolumeRegion[2];
        rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice
        rV[1] = (VolumeRegion) Imp.getRegion(1);// Region dielectrique

        double f = 1e6;//3.77E+07;
        /*
         CONSTRUCTION DU SUPPORT DES CONDUCTIVITES EQUIVALENTES
         */
        double Sigma[][] = new double[2][nDof[0] + nDof[1]];// Size = [2][nbElmttotal]
        // Conducteur
        double Re = 5.814e8;
        double Im = eps0;
        for (int i = 0; i < nDof[0]; i++) {
            Sigma[0][i] = Re;// Conductivite du cuivre
            Sigma[1][i] = Im;// Epsilon0
        }
        // Dielectrique
        Re = 0;
        Im = 4.7 * eps0;
        for (int i = nDof[0]; i < nDof[0] + nDof[1]; i++) {
            Sigma[0][i] = Re;// 0 pour le dielec
            Sigma[1][i] = Im;// Espilon
        }

        // Create the solver
        PEEC_RLMPC_ConducEquiv solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma);
        solP.setPtsDeGaussInductifs(8, 1);
        solP.setPtsDeGaussCapacitifs(4, 4);
//        solP.setMultiThread(true);
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        FaceDeg1 Omega = solP.getFunctionSpace();
        System.out.println("Activedofcount= " + Omega.getActiveDofCount());

        int nbBranches = Omega.getActiveDofCount();
        int C1 = Omega.getElementSet().getNbElement() + 1;
        int C2 = C1 + 1;
        /*
         CREER LES SOURCES
         */
        circuitPur.addSourceISimple(nbBranches, 761, 266, "Source I", 1, 0);
        /*
         RELIE LES PLAQUES
         */
//        solP.setCourtCircuit(12, 3);
        circuitPur.addSourceUSimple(nbBranches + 1, 12, 3, "cct", 0, 0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);
//solP.setMultiThread(false);
//        ElementSetHomogene esh =  new ElementSetHomogene( new Element[]{Omega.getElementSet().getElements(4)});
//        ExportGmsh exg = new ExportGmsh(esh, (defaultPath != null ? defaultPath : "Y:")+"/el4.msh");
//        esh =  new ElementSetHomogene( new Element[]{Omega.getElementSet().getElements(13)});
//        exg = new ExportGmsh(esh, (defaultPath != null ? defaultPath : "Y:")+"/el13.msh");
//        esh =  new ElementSetHomogene( new Element[]{Omega.getElementSet().getElements(762)});
//        exg = new ExportGmsh(esh, (defaultPath != null ? defaultPath : "Y:")+"/el762.msh");
//        esh =  new ElementSetHomogene( new Element[]{Omega.getElementSet().getElements(267)});
//        exg = new ExportGmsh(esh, (defaultPath != null ? defaultPath : "Y:")+"/el267.msh");
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
        System.out.println("Z(" + f + ")= " + ib[1][ib[1].length - 2] + " + j* " + ib[1][ib[1].length - 1]);
        System.out.println("|Z(" + f + ")|= " + Math.hypot(ib[1][ib[1].length - 2], ib[1][ib[1].length - 1]));
        // On reorganise les dof   
        RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
        RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));

        ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, (defaultPath != null ? defaultPath : "Y:") + "/Resultats/Dielectric/SNAKE/S_RE_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJreal.addQuantity(Jreal, "Jreal");

        ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, (defaultPath != null ? defaultPath : "Y:") + "/Resultats/Dielectric/SNAKE/S_IM_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJimag.addQuantity(Jimag, "Jimag");

        ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
        exportJreal = new ExportGmshHdiv(Omega, (defaultPath != null ? defaultPath : "Y:") + "/Resultats/Dielectric/SNAKE/S_Jmod" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJreal.addQuantityExportMod(J, "Jmod");

        File out_file = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources/Resultats/Dielectric") + "/SNAKE/S_Jmod" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        File out_file1 = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources/Resultats/Dielectric") + "/SNAKE/S_RE_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        File out_file2 = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources/Resultats/Dielectric") + "/SNAKE/S_IM_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        try {
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);

        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_ConducEquiv.class.getName()).log(Level.SEVERE, null, ex);
        }
        //GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
