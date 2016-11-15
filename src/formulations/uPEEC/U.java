/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.File;
import java.util.Scanner;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv.eps0;

/**
 *
 * @author jsiau
 */
public class U {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.out.println("Entrez la fréquence: ");
        Scanner sc = new Scanner(System.in);
        double f = sc.nextDouble();

        File fil = new File("");
        String path = fil.getAbsolutePath();
        /*        
         Nombre de regions importees : 3
         Region 0, Nom : VOLUME, type : 3, 927 elements
         Region 1, Nom : B1, type : 2, 9 elements
         Region 2, Nom : B2, type : 2, 9 elements
         */
        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_AIDED20x5.DEC";
//        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_GROS.DEC";
        ImportFlux ImF = new ImportFlux(file);
        VolumeRegion vr[] = new VolumeRegion[]{(VolumeRegion) ImF.getRegion(0)};
        double[][] cond = new double[2][vr[0].getElementSet().getNbElement()];
        double sigma =  1/ 1.72e-8;
        for (int i = 0; i < vr[0].getElementSet().getNbElement(); i++) {
            cond[0][i] = sigma;//5.814e8;
            cond[1][i] = eps0;
        }

        PEEC_RLMPC_ConducEquiv solP = new PEEC_RLMPC_ConducEquiv(vr, cond, (SurfaceRegion) ImF.getRegion(1), (SurfaceRegion) ImF.getRegion(2), false);
        
        solP.setPtsDeGaussInductifs(8, 8, 27);
        solP.setPtsDeGaussCapacitifs(4, 4);

        FaceDeg1 Omega = solP.getFunctionSpace();
        System.out.println("Activedofcount= " + Omega.getActiveDofCount());

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int nbBranches = Omega.getActiveDofCount();
        int C1 = Omega.getElementSet().getNbElement() + 1;
        int C2 = C1 + 1;
        // === CREER LES SOURCES
        int branchRef = nbBranches;
        circuitPur.addSourceISimple(branchRef, C1, C2, "Source I", 1, 0);
        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);
        solP.setVerbose(false);

        double ib[][] = solP.resolutionDirecte(f);
        Matrix res = new Matrix(2, Omega.getActiveDofCount());
        for (int i = 0; i < Omega.getActiveDofCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }

        System.out.println("I= " + ib[0][2 * branchRef] + " + i* " + ib[0][2 * branchRef + 1]);
        System.out.println("U= " + ib[1][2 * branchRef] + " + i* " + ib[1][2 * branchRef + 1]);
        System.out.println("|U|= " + Math.hypot(ib[1][2 * branchRef], ib[1][2 * branchRef + 1]));

        RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
        RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));

        ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_re.msh");
        exportJreal.addQuantity(Jreal, "Jreal");

        ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_im.msh");
        exportJimag.addQuantity(Jimag, "Jimag");

        ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
        exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/DEBUGG/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_Jmod.msh");
        exportJreal.addQuantityExportMod(J, "Jmod");

//        File out_file = new File("D:/jsiau/_Backup_Sources/Resultats/DEBUGG/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_Jmod.msh");
//        File out_file1 = new File("D:/jsiau/_Backup_Sources/Resultats/DEBUGG/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_im.msh");
//        File out_file2 = new File("D:/jsiau/_Backup_Sources/Resultats/DEBUGG/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_re.msh");
//        try {
//            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);
//        } catch (IOException ex) {
//            Logger.getLogger(PEEC_RLMPC_DEBUG_FORMULATION.class.getName()).log(Level.SEVERE, null, ex);
//        }
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
