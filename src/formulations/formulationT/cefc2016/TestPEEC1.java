/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT.cefc2016;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;
import java.io.File;

/**
 *
 * @author jsiau
 */
public class TestPEEC1 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        double f = Test_TIm.f;

        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + "/src/formulations/formulationT/PLAQUE_TROU_TETRA_4350.DEC";
//        path = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/PLAQUE_VOL.DEC";
//        path = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/LOOPANTENA.DEC";
//       
        System.out.println("Fichier maillage : " + path);
        ImportFlux mesh = new ImportFlux(path);

        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((VolumeRegion) mesh.getRegion(0), 1 / 1.68e-8, 1.0,
                null, null);
        //*/
        solP.setComputeCapa(false);
        solP.setPtsDeGaussInductifs(15, 4, 4);
//        solP.setPtsDeGaussCapacitifs(4, 4);
//        solP.setAnalyticalIntegrationCapa(true);
//        solP.setCorrectionVoisin(true);
        solP.setExportTopologie("d:/tmp");
        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        circuitPur.addSourceISimple(1000000, 6798, 7377, "source", 1.0, 0.0);
        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

        /*
        //        double[][] ib = solP.resolutionDirecte(f);
        solP.setHmatrixCompression(false);
        double[][] ib = solP.resolutionDirectePureConductor(f);
        /*/
         solP.setHmatrixCompression(false);
         double[][] ib = solP.resolutionIterativePureConductor(f);

         //*/
        Matrix res = new Matrix(2, solP.getNbLignes());
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }
        RealFaceQuantity[] J = PEEC_FULLY_SURF.exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/SERPENT/f_" + f, fd, res, solP.getDielecCell(), solP.getQ(), true);

        Pertes p = new Pertes(J[0].getElementSet());
        double loss = p.calcul(J[0], J[1], 1 / 1.68e-8, 1.0, 5);
        System.out.println("\nPertes = " + loss + "\n");
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
