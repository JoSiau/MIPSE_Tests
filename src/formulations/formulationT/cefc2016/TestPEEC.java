/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT.cefc2016;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import static formulations.formulationT.cefc2016.Test_TIm.conduc;
import static formulations.formulationT.cefc2016.Test_TIm.ep;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import formulations.uPEEC.PEEC_BIM.Zf_Serpent;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class TestPEEC {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        double f = Test_TIm.f;

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Zf_Serpent.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir = meshDir + Test_TIm.mesh;
        ImportFlux mesh = new ImportFlux(meshDir);//SNAKE_CVDS_WB  LOOPANTENNA_SURF3

        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((SurfaceRegion) mesh.getRegion(0), Test_TIm.conduc, Test_TIm.ep,
                null, null);
        //*/
        solP.setComputeCapa(false);
        solP.setPtsDeGaussInductifs(16, 9, 9);
//        solP.setPtsDeGaussCapacitifs(4, 4);
//        solP.setAnalyticalIntegrationCapa(true);
//        solP.setCorrectionVoisin(true);
        solP.setExportTopologie("d:/tmp");
        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        circuitPur.addSourceISimple(1000000, 960, 1679, "source", 1.0, 0.0);
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
        double loss = p.calcul(J[0], J[1], conduc, ep, 9);
        System.out.println("\nPertes = " + loss + "\n");
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
