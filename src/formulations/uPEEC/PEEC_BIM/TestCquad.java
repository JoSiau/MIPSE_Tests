/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.PEEC_BIM;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;
import java.io.File;

import java.util.Arrays;

/**
 *
 * @author jsiau
 */
public class TestCquad {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + "\\src\\g2elab\\mipse\\formulationInProgress\\magnetodynamic\\FormulationT\\C_QUAD.DEC";

        ImportFlux mesh = new ImportFlux(path);//SNAKE_CVDS_WB  LOOPANTENNA_SURF3

        double ep = 1e-4;
        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((SurfaceRegion) mesh.getRegion(0), 55e6, ep,
                null, null);

        solP.setPtsDeGaussInductifs(9, 4, 4);
        solP.setPtsDeGaussCapacitifs(4, 4);
        solP.setAnalyticalIntegrationCapa(true);
//        solP.setCorrectionVoisin(true);

        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int indNoeudBord = fd.getNbElement();
        int nbBranches = solP.getNbLignes();

        circuitPur.addSourceISimple(solP.getNbLignes(), 0, 61, "src", 1.0, 0.0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

//        solP.setExportTopologie("D:");
//        solP.integration();
//        solP.checkMatricesSingularities();
//        solP.getZ(10);
        /*
         //        double[][] ib = solP.resolutionDirecte(f);
         double[][] ib = solP.resolutionDirectePureConductor(f);
         /*/
        solP.setHmatrixCompression(false);
        double f = 5000;
        double[][] ib = solP.resolutionDirectePureConductor(f);

        //*/
        double I[] = new double[]{ib[0][2 * nbBranches], ib[0][2 * nbBranches + 1]};
        double U[] = new double[]{ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]};
        System.out.println("I= " + Arrays.toString(I));
        System.out.println("U= " + Arrays.toString(U));
        System.out.println("|Z|= " + (Math.hypot(U[0], U[1]) / Math.hypot(I[0], I[1])));

        Matrix res = new Matrix(2, solP.getNbLignes());
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i] / ep);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1] / ep);
        }
        System.out.println("res:\n " + res.toString());
        PEEC_FULLY_SURF.exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/SERPENT/f_" + f, fd, res, solP.getDielecCell(), solP.getQ(), true);

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
