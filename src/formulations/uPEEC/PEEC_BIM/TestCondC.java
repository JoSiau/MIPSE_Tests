/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.PEEC_BIM;

import g2elab.mipse.analytical.sourceFields.SourceField;
import g2elab.mipse.analytical.sourceFields.UniformMagneticField;
import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;
import java.io.File;

/**
 *
 * @author jsiau
 */
public class TestCondC {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + "\\src\\g2elab\\mipse\\formulationInProgress\\magnetodynamic\\U_PEEC_DIELECTRIC\\C_MED.DEC";
//       
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
        //
        double f = 1000;
        double ep = 1e-4; // 1mm
        double conduc = 55e6;
        ((SurfaceRegion) IF.getRegion(0)).setThickness(ep);

        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((SurfaceRegion) IF.getRegion(0), conduc, ep,
                null, null);
        solP.setPtsDeGaussInductifs(6, 3,3);
        solP.setPtsDeGaussCapacitifs(3, 3);
//        solP.setAnalyticalIntegrationCapa(true);

        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

        solP.setSourceField(new SourceField[]{new UniformMagneticField(0, 0, 1)});
//        solP.setExportTopologie("D:");
//        solP.integration();
//        solP.checkMatricesSingularities();
//        solP.getZ(10);
        /*
        //        double[][] ib = solP.resolutionDirecte(f);
        double[][] ib = solP.resolutionDirectePureConductor(f);
        /*/
         solP.setHmatrixCompression(false);
         double[][] ib = solP.resolutionIterativePureConductor(f);
         //*/
        Matrix res = new Matrix(2, solP.getNbLignes());
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);// / ep);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);// / ep);
        }
        System.out.println("res:\n " + res.toString());
        PEEC_FULLY_SURF.exportRes("D:/jsiau/_Backup_Sources/Resultats/C/f_" + f, fd, res, solP.getDielecCell(), solP.getQ(), true);

        RealFaceQuantity Jreal = new RealFaceQuantity(fd, res.row(0).subvector(0, fd.getActiveDofCount()));
        RealFaceQuantity Jimag = new RealFaceQuantity(fd, res.row(1).subvector(0, fd.getActiveDofCount()));
        Pertes p = new Pertes(fd);
        double pp = p.calcul(Jreal, Jimag, conduc, ep, 9);
        System.out.println("Pertes = " + pp);

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
