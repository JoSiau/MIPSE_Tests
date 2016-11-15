/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.PEEC_BIM;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class tmpTestPlaque {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        double f = 1e6;

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Zf_Serpent.class.getName()).log(Level.SEVERE, null, ex);
        }
//        meshDir = meshDir + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Volumic/WIREON.DEC";
        meshDir = meshDir + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/PLAQUE_VOL.DEC";
        ImportFlux mesh = new ImportFlux(meshDir);

        double ep = 1.0;
        double sigma = 1 / 1.68e-8;
        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((VolumeRegion) mesh.getRegion(0), sigma, ep,
                null, null);
        //*/
        solP.setPtsDeGaussInductifs(27, 8, 8);
        solP.setPtsDeGaussCapacitifs(4, 4);
//        solP.setAnalyticalIntegrationCapa(true);
//        solP.setCorrectionVoisin(true);

//        solP.setEquiPot(new int[]{512, 1212});
        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int indNoeudBord = fd.getNbElement();
        int nbBranches = solP.getNbLignes();

        circuitPur.addSourceISimple(solP.getNbLignes(), 104, 339, "Source", 1.0, 0.0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

//        solP.setExportTopologie("D:");
//        solP.integration();
//        solP.checkMatricesSingularities();
//        solP.getZ(10);
        //*
        //        double[][] ib = solP.resolutionDirecte(f);
        double[][] ib = solP.resolutionDirectePureConductor(f);
        /*/
         solP.setHmatrixCompression(false);
         double[][] ib = solP.resolutionIterativePureConductor(f);

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
        RealFaceQuantity[] J = PEEC_FULLY_SURF.exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/SERPENT/f_" + f, fd, res, solP.getDielecCell(), solP.getQ(), true);

        RealFaceQuantity Jreal = J[0];
        RealFaceQuantity Jimag = J[1];

        Pertes p = new Pertes(Jreal.getElementSet());
        double loss = p.calcul(Jreal, Jimag, sigma, ep, 8);
        System.out.println("\nPertes = " + loss + "\n");

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
