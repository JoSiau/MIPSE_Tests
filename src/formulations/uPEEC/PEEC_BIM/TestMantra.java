/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.PEEC_BIM;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMesh;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
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
public class TestMantra {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/mantaGEO.msh";
        System.out.println("Fichier maillage : " + path);
        ImportGmshMesh mesh = new ImportGmshMesh(path);
        ((ElementSurfSetHomogene) mesh.getRegion(0).getElementSet()).computeVois();
        ((ElementSurfSetHomogene) mesh.getRegion(0).getElementSet()).resetSurfOrientation();

        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((SurfaceRegion) mesh.getRegion(0), 1 / 1.68e-8, 35e-6,
                null, null
        );

        solP.setPtsDeGaussInductifs(25, 4, 4);
        solP.setPtsDeGaussCapacitifs(4, 4);
//        solP.setAnalyticalIntegrationCapa(true);
        solP.setCorrectionVoisin(false);

        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int indNoeudBord = fd.getNbElement();
        int nbBranches = solP.getNbLignes();
        if (path.endsWith("mantaGEO.msh")) {
            circuitPur.addSourceISimple(nbBranches++, 3905, 3909, "Source", 1.0, 0.0);
            circuitPur.addSourceISimple(nbBranches++, 1710, 3909, "Source", 1.0, 0.0);
            circuitPur.addSourceISimple(nbBranches++, 201, 3909, "Source", 1.0, 0.0);
            circuitPur.addSourceISimple(nbBranches++, 3900, 3909, "Source", 1.0, 0.0);
        } else if (path.endsWith("manta.msh")) {
            circuitPur.addSourceISimple(nbBranches++, 6216, 6219, "Source", 1.0, 0.0);
            circuitPur.addSourceISimple(nbBranches++, 6214, 6219, "Source", 1.0, 0.0);
            circuitPur.addSourceISimple(nbBranches++, 2743, 6219, "Source", 1.0, 0.0);
            circuitPur.addSourceISimple(nbBranches++, 376, 6219, "Source", 1.0, 0.0);
        }
        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

        solP.setExportTopologie("D:");
//        solP.integration();
//        solP.checkMatricesSingularities();
//        solP.getZ(10);
        solP.setHmatrixCompression(false);
        double f = 100;
//        double[][] ib = solP.resolutionIterativePureConductor(f);
        double[][] ib = solP.resolutionDirectePureConductor(f);

        System.err.println("Memory to store the matrices = " + (solP.getMemoryUsed() / 1024) + " ko");

//        double[][] ib = solP.resolutionIterativePureConductor(f);
//        double I[] = new double[]{ib[0][2 * nbBranches], ib[0][2 * nbBranches + 1]};
//        double U[] = new double[]{ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]};
//        System.out.println("I= " + Arrays.toString(I));
//        System.out.println("U= " + Arrays.toString(U));
//        System.out.println("|Z|= " + (Math.hypot(U[0], U[1]) / Math.hypot(I[0], I[1])));
        Matrix res = new Matrix(2, solP.getNbLignes());
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i] / 35e-6);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1] / 35e-6);
        }

        PEEC_FULLY_SURF.exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/Mantra/f_" + f, fd, res, solP.getDielecCell(), solP.getQ(), true);

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
