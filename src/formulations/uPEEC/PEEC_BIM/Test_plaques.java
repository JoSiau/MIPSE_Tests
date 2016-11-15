/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.PEEC_BIM;

import g2elab.mipse.analytical.sourceFields.CircularCoil;
import g2elab.mipse.analytical.sourceFields.RealCoilField;
import g2elab.mipse.analytical.sourceFields.SourceField;
import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF.exportRes;

/**
 *
 * @author jsiau
 */
public class Test_plaques {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        //        path = path + "/src/g2elab/mipse/formulationInProgress/electrostatic/potential/PLAQUE2000.DEC";
//        path = path + "/src/g2elab/mipse/stabilizedFormulations/formulationsGM/fichiersDEC/PLAQUE_COIL03_QUAD_1600.DEC";
        path = path + "/src/g2elab/mipse/stabilizedFormulations/formulationsGM/fichiersDEC/PLAQUE_COIL03_TRI_1020.DEC";
        ImportFlux mesh = new ImportFlux(path);
        double ep = 1e-4;
        double f = 100;
        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((SurfaceRegion) mesh.getRegion(0), 55e6, ep, null, null);
        FaceDeg1 fd = solP.getFD1();
        //
        solP.setPtsDeGaussInductifs(16, 9, 9);
        solP.setPtsDeGaussCapacitifs(4, 4);
        //
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);
        //
        /*/
        SourceField[] srcField = new SourceField[]{new UniformField(0.0, 0.0, 1.0)};
        /*/
        // Définition d'une bobine circulaire 
        double zcentre = 0.005;
        double[] c = {0.0, 0.0, zcentre};
        double[] u = {1.0, 0.0, 0.0};
        double[] v = {0.0, 1.0, 0.0};
        double[] n = {0.0, 0.0, 1.0};
        double rayon = 0.005;
        double epaisseurSpire = rayon / 20;
        double courant = 1.0;

        ColumnVector centre = new ColumnVector(c);
        ColumnVector uu = new ColumnVector(u);
        ColumnVector vv = new ColumnVector(v);
        ColumnVector nn = new ColumnVector(n);

        CircularCoil circularCoil = new CircularCoil(rayon, epaisseurSpire, centre, uu, vv, nn, mesh.getRegion(0).getElementSet().getElements()[0].getNoeuds()[0].getElementFactory());
        SourceField srcField[] = new SourceField[]{new RealCoilField(circularCoil, courant)};
        //*/
        solP.setSourceField(srcField);
        //
//        solP.setHmatrixCompression(true);
        double[][] ib = solP.resolutionIterativePureConductor(f);
        //
        Matrix current = new Matrix(2, fd.getActiveDofCount());
        System.out.println("ep = "+ep);
        for (int i = 0; i < current.getColumnCount(); i++) {
            current.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            current.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }
        System.out.println("res = \n" + current);

        Pertes p = new Pertes(fd);
        double loss = p.calcul(new RealFaceQuantity(fd, current.row(0)), new RealFaceQuantity(fd, current.row(1)), 55e6, ep, 9);

        System.out.println("\nPertes = " + loss + "\n");

        exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/TestPlaque", fd, current, null, null, true);

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

//    public static void main1(String[] args) {
//        File file = new File("");
//        String path;
//        path = file.getAbsolutePath();
//        //        path = path + "/src/g2elab/mipse/formulationInProgress/electrostatic/potential/PLAQUE2000.DEC";
//        path = path + "/src/g2elab/mipse/stabilizedFormulations/formulationForGUI/plaques/PLAQUE_625.DEC";// 625, 5625, 2500, 10000, 40000,90000,160000
//        ImportFlux mesh = new ImportFlux(path);
//        double ep = 1e-4;
//        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((SurfaceRegion) mesh.getRegion(0), 55e6, ep, null, null,
//                (SurfaceRegion) mesh.getRegion(0), 4.7 * eps0,
//                (SurfaceRegion) mesh.getRegion(0));
//        FaceDeg1 fd = solP.getFD1();
//        //
//        solP.setPtsDeGaussInductifs(9, 4, 4);
//        solP.setPtsDeGaussCapacitifs(4, 4);
//        //
//        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
//        circuitPur.finSaisie();
//        solP.setCircuitElectrique(circuitPur);
//        //
//        solP.setSourceField(new SourceField[]{new UniformField(0.0, 0.0, 1.0)});
//        //
//        solP.setHmatrixCompression(true);
//        double[][] ib = solP.resolutionIterativePureConductor(1e2);
//        //
//        Matrix current = new Matrix(2, fd.getActiveDofCount());
//        for (int i = 0; i < current.getColumnCount(); i++) {
//            current.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i] / ep);
//            current.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1] / ep);
//        }
//        System.out.println("res = \n" + current);
//
//        Pertes p = new Pertes(fd);
//        double loss = p.calcul(new RealFaceQuantity(fd, current.row(0)), new RealFaceQuantity(fd, current.row(1)), 55e6, 9);
//
//        System.out.println("\nPertes = " + loss + "\n");
//
//        exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/TestPlaque", fd, current, null, null, true);
//
//        GestionnaireTaches.getGestionnaireTaches().stop();
//    }

}
