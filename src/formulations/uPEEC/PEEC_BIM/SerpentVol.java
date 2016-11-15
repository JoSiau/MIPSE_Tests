/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.PEEC_BIM;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.Arrays;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;

/**
 *
 * @author jsiau
 */
public class SerpentVol {

    /**
     * Serpentin sans bornes.
     *
     * @param args
     */
    public static void main(String[] args) {

        System.out.println("Entrez une frequence:");
        Scanner sc = new Scanner(System.in);
        double f = sc.nextDouble();

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Zf_Serpent.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/SNAKE_CVDS.DEC");//SNAKE_CVDS_WB
        /*
         *******************************
         ***  Import du fichier Flux ***
         *******************************
         Nombre de regions importees : 3
         Region 0, Nom : CU, type : 2, 1400 elements
         Region 1, Nom : FR4, type : 2, 358 elements
         Region 2, Nom : NULL, type : 1, 716 elements
         *******************************
         ***       Fin Import        ***
         *******************************
         */
        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((VolumeRegion) mesh.getRegion(0), 1 / 1.68e-8, 35e-6,
                /*
                 (SurfaceRegion) mesh.getRegion(3), (SurfaceRegion) mesh.getRegion(4),
                 (SurfaceRegion) mesh.getRegion(2), 4.7 * eps0,
                 (SurfaceRegion) mesh.getRegion(1));
                 /*/
                null, null,
                (SurfaceRegion) mesh.getRegion(2), eps0,
                (SurfaceRegion) mesh.getRegion(1));
        //*/
        solP.setPtsDeGaussInductifs(27, 8, 8);
//        solP.setPtsDeGaussInductifs(125, 64, 64);
//        solP.setPtsDeGaussInductifs(64, 27, 27);
        solP.setPtsDeGaussCapacitifs(4, 4);
//        solP.setAnalyticalIntegrationCapa(true);
//        solP.setCorrectionVoisin(true);

//        solP.setEquiPot(new int[]{2468, 1511});
        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int indNoeudBord = fd.getNbElement();
        int nbBranches = solP.getNbLignes();

//        circuitPur.addSourceUSimple(solP.getNbLignes(), indNoeudBord, indNoeudBord + 1, "Source", 240.0, 0.0);
//        circuitPur.addSourceUSimple(nbBranches, 0, 1, "Source", 1.0, 0.0);
//        circuitPur.addSourceUSimple(nbBranches, 4916, 3168, "Source", 1.0, 0.0);
        circuitPur.addSourceUSimple(nbBranches, 1399, 703, "Source", 1.0, 0.0);
//        circuitPur.addSourceUSimple(nbBranches+1, 704, 72, "cct", 0.0, 0.0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

//        solP.setExportTopologie("D:");
//        solP.getTopologie();
//        solP.integration();
//        solP.checkMatricesSingularities();
//        solP.getZ(10);
        //*
//        double[][] ib = solP.resolutionDirecte(f);
        double[][] ib = solP.resolutionDirectePureConductor(f);
        /*/
         double[][] ib = solP.resolutionIterative(f);
        
         //*/

        double I[] = new double[]{ib[0][2 * nbBranches], ib[0][2 * nbBranches + 1]};
        double U[] = new double[]{ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]};
        System.out.println("I= " + Arrays.toString(I));
        System.out.println("U= " + Arrays.toString(U));
        System.out.println("|Z|= " + (Math.hypot(U[0], U[1]) / Math.hypot(I[0], I[1])));

        Matrix res = new Matrix(2, solP.getNbLignes());
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }
        System.out.println("res:\n " + res.toString());
        PEEC_FULLY_SURF.exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/SERPENT/f_" + f, fd, res, solP.getDielecCell(), solP.getQ(), true);

        double b[] = new double[2 * solP.getNbLignes()];
        for (int i = 0; i < b.length; i++) {
            b[i] = Math.random() * 100;
        }
        double pf[] = solP.produit(b, new double[2 * solP.getNbLignes()], f);

        Basic2D Mf = new Basic2D(solP.getZbFull(new double[solP.getNbLignes()][2 * solP.getNbLignes()], 0, solP.getNbLignes() - 1, 0, solP.getNbLignes() - 1, f));

        double rf[] = Mf.product(b, new double[2 * solP.getNbLignes()]);

        ColumnVector vf = new ColumnVector(rf);
        ColumnVector vp = new ColumnVector(pf);

        vp.sub(vf);
        System.out.println("Error rel pmv's= " + vp.norm() / vf.norm());
        System.out.println("Error abs pmv's= " + vp.norm());
        
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    /**
     * Serpentin avec des bornes sur 2 faces.
     *
     * @param args the command line arguments
     */
    public static void mainsurf(String[] args) {

        System.out.println("Entrez une frequence:");
        Scanner sc = new Scanner(System.in);
        double f = sc.nextDouble();

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_SURF.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Surf/S4_FS.DEC");
        /*
         *******************************
         ***  Import du fichier Flux ***
         *******************************
         Nombre de regions importees : 5
         Region 0, Nom : CU, type : 2, 1400 elements
         Region 1, Nom : FR4, type : 2, 358 elements
         Region 2, Nom : POS, type : 1, 4 elements
         Region 3, Nom : NEG, type : 1, 4 elements
         Region 4, Nom : NULL, type : 1, 708 elements
         *******************************
         ***       Fin Import        ***
         *******************************
         */
        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF((SurfaceRegion) mesh.getRegion(0), 1 / 1.68e-8, 35e-6, null, null,//(LineRegion) mesh.getRegion(2), (LineRegion) mesh.getRegion(3),
                (SurfaceRegion) mesh.getRegion(1), 4.7* eps0,
                (SurfaceRegion) mesh.getRegion(0));

//        solP.getTopologie();
        solP.setPtsDeGaussInductifs(64, 25, 25);
        solP.setPtsDeGaussCapacitifs(25, 25);
//        solP.setAnalyticalIntegrationCapa(true);

//        solP.setEquiPot(new int[]{512, 1212});
        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
//        circuitPur.addSourceISimple(400, 125, 130, "t", 1.0, 0.0);
        int indNoeudBord = fd.getNbElement();
        int nbBranches = solP.getNbLignes();
//        circuitPur.addSourceUSimple(nbBranches, indNoeudBord, indNoeudBord + 1, "Source", 240.0, 0.0);
        circuitPur.addSourceUSimple(nbBranches, 699, 1399, "Source", 1.0, 0.0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

//        solP.setExportTopologie("D:");
//        solP.getTopologie();
//        solP.integration();
//        solP.checkMatricesSingularities();
//        solP.getZ(10);
//        double[][] ib = solP.resolutionDirectePureConductor(f);
        double[][] ib = solP.resolutionDirecte(f);
        double I[] = new double[]{ib[0][ib[0].length - 2], ib[0][ib[0].length - 1]};
        double U[] = new double[]{ib[1][ib[1].length - 2], ib[1][ib[1].length - 1]};
        System.out.println("I= " + Arrays.toString(I));
        System.out.println("U= " + Arrays.toString(U));
        System.out.println("|Z|= " + (Math.hypot(U[0], U[1]) / Math.hypot(I[0], I[1])));

        Matrix res = new Matrix(2, fd.getActiveDofCount());
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i] / 35e-6);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1] / 35e-6);
        }

        PEEC_FULLY_SURF.exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/SERPENT/f_" + f, fd, res, solP.getDielecCell(), solP.getQ(), true);

        GestionnaireTaches.getGestionnaireTaches().stop();
    }
}