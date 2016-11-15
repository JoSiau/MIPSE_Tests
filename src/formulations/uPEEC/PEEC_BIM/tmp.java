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
import g2elab.mipse.numericalTools.matrix.MatriceIncidence;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;
import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF.exportRes;
import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Preconditioners.SPQR.Real.SPQR_c.printSparseMatrixForSPQR;

/**
 *
 * @author jsiau
 */
public class tmp {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(PEEC_FULLY_SURF.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        PEEC_FULLY_SURF solP;
        double ep = 35e-6;
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Surf/LASURF_TMP.DEC");
        //*
        solP = new PEEC_FULLY_SURF((SurfaceRegion) mesh.getRegion(0), 1 / 1.68e-8, ep, null, null,// (LineRegion) mesh.getRegion(2), (LineRegion) mesh.getRegion(3),
                (SurfaceRegion) mesh.getRegion(0), 4.7 * eps0,
                (SurfaceRegion) mesh.getRegion(0));
        solP.setPtsDeGaussInductifs(9, 4, 4);
        solP.setPtsDeGaussCapacitifs(4, 4);
        solP.computeSurfaces();
        solP.integrationConductor();
        Matrix L = solP.getL();
        Matrix P = solP.getP();

        printMatrixForMatLab(L, "L");
        printMatrixForMatLab(P, "P");

        printTopoForMatLab(solP, L.getRowCount());
        

        BlocElectriqueBasique circuitE = new BlocElectriqueBasique();
        circuitE.addSourceISimple(1000, 27, 22, "Src", 1.0, 0.0);
        circuitE.finSaisie();
        solP.setCircuitElectrique(circuitE);

        MatriceIncidence mi = solP.getSolveurCircuit().getMI();
        printSparseMatrixForSPQR(mi, "Mi");
        double[][] ib = solP.resolutionDirectePureConductor(5e9);
        System.out.println("I= " + ib[0][ib[0].length - 2] + " +j* " + ib[0][ib[0].length - 1]);
        System.out.println("U= " + ib[1][ib[0].length - 2] + " +j* " + ib[1][ib[0].length - 1]);
        System.out.println("U= " + ib[1][ib[0].length - 2] + " +j* " + ib[1][ib[0].length - 1]);

        FaceDeg1 fd = solP.getFD1();

        Matrix res = new Matrix(2, fd.getActiveDofCount());
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i] / ep);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1] / ep);
        }

        exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/Test", fd, res, solP.getDielecCell(), solP.getQ(), true);


        GestionnaireTaches.getGestionnaireTaches().stop();
        
        
        
        Double bb[];
        bb = new Double[]{(double) 20, (double) 20};
        Double[] vb = {(double) 20, (double) 20};
    }

    public static void printTopoForMatLab(PEEC_FULLY_SURF solP, int nL) throws IOException {
//        solP.setExportTopologie("d:");
        int[][] topo = solP.getTopologie();
        Ecriture f = new Ecriture("D:/MatLab/topo.out");
        for (int i = 0; i < nL; i++) {
            char type = topo[i][3] == -2 ? 'L' : (topo[i][3] == -3 ? 'C' : 'R');
            f.ecrire(type + " L" + i + " " + topo[i][0] + " N1=" + topo[i][1] + " N2=" + topo[i][2] + " 1\n");
        }
        f.close();
    }

    public static void printMatrixForMatLab(Matrix M, String name) throws IOException {
        Ecriture f = new Ecriture("D:/MatLab/" + name + ".mat");
        double t[] = new double[M.getColumnCount()];
        for (int i = 0; i < M.getRowCount(); i++) {
            M.getRow(i, t);
            f.ecrire(t, ',');
            f.ecrire(";\n");

        }
        f.close();
    }


}
