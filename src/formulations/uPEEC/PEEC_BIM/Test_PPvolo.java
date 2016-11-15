/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.PEEC_BIM;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class Test_PPvolo {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Test_PPvolo.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/PPVOL100X100X1_800E.DEC";
        ImportFlux IF = new ImportFlux(meshDir);
        //
        // Gather the meshes
        //
        // Conductor
        VolumeRegion conductor = (VolumeRegion) IF.getRegion(0);
        RegionsSet regSet = new RegionsSet(new Region[]{conductor});
        // Put the  2 dieledctric regions in an array
//        VolumeRegion dielectrics[] = new VolumeRegion[]{(VolumeRegion) mesh.getRegion(1)};
        //
        SurfaceRegion fluxPos = (SurfaceRegion) IF.getRegion(1);
        SurfaceRegion fluxNeg = (SurfaceRegion) IF.getRegion(2);

        double conduc = 55e6;

        double f = 100;

        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF(conductor, conduc, 1.0,
                fluxPos, fluxNeg);
        //*/
        solP.setPtsDeGaussInductifs(27, 8, 8);
//        solP.setPtsDeGaussInductifs(125, 64, 64);
//        solP.setPtsDeGaussInductifs(64, 27, 27);
        solP.setPtsDeGaussCapacitifs(25, 25);
        solP.setAnalyticalIntegrationCapa(true);
//        solP.setCorrectionVoisin(true);

//        solP.setEquiPot(new int[]{2468, 1511});
        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int indNoeudBord = fd.getNbElement();
        int nbBranches = solP.getNbLignes();

        circuitPur.addSourceISimple(solP.getNbLignes(), indNoeudBord, indNoeudBord + 1, "Source", 1.0, 0.0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

        /*
         double[][] ib = solP.resolutionDirectePureConductor(f);
         /*/
        double[][] ib = solP.resolutionIterativePureConductor(f);

        //*/
        double I[] = new double[]{ib[0][2 * nbBranches], ib[0][2 * nbBranches + 1]};
        double U[] = new double[]{ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]};
        System.out.println("I= " + Arrays.toString(I));
        System.out.println("U= " + Arrays.toString(U));
        double modZ = Math.hypot(U[0], U[1]) / Math.hypot(I[0], I[1]);
        System.out.println("|Z|= " + modZ);
        double omega = 2 * Math.PI * f;
        double cCom = (1.0 / (omega * modZ));
        System.out.println("C_comp = 1/(w * |Z|) = " + cCom);
        double capa = eps0 * 1e-2 / 2e-3;
        System.out.println("C_ana = eps0 * S / e = " + capa);
        System.out.println("Erreur relative = " + (Math.abs(capa - cCom) / capa));

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

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
