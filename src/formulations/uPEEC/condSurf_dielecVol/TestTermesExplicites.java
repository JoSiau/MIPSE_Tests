/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.condSurf_dielecVol;

import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.region.LineRegion;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D;
import g2elab.mipse.numericalTools.vector.sparse.SparseVectorComplex;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;

import java.io.IOException;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF.eps0;

/**
 *
 * @author jsiau
 */
public class TestTermesExplicites {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_SURF.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";

        System.out.println("Choisir le maillage:");
        System.out.println("       1: LoopAntenna");
        System.out.println("       2: Condensateur (tetraedrique)");
        System.out.println("       3: Condensateur (hexaedrique)");
        System.out.println("       4: Poutre (hexaedrique)");
        System.out.println("       5: Plaque (hexaedrique)");
        System.out.println("       6: Serpentin");
        System.out.println("       7: Plaques Paralleles format 'Europe' ");
        System.out.println("       8: Plaques Paralleles 10x10x1 mm");
        System.out.println("       9: Plaques Voisines");
        Scanner sc = new Scanner(System.in);
        String app;
        String doss;
        double sigma = 5.814 * 1e8;
        double epsilon_fr4 = 4.7;
        double epsilon_air = 1;
        double epsilon = epsilon_fr4 * eps0;
        boolean isTetra;
        double ep = 35e-6;
        // Utiliser uniquement pour calculer la capa plan.
        // Stoque les dimensions de la plaque et le gap entre elles.
        double dim[] = null;
        // Conducteur, Null, Dielectric, Borne Pos, Borne Neg.
        int ind[] = new int[]{1, 4, 0, 2, 3};
        int cct[] = null;

        app = "COND_IN_AIR";
        doss = "CondInAir";
        isTetra = false;
        epsilon = 4.7 * eps0;
        sigma = 1 / 1.72E-8;
        // Conducteur, Null, Dielectric, Borne Pos, Borne Neg.
        ind = new int[]{1, 4, 0, 2, 3};

        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/" + app + ".DEC");
        SurfaceRegion conductor = (SurfaceRegion) mesh.getRegion(ind[0]);
        LineRegion bordNull = (LineRegion) mesh.getRegion(ind[1]);
        VolumeRegion Dielec = (VolumeRegion) mesh.getRegion(ind[2]);
        LineRegion borne1 = null, borne2 = null;
        borne1 = (LineRegion) mesh.getRegion(ind[3]);
        borne2 = (LineRegion) mesh.getRegion(ind[4]);

        double f = 300;

        PEEC_RLMPC_SURF solP = new PEEC_RLMPC_SURF(conductor, bordNull, sigma, ep, Dielec, epsilon, borne1, borne2,
                3, false, false);

        solP.setPtsDeGaussInductifsVol(27, 8, 8);
        solP.setPtsDeGaussInductifsSurf(16, 9, 9);
        solP.setPtsDeGaussCapacitifs(4, 4);
//            solP.setFullAnalyticalP(true);

//        solP.setCheckNaN(true);

        boolean HmatCompression = false;
        solP.setCompression(HmatCompression);

        double b[] = new double[2 * solP.getNbLignes()];
        for (int i = 0; i < b.length; i++) {
            b[i] = Math.random() * 100;
        }
        double pf[] = solP.produit(b, new double[2 * solP.getNbLignes()], f);

        if (HmatCompression) {
            solP.setCompression(false);
            solP.integrationMono();
        }
        Basic2D Mf = new Basic2D(solP.getZbFull(new double[solP.getNbLignes()][2 * solP.getNbLignes()], 0, solP.getNbLignes() - 1, 0, solP.getNbLignes() - 1, f));

        double rf[] = Mf.product(b, new double[2 * solP.getNbLignes()]);

        ColumnVector vf = new ColumnVector(rf);
        ColumnVector vp = new ColumnVector(pf);

        vp.sub(vf);
        System.out.println("Error pmv's= " + vp.norm() / vf.norm());

        System.out.println("If nothing is printed after that, it means that the error = 0 !");
        solP.setCompression(true);
        solP.integrationMono();

        double errMax = -1, errMin = 1e32;
        for (int i = 0; i < solP.getNbLignes(); i++) {

            ColumnVector M = new ColumnVector(Mf.getArray()[i]);
            SparseVectorComplex m = solP.getTermesExplicites(i, f);
            double[] vv = m.getVv();
            int[] vi = m.getVi();
            double normM = M.norm();
            for (int j = 0; j < vi.length; j++) {
                M.setElement(2 * vi[j], M.getElement(2 * vi[j]) - vv[2 * j]);
                M.setElement(2 * vi[j] + 1, M.getElement(2 * vi[j] + 1) - vv[2 * j + 1]);
            }
            double Mnorm = M.norm();
            if (Mnorm != 0) {
                System.out.println("Erreur relative " + i + " = " + Mnorm / normM);
            }
            if (Mnorm / normM > errMax) {
                errMax = Mnorm / normM;
            }
            if (Mnorm / normM < errMin) {
                errMin = Mnorm / normM;
            }
        }
        System.out.println("\n\n");
        System.out.println("Erreur relative max = "+errMax);
        System.out.println("Erreur relative min = "+errMin);
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
