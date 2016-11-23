/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.ArticlesCodesJS.soutenance;

import static formulations.uPEEC.ArticlesCodesJS.soutenance.Comparison_CompTechs.loadMesh;
import static formulations.uPEEC.ArticlesCodesJS.soutenance.Comparison_CompTechs.mesh;
import static formulations.uPEEC.ArticlesCodesJS.soutenance.DebugFMM.getErr;
import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.mipseCore.blockAssembly.real.RealMatrixBlocMult;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.numericalTools.matrix.real.AbstractMatrixReal;
import g2elab.mipse.numericalTools.matrix.real.dense.gotMatrix.GOTMatrix;
import got.matrix.Matrix;

/**
 *
 * @author jsiau
 */
public class DebugFMM_matP {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        AbstractMatrixReal[] pfmm = getSolver(Compression.FMM).integP();
        AbstractMatrixReal[] p = getSolver(Compression.No).integP();

        System.out.println("Erreur Pe = " + getErr(p[0], pfmm[0]));
        System.out.println("Erreur Pi = " + getErr(p[1], pfmm[1]));

        Matrix pi = ((StorageFull) p[1]).getMatrix();
        AbstractMatrixReal pifmm = ((RealMatrixBlocMult) pfmm[1]).getMatrice2();
        System.out.println("Erreur Pi2 = " + getErr(new GOTMatrix(pi), pifmm));
    }

    public static PEEC_RLMPC_Volume getSolver(Compression compression) {

        if (mesh == null) {
            loadMesh();
        }
        ///////////// Gather the meshes
        // Conductor
        VolumeRegion conductor = (VolumeRegion) mesh.getRegion(0);
        // Put the  2 dieledctric regions in an array
        VolumeRegion dielectrics[] = new VolumeRegion[]{(VolumeRegion) mesh.getRegion(1)};
        // Interface between the conductor and the dielectrics
        SurfaceRegion interfaceCondDielec = (SurfaceRegion) mesh.getRegion(2);
        // And compute their union
        SurfaceRegion dielecAir = (SurfaceRegion) mesh.getRegion(3);
        SurfaceRegion fluxPos = null;
        SurfaceRegion fluxNeg = null;
        //
        PEEC_RLMPC_Volume solP = new PEEC_RLMPC_Volume(conductor, 1 / 1.68e-8,
                dielectrics, new double[]{4.7 * eps0},
                fluxPos, fluxNeg, interfaceCondDielec, dielecAir);
        //
        solP.setPtsDeGaussInductifsVol(27, 8, 8);
        solP.setPtsDeGaussCapacitifs(4, 4);
        solP.setFullAnalyticalP(true);
        //
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
        int nbBranches = solP.getNbLignes();
        if (mesh.equals("SERP_3EP")) {
            circuitPur.addSourceISimple(nbBranches, 14351, 14271, "Source I", 1.0, 0.0);
        } else {

        }
        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);
        solP.setExportTopologie("D:");
        solP.setCompression(compression);
        //
        return solP;
    }

}
