/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.PEEC_BIM;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMesh;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.numericalTools.matrix.MatriceIncidence;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;
import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.PEEC_FULLY_SURF.exportRes;

/**
 *
 * @author jsiau
 */
public class TestBobZach {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(TestBobZach.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportGmshMesh mesh = new ImportGmshMesh(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/bob2.msh");
        mesh.meshSummary();
        //
        // Gather the meshes
        //
        // Conductor
        VolumeRegion conductor = (VolumeRegion) mesh.getRegion(6);
        // Put the  2 dieledctric regions in an array
        VolumeRegion dielectrics[] = new VolumeRegion[]{(VolumeRegion) mesh.getRegion(7), (VolumeRegion) mesh.getRegion(8)};
        // Interface between the conductor and the dielectrics
        SurfaceRegion interfaceCondDielec = (SurfaceRegion) mesh.getRegion(1);
        // Get the dielectrics border
        SurfaceRegion bobAir = (SurfaceRegion) mesh.getRegion(2);
        SurfaceRegion formerAir = (SurfaceRegion) mesh.getRegion(3);
        // And compute their union
        SurfaceRegion dielecAir = (SurfaceRegion) bobAir.generateUnionWith(formerAir);
        dielecAir.exportGmsh("d:/dielecAir.msh");
        SurfaceRegion fluxPos = (SurfaceRegion) mesh.getRegion(4);
        SurfaceRegion fluxNeg = (SurfaceRegion) mesh.getRegion(5);

        double ep = 1.0;
        PEEC_FULLY_SURF solP = new PEEC_FULLY_SURF(conductor, 5.8502e7, ep, fluxPos, fluxNeg, bobAir, 3.2 * eps0, interfaceCondDielec);

        double f = 8376776.401;//8e6;

        solP.setPtsDeGaussInductifs(27, 8, 8);
//        solP.setPtsDeGaussInductifsVol(125, 64, 64);
        solP.setPtsDeGaussCapacitifs(4, 4);
//        solP.setFullAnalyticalP(true);
        FaceDeg1 fd = solP.getFD1();
        //
        //
        //
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        int indNoeudBord = fd.getNbElement();
        circuitPur.addSourceUSimple(solP.getNbLignes(), indNoeudBord, indNoeudBord + 1, "Source", 1.0, 0.0);

        circuitPur.finSaisie();

        solP.setCircuitElectrique(circuitPur);

        MatriceIncidence mi = solP.getSolveurCircuit().getMI().transpose();
        mi.saveMTX("d:/tmp/Mi.mtx");

        solP.setExportTopologie("D:");
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // Resolution
        /*
         boolean HmatComp = true;
         solP.setCompression(HmatComp?Compression.HCA:Compression.No);// Enleve la com pression
         solP.setParamIterativeSolver(1, 10000, -1e-8, 1, -150);
         solP.setParamPreconditionner(1, 500, 0, new double[]{500, -3e-1, 500});
         double ib[][] = solP.resolutionIterative(f);
         /*/
        double ib[][] = solP.resolutionDirecte(f);
        //*/
        Matrix res = new Matrix(2, fd.getActiveDofCount());
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i] / ep);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1] / ep);
        }

        exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/Test", fd, res, solP.getDielecCell(), solP.getQ(), true);

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
