/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.condSurf_dielecVol;

import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.contraints.ExternalDriveFaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceConstraint;
import g2elab.mipse.meshCore.elements.ElementVolSetHomogene;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.VolumeRegion;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class tmp {

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
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/2D.DEC");

        RegionsSet regSet = new RegionsSet(mesh);
        VolumeRegion d = (VolumeRegion) regSet.generateUnion(new int[]{0, 1}, "union");
        regSet.addRegion(d);

        FaceDeg1 FD = new FaceDeg1((ElementVolSetHomogene) d.getElementSet(), new FaceConstraint[]{new ExternalDriveFaceConstraint(mesh.getRegion(0)), new ExternalDriveFaceConstraint(mesh.getRegion(1))});
        System.out.println(FD.getActiveDofCount());
        System.out.println(FD.getImplicitConstraintCount()[0] + " , " + FD.getImplicitConstraintCount()[1]);
        System.out.println("FD= " + FD.getElementSet().getElements()[0]);
        
        FaceDeg1 f = new FaceDeg1((ElementVolSetHomogene) mesh.getRegion(0).getElementSet(), new ExternalDriveFaceConstraint(mesh.getRegion(2)));
        System.out.println(f.getActiveDofCount());
        System.out.println(f.getImplicitConstraintAdresses()[0]);
        System.out.println(f.getImplicitConstraintCount()[0]);
    }

}
