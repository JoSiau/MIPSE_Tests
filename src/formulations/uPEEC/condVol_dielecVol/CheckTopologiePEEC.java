/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.condVol_dielecVol;

import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.Node;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.meshCore.tree.faceTree.Face;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv.eps0;

/**
 *
 * @author jsiau
 */
public class CheckTopologiePEEC {

    public static void main(String[] args) {
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(CheckTopologiePEEC.class.getName()).log(Level.SEVERE, null, ex);
        }

        String file = meshDir + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/PP_D16.DEC";
        ImportFlux Imp = new ImportFlux(file);
        int nDof[] = new int[Imp.getRegions().length];
        for (int i = 0; i < nDof.length; i++) {
            nDof[i] = Imp.getRegion(i).getElementSet().getNbElement();
        }
        VolumeRegion rV[] = new VolumeRegion[2];
        rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice
        rV[1] = (VolumeRegion) Imp.getRegion(1);// Region dielectrique

        double f = 100;
        /*
         CONSTRUCTION DU SUPPORT DES CONDUCTIVITES EQUIVALENTES
         */
        double Sigma[][] = new double[2][nDof[0] + nDof[1]];// Size = [2][nbElmttotal]
        // Conducteur
        double Re = 5.814e8;
        double Im = eps0;
        for (int i = 0; i < nDof[0]; i++) {
            Sigma[0][i] = Re;// Conductivite du cuivre
            Sigma[1][i] = Im;// Epsilon0
        }
        // Dielectrique
        Re = 0;
        double epsilon = 4.7 * eps0;
        for (int i = nDof[0]; i < nDof[0] + nDof[1]; i++) {
            Sigma[0][i] = Re;// 0 pour le dielec
            Sigma[1][i] = epsilon;// Espilon
        }

        // Create the solver
        PEEC_RLMPC_ConducEquiv solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma, (SurfaceRegion) Imp.getRegion(2), (SurfaceRegion) Imp.getRegion(3), true);
        solP.setPtsDeGaussInductifs(27, 8);
        solP.setPtsDeGaussCapacitifs(4, 4);

        FaceDeg1 FSd = solP.getFunctionSpace();
        Cell gammaExt = solP.getGammaExt();
        //
        // Verifie la numerotation des faces de bords avec celle de la renumerotation du dielectric
        //
        System.out.println("=================================================");
        System.out.println("============  CHECK DIELEC EXTERIOR  ============");
        System.out.println("=================================================");
        for (int i = 0; i < gammaExt.getActiveDofCount(); i++) {
            Face fd = FSd.getFacesSet().getListFaces()[FSd.getImplicitConstraintAdresses()[1] + i];
            Face fe = new Face(gammaExt.getElementSet().getElements(i));

            if (!fd.compareFace(fe)) {
                System.out.println("faces differentes #" + i + ": ");
                System.out.println("fd= \n" + fd);
                Node n[] = fd.getNodes();
                for (int j = 0; j < n.length; j++) {
                    System.out.println("n[" + n[j].getGlobalNum() + "] = {" + n[j].getCoord(0) + " , " + n[j].getCoord(1) + " , " + n[j].getCoord(2));
                }

                System.out.println("fe= \n" + fe);
                n = fe.getNodes();
                for (int j = 0; j < n.length; j++) {
                    System.out.println("n[" + n[j].getGlobalNum() + "] = {" + n[j].getCoord(0) + " , " + n[j].getCoord(1) + " , " + n[j].getCoord(2));
                }

                System.out.println("");
            }
        }

        System.out.println("==================================================");
        System.out.println("============  CHECK DIELEC INTERFACE  ============");
        System.out.println("==================================================");
        Cell FSc = solP.getGammaInt();
        for (int i = 0; i < FSc.getNbElement(); i++) {
            Face fd = FSd.getFacesSet().getListFaces()[FSd.getImplicitConstraintAdresses()[0] + i];
            Face fe = new Face(FSc.getElementSet().getElements(i));

            if (!fd.compareFace(fe)) {
                System.out.println("faces differentes #" + i + ": ");
                System.out.println("fd= \n" + fd);
                Node n[] = fd.getNodes();
                for (int j = 0; j < n.length; j++) {
                    System.out.println("n[" + n[j].getGlobalNum() + "] = {" + n[j].getCoord(0) + " , " + n[j].getCoord(1) + " , " + n[j].getCoord(2));
                }

                System.out.println("fe= \n" + fe);
                n = fe.getNodes();
                for (int j = 0; j < n.length; j++) {
                    System.out.println("n[" + n[j].getGlobalNum() + "] = {" + n[j].getCoord(0) + " , " + n[j].getCoord(1) + " , " + n[j].getCoord(2));
                }

                System.out.println("");
            }
        }
        System.out.println("============ Fin du check ============");

    }

}
