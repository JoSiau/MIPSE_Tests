/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.condSurf_dielecVol;

import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.*;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.LineRegion;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.meshCore.tree.faceTree.Face;

import java.io.IOException;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_Surf_Dielec_Vol.PEEC_RLMPC_SURF.eps0;

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
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/COND_IN_AIR.DEC");//SNAKE_DIELEC_SURF_4

        double ep = 35e-6;
        double sigma = 1 / 2.836E-8;
        double epsilon = 3.4 * eps0;
        double f = 1e9;

//        ////////////////////////////////////       Test L dielec/dielec
//        FaceDeg1 FSd = new FaceDeg1((ElementSetHomogene) mesh.getRegion(0).getElementSet(), new ExternalDriveFaceConstraint(mesh.getRegion(1)));
//        RowVector q = new RowVector(FSd.getNbElement());
//        double omega = 2 * Math.PI * f;
//        double cst = mu0 * omega / (4 * Math.PI);
//        q.setAllElements(cst * (epsilon - eps0) / epsilon);
//        Cell support = new Cell(FSd.getElementSet());
//        RealScalarCellQuantity quant = new RealScalarCellQuantity(support, q);
//        GalerkinIntegralFormulationFull IF = new GalerkinIntegralFormulationFull(FSd, FSd, new MultGvect(), new NoCorrection(new FixedGaussSource(8)), 27);
//        IF.assembly(quant);
//        Matrix M = ((StorageFull) IF.getStore()).getMatrix();
//        System.out.println("M=\n"+M);
//        ////////////////////////////////////       Test R
//        Region Conductor = mesh.getRegion(1);
//        Region Dielectric = mesh.getRegion(0);
//
//        FaceDeg1 FSc = new FaceDeg1((ElementSetHomogene) Conductor.getElementSet(), new FaceConstraint[]{new FaceRealDirichlet(mesh.getRegion(4), 0.0),
//            new ExternalDriveFaceConstraint(mesh.getRegion(2)), new ExternalDriveFaceConstraint(mesh.getRegion(3))});
//
//        // Calcul du bord
//        SurfaceRegion DielecBorder = (SurfaceRegion) mesh.getRegion(0).generateBorder();
//        // Enleve le conducteur
//        ElementSurfSetHomogene DBwoCond = (ElementSurfSetHomogene) DielecBorder.getComplement(new SurfaceRegion(new ElementSurfSetHomogene(duplicate((ElementSetHomogene) Conductor.getElementSet()))), 2);
//        SurfaceRegion gamma = new SurfaceRegion(DBwoCond);
//
//        FaceDeg1 FSdielec = new FaceDeg1((ElementVolSetHomogene) Dielectric.getElementSet(), new FaceConstraint[]{new ExternalDriveFaceConstraint(gamma), new ExternalDriveFaceConstraint(Conductor)});
//
//        Cell GammaExt = new Cell((ElementSurfSetHomogene) gamma.generateUnionWith(Conductor).getElementSet());
//
//
//        RowVector q = new RowVector(FSc.getNbElement());
//        q.setAllElements(2.836e-8 / 35e-6);
//        Cell support = new Cell(FSc.getElementSet());
//        RealScalarCellQuantity iSigma = new RealScalarCellQuantity(support, q);
//        FiniteElementFormulation EF = new FiniteElementFormulation(FSc);
//        EF.assembly(4, iSigma);
//        System.out.println("EF=\n" + EF.getStore().getMatrixPrecond(null).getFullMatrix().toString());
//
//        return;
//            this.FScond = new FaceDeg1((ElementSurfSetHomogene) Conductor.getElementSet(), new FaceConstraint[]{new FaceRealDirichlet(borderCond, 0.0),
//                new ExternalDriveFaceConstraint(Borne1), new ExternalDriveFaceConstraint(Borne2)});
        PEEC_RLMPC_SURF solP = new PEEC_RLMPC_SURF((SurfaceRegion) mesh.getRegion(1), (LineRegion) mesh.getRegion(4), sigma, ep, (VolumeRegion) mesh.getRegion(0), epsilon, (LineRegion) mesh.getRegion(2), (LineRegion) mesh.getRegion(3),
                3, false,  false);

        FaceDeg1 FSd = solP.getFD1_Dielectric();
        Cell gammaExt = solP.getGammaExt();
        FaceDeg1 FSc = solP.getFD1_Conductor();

        //
        // Verifie la numerotation des faces de bords avec celle de la renumerotation du dielectric
        //
        System.out.println("=================================================");
        System.out.println("============  CHECK DIELEC EXTERIOR  ============");
        System.out.println("=================================================");
        System.out.println("NbDielecFace= " + FSd.getImplicitConstraintCount()[0]);
//        System.out.println("NbConductorFace= " + FSd.getImplicitConstraintCount()[1]);
        System.out.println("NbElmentExt= " + gammaExt.getActiveDofCount());
        for (int i = 0; i < gammaExt.getActiveDofCount() - solP.getNbElementInAir(); i++) {
            Face fd = FSd.getFacesSet().getListFaces()[FSd.getImplicitConstraintAdresses()[0] + i];
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
        System.out.println("============  CHECK DIELEC CONDUCTOR  ============");
        System.out.println("==================================================");
        for (int i = 0; i < FSc.getNbElement() - solP.getNbElementInAir(); i++) {
            Face fd = FSd.getFacesSet().getListFaces()[FSd.getImplicitConstraintAdresses()[1] + i];
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

        if (solP.getNbElementInAir() != 0) {

            System.out.println("====================================================");
            System.out.println("============  CHECK CONDUCTOR EXTERIOR  ============");
            System.out.println("====================================================");

            int offSetE = gammaExt.getNbElement() - solP.getNbElementInAir();
            int offSetC = FSc.getNbElement() - solP.getNbElementInAir();
            for (int i = 0; i < solP.getNbElementInAir(); i++) {
                Face fd = new Face(FSc.getElementSet().getElements(offSetC + i));
                Face fe = new Face(gammaExt.getElementSet().getElements(i + offSetE));

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

        }

        System.out.println("============ Fin du check ============");

        System.out.println("Tapez ce que vous voulez pour continuer ! Mais pas les femmes !");
        new Scanner(System.in).next();

        int matTopo[][] = solP.getTopologie();

        System.out.println("FSc.ActiveDof= " + FSc.getActiveDofCount());
        Cell C = new Cell(FSc.getElementSet());

        int lienFaceElement[][] = FSc.getFacesSet().getPosElementFace();

        for (int i = 0; i < lienFaceElement.length; i++) {
            for (int j = 0; j < lienFaceElement[i].length; j++) {
                System.out.print(lienFaceElement[i][j] + " , ");
            }
            System.out.println("");
        }

        int ind[] = solP.getIndexNodeBorne1();
        System.out.print("indB1= ");
        for (int i = 0; i < ind.length; i++) {
            System.out.print(ind[i] + " ,");
        }
        System.out.println("");

        ind = solP.getIndexNodeBorne2();
        System.out.print("indB2= ");
        for (int i = 0; i < ind.length; i++) {
            System.out.print(ind[i] + " ,");
        }
        System.out.println("");

        Face Fd[] = FSd.getFacesSet().getListFaces();
        System.out.println("Fd.length= " + Fd.length + "\t C.dof= " + C.getActiveDofCount());
        int indF = 35;
        int indT = 3;
        System.out.println("T0= \n" + FSc.getElementSet().getElements(indT).getNoeuds()[0] + "\n" + FSc.getElementSet().getElements(indT).getNoeuds()[1] + "\n" + FSc.getElementSet().getElements(indT).getNoeuds()[2] + "\n" + FSc.getElementSet().getElements(indT).getNoeuds()[3] + "\t");
        System.out.println("Find= \n" + Fd[indF].getNodes()[0] + "\n" + Fd[indF].getNodes()[1] + "\n" + Fd[indF].getNodes()[2] + "\n" + Fd[indF].getNodes()[3] + "\n");

    }
}
