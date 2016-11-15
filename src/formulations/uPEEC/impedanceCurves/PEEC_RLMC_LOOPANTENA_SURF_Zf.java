/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.impedanceCurves;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv;
import formulations.uPEEC.PEECwithCapa_Full;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.IO.paraview.ExportVtkHdiv;
import g2elab.mipse.meshCore.contraints.ExternalDriveFaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceRealDirichlet;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class PEEC_RLMC_LOOPANTENA_SURF_Zf {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        double epaisseur = 2 * 1e-5;
        int nbF = 100;
        double pas = 10 / (double) nbF;
        System.out.println("pas= " + pas);
        double f_deb = 1e9;

        System.out.println("Nombre de frequences a analyser= " + nbF);
        double ff[] = new double[nbF];
        for (int i = 0; i < nbF; i++) {
            ff[i] = f_deb + i * pas;
        }

        Ecriture save = new Ecriture("D:/Zf_surf.txt");

        double ib[][];
        Matrix res;
        PEEC_RLMPC_ConducEquiv solP;
        VolumeRegion rV[] = new VolumeRegion[2];
        for (int j = 0; j < nbF; j++) {
            System.out.println((j + 1) + " / " + nbF);
            System.out.println("f= " + ff[j]);
            double f = ff[j];
            File fic = new File("");
            String path = fic.getAbsolutePath();
            String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/LA_SURF.DEC";

            ImportFlux IF = new ImportFlux(file);
            Region reg = IF.getRegion(0);
            reg.exportGmsh("d:/tmp.msh");
            ElementSurfSetHomogene mesh = new ElementSurfSetHomogene(reg.getElementSet().getElements());

            FaceRealDirichlet FluxNull = new FaceRealDirichlet(IF.getRegion(1), 0);
            FaceConstraint cons[] = new FaceConstraint[]{FluxNull, new ExternalDriveFaceConstraint(IF.getRegion(2)), new ExternalDriveFaceConstraint(IF.getRegion(3))};
            FaceDeg1 FS = new FaceDeg1(mesh, cons);

            PEECwithCapa_Full solPEEC = new PEECwithCapa_Full(FS, 2);
            // Affectation despoints de Gauss
            solPEEC.setPtsDeGaussInductifs(16, 4);
            solPEEC.setPtsDeGaussCapacitifs(4, 4);
            solPEEC.setRho(2.836e-8);
            solPEEC.setThickness(epaisseur);

            // Elmt hg = 2069 , bg = 724
            // Elmt hd = 35 , bd = 10
            BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

            int[] num = FS.getImplicitConstraintCount();
            int nbE = mesh.getNbElement();
            int num0 = num[0];
            int num1 = num[1];
            int nbBranches = solPEEC.getNbLignes();

            int indCapa = 1;
            int C1 = nbE + num0 + num1 + 1 + indCapa;
            int C2 = C1 + 1;
////        
            for (int in = 0; in < num0; in++) {
                circuitPur.addSourceUSimple(nbBranches + in, nbE + 1 + in, C1, "CommunS1", 0, 0.0);
            }
//
            for (int out = 0; out < num1; out++) {
                circuitPur.addSourceUSimple(nbBranches + num0 + out, nbE + 1 + num0 + out, C2, "CommunS2", 0, 0.0);
            }

            int iRef = nbBranches + num0 + num1;
            circuitPur.addSourceISimple(nbBranches + num0 + num1, C1, C2, "Source I", 1, 0);

            circuitPur.finSaisie(); // Methode a appeler necessairement a la fin de la saisie
            solPEEC.setCircuitElectrique(circuitPur);

            int nbActiveDof = FS.getActiveDofCount();

            // Resolution
            long tres0 = System.currentTimeMillis();
            ib = solPEEC.resolutionDirecte(f);
            long tres1 = System.currentTimeMillis();
            System.out.println("Temps de resolution:" + (tres1 - tres0));
            // On reorganise les dof   
            res = new Matrix(2, FS.getActiveDofCount());
            for (int i = 0; i < FS.getActiveDofCount(); i++) {
                res.setElement(0, solPEEC.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i] / epaisseur);
                res.setElement(1, solPEEC.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1] / epaisseur);
            }
            System.out.println("Z= " + ib[1][2 * iRef] + " + i* " + ib[1][2 * iRef + 1]);
            System.out.println("I= " + ib[0][2 * iRef] + " + i* " + ib[0][2 * iRef + 1]);
            save.ecrire(Math.hypot(ib[1][2 * iRef], ib[1][2 * iRef + 1]) + "\n");

            if (j <=2) {
                // On reorganise les dof   
                RealFaceQuantity Jreal = new RealFaceQuantity(FS, res.row(0));
                RealFaceQuantity Jimag = new RealFaceQuantity(FS, res.row(1));

                ExportGmshHdiv exportJreal = new ExportGmshHdiv(FS, "D:/jsiau/_Backup_Sources/Resultats/SNAKE/PS_RE_f" + f + "_nDof" + nbActiveDof + ".msh");
                exportJreal.addQuantity(Jreal, "Jreal");

                ExportGmshHdiv exportJimag = new ExportGmshHdiv(FS, "D:/jsiau/_Backup_Sources/Resultats/SNAKE/PS_IM_f" + f + "_nDof" + nbActiveDof + ".msh");
                exportJimag.addQuantity(Jimag, "Jimag");

                ComplexFaceQuantity J = new ComplexFaceQuantity(FS, res);
                exportJreal = new ExportGmshHdiv(FS, "D:/jsiau/_Backup_Sources/Resultats/SNAKE/PS_Jmod_f" + f + ".msh");
                exportJreal.addQuantityExportMod(J, "Jmod");

                ExportVtkHdiv ez = new ExportVtkHdiv(FS, "D:/jsiau/_Backup_Sources/Resultats/SNAKE/PS_f" + f + "_nDof" + nbActiveDof + ".vtk");
                ez.addQuantity(J, "J");

                String output_file = "D:/jsiau/_Backup_Sources/Resultats/SNAKE/PS_Jmod_f" + f + ".msh";
                File out_file = new File(output_file);
                File out_file1 = new File("D:/jsiau/_Backup_Sources/Resultats/SNAKE/PS_RE_f" + f + "_nDof" + nbActiveDof + ".msh");
                File out_file2 = new File("D:/jsiau/_Backup_Sources/Resultats/SNAKE/PS_IM_f" + f + "_nDof" + nbActiveDof + ".msh");
                try {
                    Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);
                } catch (IOException ex) {
                    Logger.getLogger(PEEC_RLMC_LOOPANTENA_SURF_Zf.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
        save.close();
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
