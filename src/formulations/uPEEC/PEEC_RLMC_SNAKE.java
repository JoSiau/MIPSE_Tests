/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.IO.paraview.ExportVtkHdiv;
import g2elab.mipse.meshCore.contraints.ExternalDriveFaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceConstraint;
import g2elab.mipse.meshCore.elements.ElementVolSetHomogene;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class PEEC_RLMC_SNAKE {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        double f = 100;

        File fic = new File("");
        String path = fic.getAbsolutePath();
        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/SERPENT_3C.DEC";

        ImportFlux IF = new ImportFlux(file);
        Region reg = IF.getRegion(0);
        reg.exportGmsh("d:/tmp.msh");
        ElementVolSetHomogene mesh = new ElementVolSetHomogene(reg.getElementSet().getElements());

        RegionsSet regSet = new RegionsSet(IF.getRegions());
        SurfaceRegion border = (SurfaceRegion) regSet.generateBorder(IF.getRegion(0));
        regSet.addRegion(border);
        border = regSet.generateSurfaceComplement(border, IF.getRegion(1), "tmp");
        regSet.addRegion(border);
        border = regSet.generateSurfaceComplement(border, IF.getRegion(2), "border w/o borne");

        FaceConstraint cons[] = new FaceConstraint[]{new ExternalDriveFaceConstraint(border), new ExternalDriveFaceConstraint(IF.getRegion(1)), new ExternalDriveFaceConstraint(IF.getRegion(2))};
        FaceDeg1 FS = new FaceDeg1(mesh, cons);

        PEECwithCapa_Full solPEEC = new PEECwithCapa_Full(FS, 2);
        // Affectation despoints de Gauss
        solPEEC.setPtsDeGaussInductifs(8, 8);
        solPEEC.setPtsDeGaussCapacitifs(4, 4);

        // Elmt hg = 2069 , bg = 724
        // Elmt hd = 35 , bd = 10
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        int[] num = FS.getImplicitConstraintCount();
        int nbE = mesh.getNbElement();
        int nbCapa = num[0]; //solPEEC.getSizeC();
        int num0 = num[1];
        int num1 = num[2];

        int nbBranches = solPEEC.getNbLignes();
////     noeud commun des capas (noeud GND)   
        int indCapa = 1;
        int C1 = (nbE + 1) + nbCapa + num0 + num1 + indCapa;
        int C2 = C1 + 1;
////        
        for (int in = 0; in < num0; in++) {
            circuitPur.addSourceUSimple(nbBranches + in, nbE + nbCapa + 1 + in, C1, "CommunS1", 0, 0.0);
        }
//
        for (int out = 0; out < num1; out++) {
            circuitPur.addSourceUSimple(nbBranches + num0 + out, nbE + nbCapa + 1 + num0 + out, C2, "CommunS2", 0, 0.0);
        }

        circuitPur.addSourceISimple(nbBranches + num0 + num1, C2, C1, "SourceI", 1, 0.0);
        int iRef = nbBranches + num0 + num1;

        circuitPur.finSaisie(); // Methode a appeler necessairement a la fin de la saisie

        circuitPur.toString();
        // Affectation du circuit exterieur
        solPEEC.setCircuitElectrique(circuitPur);

        int nbActiveDof = FS.getActiveDofCount();

        // Resolution
        long tres0 = System.currentTimeMillis();
        double ib[][] = solPEEC.resolutionDirecte(f);
        long tres1 = System.currentTimeMillis();
        System.out.println("Temps de resolution:" + (tres1 - tres0));
        // On reorganise les dof   
        Matrix res = new Matrix(2, FS.getActiveDofCount());
        for (int i = 0; i < FS.getActiveDofCount(); i++) {
            res.setElement(0, solPEEC.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solPEEC.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }
        System.out.println("Z= " + ib[1][2 * iRef] + " + i* " + ib[1][2 * iRef + 1]);
        System.out.println("I= " + ib[0][2 * iRef] + " + i* " + ib[0][2 * iRef + 1]);

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
            Logger.getLogger(PEEC_RLMC_SNAKE.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

}
