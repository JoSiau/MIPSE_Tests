/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.FormulationTJ;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.FormulationTJ_alpha;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.MagnetodynamicTFormulation;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.material.AbstractMaterial;
import g2elab.mipse.material.conductor.LinearConductor;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshCell;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.Quantity;
import g2elab.mipse.meshCore.quantity.RealVectorCellQuantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.numericalTools.iterativeSolver.Solver;
import g2elab.mipse.tools.files.Ecriture;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The type Test volumic.
 * @author jsiau
 */
public class TestVolumic {

    /**
     * The entry point of application.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        System.out.println("Test volumique avec trou !");
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + "/src/formulations/formulationT/PLAQUE_TROU_TETRA_4350.DEC";
//        path = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/PLAQUE_VOL.DEC";
//        path = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/LOOPANTENA.DEC";
//       
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
//        ArrayList<MultiImport> arrIF = new ArrayList();
//        arrIF.add(IF);
//        ElementFactory efactory = new ElementFactory(arrIF, 1e-9, 1e-9, 1e-9);

//        RegionsSet regSet = new RegionsSet(IF);
//        ((ElementSurfSetHomogene) IF.getRegion(0).getElementSet()).computeVois();
//        ((ElementSurfSetHomogene) IF.getRegion(0).getElementSet()).resetSurfOrientation();
        RegionsSet regSet = new RegionsSet(new Region[]{IF.getRegion(0)});
        double ep = 1.0;
//        double ep = 1e-4;
        double conduc = 1 / 1.68e-8;
//        double conduc = 5.5e7;
        // Champ source

        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
        //
        FormulationTJ sol = new FormulationTJ();
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
        sol.setCompression(Compression.No);
        sol.setCapacitiveComputation(false);
        sol.setSolver(Solver.FGMRes);
//        sol.setCompression("HCA");
        sol.setFrequency(100);
        //
        // ANTENNA
//        sol.addEquiPotentialRegion((SurfaceRegion) IF.getRegion(2));
//        sol.addEquiPotentialRegion((SurfaceRegion) IF.getRegion(3));
        //
//        /*
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
        circuit.addSourceISimple(1000000000, 3552, 4348, "Source", 1.0, 0.0);// PLAQUE_TROU_TETRA_4350
//        circuit.addSourceISimple(1000000000, 422, 657, "Source", 1.0, 0.0);// PLAQUE VOL
//        circuit.addSourceISimple(1000000000, regSet.getElTotalCount(), regSet.getElTotalCount() + 1, "Source", 1.0, 0.0);// ANTENNA
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        /*/
         sol.setSourceField(new  SourceField[]{new UniformField(0, 0, 1)});
         //*/
        sol.plotCircuit();
        //
        //
        //
        ArrayList<Quantity> sortie = sol.solve();
//        RealNodalQuantity Treal = null, Timag = null;
        RealVectorCellQuantity Jreal = null, Jimag = null;
        for (Quantity sortie1 : sortie) {
            switch (sortie1.getName()) {
                case "Treal":
//                    Treal = (RealNodalQuantity) sortie1;
                    break;
                case "Timag":
//                    Timag = (RealNodalQuantity) sortie1;
                    break;
                case "Jreal":
                    Jreal = (RealVectorCellQuantity) sortie1;
                    break;
                case "Jimag":
                    Jimag = (RealVectorCellQuantity) sortie1;
                    break;
                default:
                    System.err.println("Problem, quantity not identified !");
                    break;
            }
        }
        //
        String pathProspro = "D:/tmp/";
//        ExportGmshHgrad expTR = new ExportGmshHgrad((Hgrad) Treal.getFS(), pathProspro + "Ttotreal.msh");
//        expTR.addQuantity(Treal, "Ttotreal");
//        ExportGmshHgrad expTI = new ExportGmshHgrad((Hgrad) Timag.getFS(), pathProspro + "Ttotimag.msh");
//        expTI.addQuantity(Timag, "Ttotimag");
        ExportGmshCell expJR = new ExportGmshCell((Cell) Jreal.getFS(), pathProspro + "Jreal.msh");
        expJR.addQuantity(Jreal, "Jreal");
        ExportGmshCell expJI = new ExportGmshCell((Cell) Jimag.getFS(), pathProspro + "Jimag.msh");
        expJI.addQuantity(Jimag, "Jimag");
        //
        try {
            Ecriture merge = new Ecriture(pathProspro + "merge.msh");
//            merge.ecrire("Merge '" + pathProspro + "Ttotreal.msh';\n");
//            merge.ecrire("Merge '" + pathProspro + "Ttotimag.msh';\n");
            merge.ecrire("Merge '" + pathProspro + "Jreal.msh';\n");
            merge.ecrire("Merge '" + pathProspro + "Jimag.msh';\n");
            if (MagnetodynamicTFormulation.DEBUGG && !sol.isIsCapa()) {
                for (int i = 0; i < sol.getNbMILoops(); i++) {
                    merge.ecrire("Merge '" + pathProspro + "Jloop" + i + "Real.msh';\n");
                    merge.ecrire("Merge '" + pathProspro + "Jloop" + i + "Imag.msh';\n");
                }
            }
            merge.close();
            File out_file = new File(pathProspro + "merge.msh");
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file);
        } catch (IOException ex) {
            Logger.getLogger(MagnetodynamicTFormulation.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        //
        Pertes p = new Pertes(Jreal.getElementSet());
//        double loss = p.calcul(Jreal, Jimag, conduc, ep, 8);
        double loss = p.calcul(Jreal, Jimag, conduc, ep, 1);
        System.out.println("\nPertes = " + loss + "\n");
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
