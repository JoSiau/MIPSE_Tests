/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT.cefc2016;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.FormulationTJ_alpha;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.MagnetodynamicTFormulation;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.material.AbstractMaterial;
import g2elab.mipse.material.conductor.LinearConductor;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshCell;
import g2elab.mipse.meshCore.elements.ElementFactory;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.Quantity;
import g2elab.mipse.meshCore.quantity.RealVectorCellQuantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.numericalTools.iterativeSolver.Solver;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.files.Exec;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class Test_TIm {

    public static double f = 100;
    public static double ep = 1e-4;
    public static double conduc = 55e6;
//    public static String mesh = "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/4HOLES_TRI.DEC";
    public static String mesh = "/src/formulations/formulationT/cefc2016/TOR_SURF.DEC";

    public static void main(String[] args) {
        ElementFactory ef = new ElementFactory(ElementFactory.STEP, ElementFactory.STEP, ElementFactory.STEP);
        System.out.println("Test avec trou !");
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + mesh;
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
        RegionsSet regSet = new RegionsSet(new Region[]{IF.getRegions()[0]});
        ((SurfaceRegion) regSet.getRegion(0)).setThickness(ep);

        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
        //
        ////////////////////////////////////////////////////////////////////////
        FormulationTJ_alpha sol = new FormulationTJ_alpha();
        sol.DEBUGG_EXPORTMATRICES = false;
        sol.DEBUGG = false;
        sol.setPrecond(1);
        //
        //
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
        sol.setCompression(Compression.HCA);
        sol.setCapacitiveComputation(false);

        sol.setFrequency(f);

        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
        // tor
        circuit.addSourceISimple(1000000, 960, 1679, "source", 1.0, 0.0);
        // 4 trous
//        circuit.addSourceISimple(10000, 515, 509, "source", 1.0, 0.0);
//        circuit.addSourceISimple(10001, 514, 219, "source", 1.0, 0.0);

        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);

        sol.setSolver(Solver.FGMRes);
//        sol.setSolver(Solver.LU);
        ArrayList<Quantity> sortie = sol.solve();
        RealVectorCellQuantity Jreal = null, Jimag = null;
        for (Quantity sortie1 : sortie) {
            switch (sortie1.getName()) {
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
        Pertes p = new Pertes(Jreal.getElementSet());
//        Jreal.scale(ep);
//        Jreal.scale(ep);
        double loss = p.calcul(Jreal, Jimag, conduc, ep, 9);
        System.out.println("\nPertes = " + loss + "\n");
        //
        String pathProspro = "D:/tmp/";
        ExportGmshCell expJR = new ExportGmshCell((Cell) Jreal.getFS(), pathProspro + "Jreal.msh");
        expJR.addQuantity(Jreal, "Jreal");
        ExportGmshCell expJI = new ExportGmshCell((Cell) Jimag.getFS(), pathProspro + "Jimag.msh");
        expJI.addQuantity(Jimag, "Jimag");
        //
        try {
            Ecriture merge = new Ecriture(pathProspro + "merge.msh");
            merge.ecrire("Merge '" + pathProspro + "Jreal.msh';\n");
            merge.ecrire("Merge '" + pathProspro + "Jimag.msh';\n");
            merge.close();
            Exec.openFile(pathProspro + "merge.msh");
        } catch (IOException ex) {
            Logger.getLogger(MagnetodynamicTFormulation.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
