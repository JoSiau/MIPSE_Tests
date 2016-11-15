/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.FormulationLoopStar;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.FormulationTJ;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.MagnetodynamicTFormulation;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.material.AbstractMaterial;
import g2elab.mipse.material.conductor.LinearConductor;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshCell;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHgrad;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.Hgrad;
import g2elab.mipse.meshCore.quantity.Quantity;
import g2elab.mipse.meshCore.quantity.RealNodalQuantity;
import g2elab.mipse.meshCore.quantity.RealVectorCellQuantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.numericalTools.iterativeSolver.Solver;
import g2elab.mipse.tools.files.Ecriture;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class TestCircuit_LS {

    /**
     * Test uniquement le couplage circuit
     *
     * @param args
     */
    public static void main(String[] args) {
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + "\\src\\g2elab\\mipse\\formulationInProgress\\magnetodynamic\\FormulationT\\C_QUAD.DEC";
//       
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
//        ArrayList<MultiImport> arrIF = new ArrayList();
//        arrIF.add(IF);
//        ElementFactory efactory = new ElementFactory(arrIF, 1e-9, 1e-9, 1e-9);

//        RegionsSet regSet = new RegionsSet(IF);
//        ((ElementSurfSetHomogene) IF.getRegion(0).getElementSet()).computeVois();
        ((ElementSurfSetHomogene) IF.getRegion(0).getElementSet()).resetSurfOrientation();
        RegionsSet regSet = new RegionsSet(new Region[]{IF.getRegion(0)});
        double ep = 1e-4; // 0.1mm
        double conduc = 55e6;
        ((SurfaceRegion) regSet.getRegion(0)).setThickness(ep);
        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
        //
        FormulationLoopStar sol = new FormulationLoopStar();
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
//        /*
        sol.setCompression(Compression.No);
        /*/
         sol.setCompression("HCA");
         //*/
        sol.setFrequency(50);
        sol.DEBUGG_EXPORTMATRICES = false;
//        sol.plotCircuit();
//
        //
//        /*     
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
         //
        //
        circuit.addSourceISimple(100000, 0, 61, "src", 1.0, 0.0);
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        /*/
         UniformMagneticField srcField = new UniformMagneticField(0.0, 0.0, 1.0); //1.0);
         sol.setSources(new SourceField[]{srcField});
         //*/
        //
        //
        //
        sol.setSolver(Solver.LU);
        ArrayList<Quantity> sortie = sol.solve();
        RealNodalQuantity Treal = null, Timag = null;
        RealVectorCellQuantity Jreal = null, Jimag = null;
        for (Quantity sortie1 : sortie) {
            switch (sortie1.getName()) {
                case "Treal":
                    Treal = (RealNodalQuantity) sortie1;
                    break;
                case "Timag":
                    Timag = (RealNodalQuantity) sortie1;
                    break;
                case "JtotReal":
                    Jreal = (RealVectorCellQuantity) sortie1;
                    break;
                case "JtotImag":
                    Jimag = (RealVectorCellQuantity) sortie1;
                    break;
                default:
                    System.err.println("Problem, quantity not identified !");
                    break;
            }
        }
        //
        Pertes p = new Pertes(Jreal.getElementSet());
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
            merge.ecrire("Plugin(MathEval).View = 0;\n"
                    + "Plugin(MathEval).Run;\n"
                    + "Plugin(MathEval).View = 1;\n"
                    + "Plugin(MathEval).Run;");
            merge.close();
            File out_file = new File(pathProspro + "merge.msh");
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file);
        } catch (IOException ex) {
            Logger.getLogger(MagnetodynamicTFormulation.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

    public static void mains(String[] args) {
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere/SPHERE_572.DEC";
//       
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
//        ArrayList<MultiImport> arrIF = new ArrayList();
//        arrIF.add(IF);
//        ElementFactory efactory = new ElementFactory(arrIF, 1e-9, 1e-9, 1e-9);

//        RegionsSet regSet = new RegionsSet(IF);
//        ((ElementSurfSetHomogene) IF.getRegions(0).getElementSet()).computeVois();
        ((ElementSurfSetHomogene) IF.getRegion(0).getElementSet()).resetSurfOrientation();
        RegionsSet regSet = new RegionsSet(new Region[]{IF.getRegion(0)});
        double ep = 1e-4; // 0.1mm
        double conduc = 55e6;
        ((SurfaceRegion) regSet.getRegion(0)).setThickness(ep);
        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
        //
        FormulationLoopStar sol = new FormulationLoopStar();
//        FormulationLoopFace sol = new FormulationLoopFace();
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
//        /*
        sol.setCompression(Compression.No);
        /*/
         sol.setCompression("HCA");
         //*/
        sol.setFrequency(1e3);
//        sol.plotCircuit();
//
        //
//        /*     
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
         //
        //
        circuit.addSourceISimple(100000, 71, 499, "src", 1.0, 0.0);
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        /*/
         UniformMagneticField srcField = new UniformMagneticField(0.0, 0.0, 1.0); //1.0);
         sol.setSources(new SourceField[]{srcField});
         //*/
        //
        //
        //
        sol.setSolver(Solver.LU);
        ArrayList<Quantity> sortie = sol.solve();
        RealNodalQuantity Treal = null, Timag = null;
        RealVectorCellQuantity Jreal = null, Jimag = null;
        for (Quantity sortie1 : sortie) {
            switch (sortie1.getName()) {
                case "Treal":
                    Treal = (RealNodalQuantity) sortie1;
                    break;
                case "Timag":
                    Timag = (RealNodalQuantity) sortie1;
                    break;
                case "JtotReal":
                    Jreal = (RealVectorCellQuantity) sortie1;
                    break;
                case "JtotImag":
                    Jimag = (RealVectorCellQuantity) sortie1;
                    break;
                default:
                    System.err.println("Problem, quantity not identified !");
                    break;
            }
        }
        //
        Pertes p = new Pertes(Jreal.getElementSet());
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
            merge.ecrire("Plugin(MathEval).View = 0;\n"
                    + "Plugin(MathEval).Run;\n"
                    + "Plugin(MathEval).View = 1;\n"
                    + "Plugin(MathEval).Run;");
            merge.close();
            File out_file = new File(pathProspro + "merge.msh");
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file);
        } catch (IOException ex) {
            Logger.getLogger(MagnetodynamicTFormulation.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

    /**
     * Test le circuit coupled avec les loops
     *
     * @param args
     */
    public static void mainqq(String[] args) {
        System.out.println("Test avec trou !");
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/4HOLES_QUAD.DEC";
//       
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
//        ArrayList<MultiImport> arrIF = new ArrayList();
//        arrIF.add(IF);
//        ElementFactory efactory = new ElementFactory(arrIF, 1e-9, 1e-9, 1e-9);

        RegionsSet regSet = new RegionsSet(new Region[]{IF.getRegion(0)});
        double ep = 1e-4; // 1mm
        double conduc = 1 / 1.72e-8;
        ((SurfaceRegion) regSet.getRegion(0)).setThickness(ep);
        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
        //
        FormulationTJ sol = new FormulationTJ();
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
        sol.setCompression(Compression.No);
//        sol.setCompression("HCA");
        sol.setFrequency(100);
        //
        //
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
        // 4 trous
        circuit.addSourceUSimple(10000, 390, 381, "source", 1.0, 0.0);
        circuit.addSourceUSimple(10001, 399, 372, "source", 1.0, 0.0);
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        //
        //
        //
        ArrayList<Quantity> sortie = sol.solve();
        RealNodalQuantity Treal = null, Timag = null;
        RealVectorCellQuantity Jreal = null, Jimag = null;
        for (Quantity sortie1 : sortie) {
            switch (sortie1.getName()) {
                case "Treal":
                    Treal = (RealNodalQuantity) sortie1;
                    break;
                case "Timag":
                    Timag = (RealNodalQuantity) sortie1;
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
        Pertes p = new Pertes(Jreal.getElementSet());
        double loss = p.calcul(Jreal, Jimag, conduc, ep, 9);
        System.out.println("\nPertes = " + loss + "\n");
        //
        String pathProspro = "D:/tmp/";
        ExportGmshHgrad expTR = new ExportGmshHgrad((Hgrad) Treal.getFS(), pathProspro + "Ttotreal.msh");
        expTR.addQuantity(Treal, "Ttotreal");
        ExportGmshHgrad expTI = new ExportGmshHgrad((Hgrad) Timag.getFS(), pathProspro + "Ttotimag.msh");
        expTI.addQuantity(Timag, "Ttotimag");
        ExportGmshCell expJR = new ExportGmshCell((Cell) Jreal.getFS(), pathProspro + "Jreal.msh");
        expJR.addQuantity(Jreal, "Jreal");
        ExportGmshCell expJI = new ExportGmshCell((Cell) Jimag.getFS(), pathProspro + "Jimag.msh");
        expJI.addQuantity(Jimag, "Jimag");
        //
        try {
            Ecriture merge = new Ecriture(pathProspro + "merge.msh");
            merge.ecrire("Merge '" + pathProspro + "Ttotreal.msh';\n");
            merge.ecrire("Merge '" + pathProspro + "Ttotimag.msh';\n");
            merge.ecrire("Merge '" + pathProspro + "Jreal.msh';\n");
            merge.ecrire("Merge '" + pathProspro + "Jimag.msh';\n");
            if (MagnetodynamicTFormulation.DEBUGG) {
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
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
