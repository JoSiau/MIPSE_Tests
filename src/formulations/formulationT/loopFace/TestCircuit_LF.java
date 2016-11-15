/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT.loopFace;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.LoopFace.FormulationLoopFaceBrutBeta;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.material.AbstractMaterial;
import g2elab.mipse.material.conductor.LinearConductor;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.quantity.Quantity;
import g2elab.mipse.meshCore.quantity.RealNodalQuantity;
import g2elab.mipse.meshCore.quantity.RealVectorCellQuantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.numericalTools.iterativeSolver.Solver;

import java.io.File;
import java.util.ArrayList;

/**
 *
 * @author jsiau
 */
public class TestCircuit_LF {

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
//        FormulationLoopStar sol = new FormulationLoopStar();
        FormulationLoopFaceBrutBeta sol = new FormulationLoopFaceBrutBeta();
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
        sol.setSolver(Solver.FGMRes);
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
        FormulationLoopFaceBrutBeta.exportRes(sortie, "d:/tmp/");
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
