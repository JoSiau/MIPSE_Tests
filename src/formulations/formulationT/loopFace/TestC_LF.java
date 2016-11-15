/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT.loopFace;

import g2elab.mipse.analytical.sourceFields.SourceField;
import g2elab.mipse.analytical.sourceFields.UniformMagneticField;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.LoopFace.FormulationLoopFaceBrutBeta;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.material.AbstractMaterial;
import g2elab.mipse.material.conductor.LinearConductor;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.quantity.Quantity;
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
public class TestC_LF {

    /**
     * Test uniquement le couplage circuit
     *
     * @param args
     */
    public static void main(String[] args) {
        File file = new File("");
        String path;
        path = file.getAbsolutePath();

        path = path + "\\src\\g2elab\\mipse\\formulationInProgress\\magnetodynamic\\U_PEEC_DIELECTRIC\\C_MED.DEC";
//       
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
        
        RegionsSet regSet = new RegionsSet(new Region[]{IF.getRegion(0)});
        double ep = 1e-4; // 0.1mm
        double conduc = 55e6;
        ((SurfaceRegion) regSet.getRegion(0)).setThickness(ep);
        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
        //
        FormulationLoopFaceBrutBeta sol = new FormulationLoopFaceBrutBeta();
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
//        /*
        sol.setCompression(Compression.No);
        /*/
         sol.setCompression("HCA");
         //*/
        sol.DEBUGG_EXPORTMATRICES = false;
        
        double f = 1e3;
        sol.setFrequency(f);
//        sol.plotCircuit();
//
        //
        /*     
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
        sol.setSolver(Solver.SOR);
        ArrayList<Quantity> sortie = sol.solve();
        RealVectorCellQuantity Jreal = null, Jimag = null;
        for (Quantity sortie1 : sortie) {
            switch (sortie1.getName()) {
                case "JtotReal":
                    Jreal = (RealVectorCellQuantity) sortie1;
                    break;
                case "JtotImag":
                    Jimag = (RealVectorCellQuantity) sortie1;
                    break;
                default:
                    break;
            }
        }
        //
        Pertes p = new Pertes(Jreal.getElementSet());
        double loss = p.calcul(Jreal, Jimag, conduc, ep, 9);
        System.out.println("\nPertes = " + loss + "\n");
        //
        FormulationLoopFaceBrutBeta.exportRes(sortie, "d:/tmp/C/" + f);
//        sol.exportRes(f, "d:/tmp/C/b");
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
