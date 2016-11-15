/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT.loopFace;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.LoopFace.FormulationLoopFaceBrutBeta;
import g2elab.mipse.material.AbstractMaterial;
import g2elab.mipse.material.conductor.LinearConductor;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.quantity.Quantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.numericalTools.iterativeSolver.Solver;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author jsiau
 */
public class TestLoopAntenna_LF {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().setNbCPU(4);
        // Lecture du fichier DEC puis recuperation du maillage
        File file = new File("");
        String path = file.getAbsolutePath();
        String nom = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/LOOPANTENNA_SURF3.DEC";

        ImportFlux imp = new ImportFlux(nom);
        Region[] regs = imp.getRegions();

        Region cond = regs[0];
        double ep = 2 * 1E-5;
//        double ep = 35e-6;
        ((ElementSurfSetHomogene) cond.getElementSet()).setThickness(ep);
        ((SurfaceRegion) cond).setThickness(ep);
        //
        double conduc = 1 / 1.72e-8;
        // Champ source

        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
        //
        FormulationLoopFaceBrutBeta sol = new FormulationLoopFaceBrutBeta();
        sol.setMeshGroup(new RegionsSet(new Region[]{cond}));
        sol.setMaterials(materials);
        sol.setCompression(Compression.No);
//        sol.setCompression("HCA");
        //
        //
//        sol.plotCircuit();
        //*
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
//        /*
        int b1[] = {273, 120, 233};
        int b2[] = {175, 140, 259};
        // Equipot 
        int nb1 = 10000;
        for (int i = 0; i < b1.length; i++) {
            circuit.addSourceUSimple(100000 + i, b1[i], nb1, "cct", 0.0, 0.0);
        }
        int nb2 = nb1 + 1;
        for (int i = 0; i < b2.length; i++) {
            circuit.addSourceUSimple(100000 + b1.length + i, b2[i], nb2, "cct", 0.0, 0.0);
        }
        //
        circuit.addSourceISimple(100000 + b1.length + b2.length, nb1, nb2, "src", 1.0, 0.0);
        /*/
         circuit.addSourceISimple(100000 , 273, 175, "src", 1.0, 0.0);
        
         //*/
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        /*/
         sol.setSourceField(new  SourceField[]{new UniformField(0, 0, 1)});
         //*/
        //
//        double[] f = Zf_PEEC_RLMPC_SURF_SERPENT.compteFrequenciesLog(9, 10, 2);
//        double[] f = new double[]{2.310129700083158E9};
//        sol.setFrequency(f);
        sol.setFrequency(new double[]{1e9});
//        sol.setFrequency(new double[]{1e9, 2e9, 3e9, 4e9, 5e9, 6e9, 7e9, 8e9, 9e9, 1e10});
        //
        sol.setSolver(Solver.FGMRes);
        ArrayList<Quantity> sortie = sol.solve();
        FormulationLoopFaceBrutBeta.exportRes(sortie, "d:/tmp/loopantenna/");
        sol.printCircuitInfos();

        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();

    }

}
