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
import g2elab.mipse.meshCore.quantity.Quantity;
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
public class TestSerpent_LF {

    public static void main(String[] args) {
        System.out.println("Test avec trou !");
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path += "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Surf/S3_NB_FS.DEC";
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
        RegionsSet regSet = new RegionsSet(new Region[]{IF.getRegion(0)});
        double ep = 35e-6;
        ((SurfaceRegion) regSet.getRegion(0)).setThickness(ep);
        // Matériaux
        double conduc = 1 / 1.68e-8;
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
        //
        FormulationLoopFaceBrutBeta sol = new FormulationLoopFaceBrutBeta();
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
        sol.setCompression(Compression.No);
//        sol.setCompression("HCA");
//        double f = 4.8e7;
        double f[] = getF();
        sol.setFrequency(f);
        sol.setNormalization(false);
//        sol.plotCircuit();
        //
        //
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
        circuit.addSourceISimple(1000000, 404, 809, "Source", 1.0, 0.0);
        circuit.addSourceUSimple(1000001, 297, 702, "cct", 0, 0);
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        //
        //
        //
        sol.setSolver(Solver.LU);
        sol.setPrecond(1);
        ArrayList<Quantity> sortie = sol.solve();
        String pathProspro = "D:/tmp/";
        sol.exportRes(sortie, pathProspro);
        sol.printCircuitInfos();
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

    public static double[] getF() {
        return new double[]{1000000, 1291549.665, 1668100.537, 2154434.69, 2782559.402,
            3593813.664, 4641588.834, 5994842.503, 7742636.827, 1.00E+07, 1.02E+07, 1.05E+07,
            1.07E+07, 1.10E+07, 1.12E+07, 1.15E+07, 1.18E+07, 1.20E+07, 1.23E+07, 1.26E+07,
            1.29E+07, 1.32E+07, 1.35E+07, 1.38E+07, 1.42E+07, 1.45E+07, 1.48E+07, 1.52E+07,
            1.56E+07, 1.59E+07, 1.63E+07, 1.67E+07, 1.71E+07, 1.75E+07, 1.79E+07, 1.83E+07,
            1.87E+07, 1.92E+07, 1.96E+07, 2.01E+07, 2.06E+07, 2.10E+07, 2.15E+07, 2.21E+07,
            2.26E+07, 2.31E+07, 2.36E+07, 2.42E+07, 2.48E+07, 2.54E+07, 2.60E+07, 2.66E+07,
            2.72E+07, 2.78E+07, 2.85E+07, 2.92E+07, 2.98E+07, 3.05E+07, 3.13E+07, 3.20E+07,
            3.27E+07, 3.35E+07, 3.43E+07, 3.51E+07, 3.59E+07, 3.68E+07, 3.76E+07, 3.85E+07,
            3.94E+07, 4.04E+07, 4.13E+07, 4.23E+07, 4.33E+07, 4.43E+07, 4.53E+07, 4.64E+07,
            4.75E+07, 4.86E+07, 4.98E+07, 5.09E+07, 5.21E+07, 5.34E+07, 5.46E+07, 5.59E+07,
            5.72E+07, 5.86E+07, 5.99E+07, 6.14E+07, 6.28E+07, 6.43E+07, 6.58E+07, 6.73E+07,
            6.89E+07, 7.05E+07, 7.22E+07, 7.39E+07, 7.56E+07, 7.74E+07, 7.92E+07, 8.11E+07,
            8.30E+07, 8.50E+07, 8.70E+07, 8.90E+07, 9.11E+07, 9.33E+07, 9.55E+07, 9.77E+07, 1.00E+08,
            2e8, 3e8, 4e8, 5e8, 6e8, 7e8, 8e8, 9e8, 1e9
        };
    }

}
//