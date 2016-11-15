/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.FormulationTJ_alpha;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.MagnetodynamicTFormulation;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.material.AbstractMaterial;
import g2elab.mipse.material.conductor.LinearConductor;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshCell;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHgrad;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.Hgrad;
import g2elab.mipse.meshCore.quantity.Quantity;
import g2elab.mipse.meshCore.quantity.RealNodalQuantity;
import g2elab.mipse.meshCore.quantity.RealVectorCellQuantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.tools.files.Ecriture;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv.eps0;
import g2elab.mipse.mipseCore.matrixCompression.Compression;

/**
 * The type Test surf pP.
 * @author jsiau
 */
public class TestSurfPP {

    /**
     * Test sur une capa pour le calcul de la capa plan.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.out.println("Test avec trou !");
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/PP100X100X1_800E.DEC";
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
        RegionsSet regSet = new RegionsSet(new Region[]{IF.getRegion(0)});
        double ep = 1e-6; // 1mm
        double conduc = 55e6;
        ((SurfaceRegion) regSet.getRegion(0)).setThickness(ep);
        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
//        /*
        FormulationTJ_alpha sol = new FormulationTJ_alpha();
        /*/
         FormulationTJ sol = new FormulationTJ();
         //*/
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
        sol.setCompression(Compression.No);
        //
        double f = 100;
        sol.setFrequency(f);
        sol.plotCircuit();
        //
        //
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
        int b1[] = {112, 116, 126, 99, 129, 137, 136, 194, 179, 188};
        int b2[] = {602, 596, 575, 579, 578, 580, 586, 591, 553, 594};
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
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
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
            merge.ecrire("Merge '" + pathProspro + "Jreal.msh';\n");
            merge.ecrire("Merge '" + pathProspro + "Jimag.msh';\n");
            if (MagnetodynamicTFormulation.DEBUGG) {
                merge.ecrire("Merge '" + pathProspro + "JloopReal.msh';\n");
                merge.ecrire("Merge '" + pathProspro + "JloopImag.msh';\n");
            }
            merge.close();
            File out_file = new File(pathProspro + "merge.msh");
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file);
        } catch (IOException ex) {
            Logger.getLogger(MagnetodynamicTFormulation.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        //
        double modZ = Math.hypot(sol.getResCircuit().get(f).getElement(0, 0), sol.getResCircuit().get(f).getElement(0, 1));
        System.out.println("mod(Z) = " + modZ);
        double omega = 2 * Math.PI * f;
        double cCom = (1.0 / (omega * modZ));
        System.out.println("C_comp = 1/(w * |Z|) = " + cCom);
        double C = eps0 * 1e-2 / 1e-3;
        System.out.println("C_ana = eps0 * S / e = " + C);
        //
        System.out.println("Erreur relative = " + (Math.abs(C - cCom) / C));

        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

    /**
     * Test sur une capa avec le maillage minimal.
     *
     * @param args the command line arguments
     */
    public static void main4(String[] args) {
        System.out.println("Test avec trou !");
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/PP100X100X10_SIMPLE.DEC";
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
        RegionsSet regSet = new RegionsSet(new Region[]{IF.getRegion(0)});
        double ep = 1e-3; // 1mm
        double conduc = 55e6;
        ((SurfaceRegion) regSet.getRegion(0)).setThickness(ep);
        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
//        /*
        FormulationTJ_alpha sol = new FormulationTJ_alpha();
        /*/
         FormulationTJ sol = new FormulationTJ();
         //*/
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
        sol.setCompression(Compression.No);
        sol.setFrequency(1e9);
//        sol.plotCircuit();
        //
        //
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
        circuit.addSourceISimple(1000, 0, 4, "src", 1.0, 0.0);
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
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
            merge.ecrire("Merge '" + pathProspro + "Jreal.msh';\n");
            merge.ecrire("Merge '" + pathProspro + "Jimag.msh';\n");
            if (MagnetodynamicTFormulation.DEBUGG) {
                merge.ecrire("Merge '" + pathProspro + "JloopReal.msh';\n");
                merge.ecrire("Merge '" + pathProspro + "JloopImag.msh';\n");
            }
            merge.close();
            File out_file = new File(pathProspro + "merge.msh");
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file);
        } catch (IOException ex) {
            Logger.getLogger(MagnetodynamicTFormulation.class
                    .getName()).log(Level.SEVERE, null, ex);
        }
        System.out.println("Capacitance currents:");
        double[] ib = sol.getCapaCurrents();
        for (int i = 0; i < ib.length / 2; i++) {
            System.out.println("capa " + i + " = " + ib[2 * i] + " + j * " + ib[2 * i + 1]);
        }
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
