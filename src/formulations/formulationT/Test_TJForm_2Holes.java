/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT;

import g2elab.mipse.analytical.sourceFields.SourceField;
import g2elab.mipse.analytical.sourceFields.UniformMagneticField;
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
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.tools.files.Ecriture;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * The type Test _ tJ form _ 2 holes.
 * @author jsiau
 */
public class Test_TJForm_2Holes {

    /**
     * The entry point of application.
     *
     * @param args the input arguments
     */
    public static void main(String[] args) {
        System.out.println("Test avec trou !");
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
//        path = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/2HOLES_COIL03_TRI_B.DEC";
        path = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/FormulationT/4HOLES_TRI.DEC";
//        path = path + "/src/g2elab/mipse/stabilizedFormulations/formulationsGM/TOR_BIGMESH_IMPSURF_COIL12.DEC";CONTACTEUR_MAGNET_TRI_2094.DEC
//        path = path + "/src/g2elab/mipse/stabilizedFormulations/formulationsGM/fichiersDEC/PLAQUE_COIL03_TRI_1020.DEC";
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
        double ep = 1e-4;
        double conduc = 55e6;
        ((SurfaceRegion) regSet.getRegion(0)).setThickness(ep);
        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //*
        FormulationTJ_alpha sol = new FormulationTJ_alpha();
        /*/
        FormulationTJ sol = new FormulationTJ();
        //*/
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
        sol.setCompression(Compression.No);
//        sol.setCompression("HCA");
        sol.setFrequency(100);    
        //* 
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
        // 2 trous
//        circuit.addSourceISimple(10000, 416, 9, "source", 1.0, 0.0);
//        circuit.addSourceISimple(10001, 417, 414, "source", 1.0, 0.0);
        // 4 trous
        circuit.addSourceISimple(10000, 515, 509, "source", 1.0, 0.0);
        circuit.addSourceISimple(10001, 514, 219, "source", 1.0, 0.0);
        
        
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        /*/
        sol.setSources(new SourceField[]{srcField});
        //*/
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

//        ComplexMatrixAssembly b = (ComplexMatrixAssembly) sol.matrixZ.getBlock(1, 1);
//
//        double mat[][] = new double[2][4];
//
//        mat[0] = b.product(new double[]{1.0, 0.0, 0.0, 0.0}, mat[0]);
//        mat[1] = b.product(new double[]{0.0, 0.0, 1.0, 0.0}, mat[1]);
//
//        System.out.println("mf:\n" + (new Basic2D(mat)).transpose());
//
//        Basic2D magg = TJ_Preconditioner.agglomerate(b);
//        System.out.println("magg:\n" + magg);
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

}