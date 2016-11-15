/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.FormulationTJ_alpha;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.MagnetodynamicTFormulation;
import g2elab.mipse.material.AbstractMaterial;
import g2elab.mipse.material.conductor.LinearConductor;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshCell;
import g2elab.mipse.meshCore.IO.paraview.ExportVtkCell;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.Quantity;
import g2elab.mipse.meshCore.quantity.RealVectorCellQuantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.numericalTools.iterativeSolver.Solver;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class TestBobineRC {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        GestionnaireTaches.getGestionnaireTaches().setNbCPU(8);

        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(TestBobineRC.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportFlux IF = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/ArticlesCodesJS/FormulationVolumique/RC_1SPIRE.DEC");
        //
        // Gather the meshes
        //
        // Conductor
        VolumeRegion conductor = (VolumeRegion) IF.getRegion(0);
        RegionsSet regSet = new RegionsSet(new Region[]{conductor});
        // Put the  2 dieledctric regions in an array
//        VolumeRegion dielectrics[] = new VolumeRegion[]{(VolumeRegion) mesh.getRegion(1)};
        //
        SurfaceRegion fluxPos = (SurfaceRegion) IF.getRegion(2);
        SurfaceRegion fluxNeg = (SurfaceRegion) IF.getRegion(3);

        double sigma = 1 / 1.8e-5;
        double conduc = 55555.55555555;
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
        
        
        
        sol.setSolver(Solver.FGMRes);
        sol.setPrecond(1);
        
        
        
        
        //
        double f = 1e9;
        sol.setFrequency(f);
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();

//        //
        sol.addEquiPotentialRegion((SurfaceRegion) IF.getRegion(2));
        sol.addEquiPotentialRegion((SurfaceRegion) IF.getRegion(3));
        circuit.addSourceISimple(1000000, IF.getRegion(0).getEltCount(), IF.getRegion(0).getEltCount() + 1, "src", 1.0, 0.0);

//        sol.addEquiPotentialRegion((SurfaceRegion) IF.getRegion(1));
//        sol.addEquiPotentialRegion((SurfaceRegion) IF.getRegion(2));
//        sol.addEquiPotentialRegion((SurfaceRegion) IF.getRegion(3));
//        sol.addEquiPotentialRegion((SurfaceRegion) IF.getRegion(4));
//        circuit.addSourceISimple(1000000, IF.getRegion(0).getEltCount(), IF.getRegion(0).getEltCount() + 2, "src", 1.0, 0.0);
//        circuit.addSourceISimple(1000001, IF.getRegion(0).getEltCount()+1, IF.getRegion(0).getEltCount() + 3, "src", 1.0, 0.0);
        //
//        sol.plotCircuit();
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        //
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
            }
        }
        //
//        Pertes p = new Pertes(Jreal.getElementSet());
//        double loss = p.calcul(Jreal, Jimag, conduc, 1.0, 8);
//        System.out.println("\nPertes = " + loss + "\n");
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

        ExportVtkCell expJRvtk = new ExportVtkCell((Cell) Jreal.getFS(), pathProspro + "Jreal.vtk");
        expJRvtk.addQuantity(Jreal, "Jreal");
        ExportVtkCell expJIvtk = new ExportVtkCell((Cell) Jimag.getFS(), pathProspro + "Jimag.vtk");
        expJIvtk.addQuantity(Jimag, "Jimag");

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

        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();

    }

}
