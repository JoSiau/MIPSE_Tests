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
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.Quantity;
import g2elab.mipse.meshCore.quantity.RealVectorCellQuantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.tools.files.Ecriture;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import formulations.uPEEC.impedanceCurves.Zf_PEEC_RLMPC_SURF_SERPENT;

/**
 * The type Test loop antenna.
 * @author jsiau
 */
public class TestLoopAntenna {

    /**
     * The entry point of application.
     *
     * @param args the command line arguments
     * @throws IOException the iO exception
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
        FormulationTJ_alpha sol = new FormulationTJ_alpha();
        sol.setMeshGroup(new RegionsSet(new Region[]{cond}));
        sol.setMaterials(materials);
        sol.setCompression(Compression.No);
//        sol.setCompression("HCA");
        //
        //
//        sol.plotCircuit();
        //*
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();

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

        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        /*/
         sol.setSourceField(new  SourceField[]{new UniformField(0, 0, 1)});
         //*/
        //
        double[] f = Zf_PEEC_RLMPC_SURF_SERPENT.compteFrequenciesLog(9, 10, 100);
//        double[] f = new double[]{2.310129700083158E9};
//        sol.setFrequency(f);
        sol.setFrequency(new double[]{1e9});
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
        System.out.println("U_src I = ");
        Map<Double, Matrix> m = sol.getResMap();
        int ind = sol.getRotAlpha().getActiveDofCount() + sol.getNbMILoops();
        for (int i = 0; i < m.size(); i++) {
            Matrix mat = m.get(f[i]);
            if (mat != null) {
                System.out.println(f[i] + "   " + mat.getElement(ind, 0) + " +j* " + mat.getElement(ind, 1));
            } else {
                System.out.println(f[i] + "   -1");
            }
        }

        System.out.println("Nb Iterations :");
        for (int i = 0; i < f.length; i++) {
            System.out.println(f[i] + "   " + sol.getNbIt()[i]);
        }

        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();

    }

}
