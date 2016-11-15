/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.formulationT.loopFace;

import g2elab.mipse.analytical.sourceFields.SourceField;
import g2elab.mipse.analytical.sourceFields.UniformMagneticField;
import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.LoopFace.FormulationLoopFaceBrutBeta;
import g2elab.mipse.formulationInProgress.magnetodynamic.Pertes;
import g2elab.mipse.material.AbstractMaterial;
import g2elab.mipse.material.conductor.LinearConductor;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.Quantity;
import g2elab.mipse.meshCore.quantity.RealVectorCellQuantity;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.numericalTools.iterativeSolver.Solver;
import g2elab.mipse.tools.debugTools.VectorCheckTools;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.files.Lecture;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

/**
 *
 * @author jsiau
 */
public class TestScaling {

    static FaceDeg1 faces;

    /**
     * Test uniquement le couplage circuit
     *
     * @param args
     */
    public static void main(String[] args) throws IOException {
        Matrix res, jF, res1, jF1;
        int nbFace;
        String saveingPath = "d:/tmp/";
        //
        Lecture fic;
        try {
            fic = new Lecture(saveingPath + "res.out");
            System.out.println("Reading res.out");
            int n = fic.lireLigneInt();
            nbFace = fic.lireLigneInt();
            double[] t = fic.lireLigneDoubleVecteur(' ');
            res = new Matrix(n, 2, t);
            fic.close();
        } catch (IOException ex) {
            System.out.println("Computing res.out");
            res = solve(false, null).copy();
            Ecriture ec = new Ecriture(saveingPath + "res.out");
            ec.ecrire(res.getRowCount() + "\n");
            ec.ecrire(faces.getActiveDofCount() + "\n");
            double t[] = new double[2 * res.getRowCount()];
            res.get(t);
            ec.ecrire(t, ' ');
            ec.close();
            nbFace = faces.getActiveDofCount();
        }
        System.out.println("Nb Face = " + nbFace);
        jF = res.submatrix(res.getRowCount() - nbFace, 0, nbFace, 2);
        //
        //
        //
        try {
            fic = new Lecture(saveingPath + "resNorm.out");
            System.out.println("Reading resNorm.out");
            int n = fic.lireLigneInt();
            nbFace = fic.lireLigneInt();
            double[] t = fic.lireLigneDoubleVecteur(' ');
            res1 = new Matrix(n, 2, t);
            fic.close();
        } catch (IOException ex) {
            System.out.println("Computing resNorm.out");
            res1 = solve(true, null).copy();
            Ecriture ec = new Ecriture(saveingPath + "resNorm.out");
            ec.ecrire(res.getRowCount() + "\n");
            ec.ecrire(faces.getActiveDofCount() + "\n");
            double t[] = new double[2 * res1.getRowCount()];
            res1.get(t);
            ec.ecrire(t, ' ');
            ec.close();
            nbFace = faces.getActiveDofCount();
        }
        System.out.println("Nb Face = " + nbFace);
        jF1 = res1.submatrix(res.getRowCount() - nbFace, 0, nbFace, 2);

        //
        //
        //
        System.out.println("Error global = " + getError(res, res1));
        //
        //
        //
        System.out.println("Error faces = " + getError(jF, jF1));

        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

    private static Matrix solve(boolean normalize, ImportFlux IF) {
        if (IF == null) {
            File file = new File("");
            String path;
            path = file.getAbsolutePath();
//        path = path + "\\src\\g2elab\\mipse\\formulationInProgress\\magnetodynamic\\U_PEEC_DIELECTRIC\\C_MED.DEC";
            path += "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Surf/S3_NB_FS.DEC";
            IF = new ImportFlux(path);
        }

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

//        sol.plotCircuit();
//
        //
//        /*     
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
        //
        //
        circuit.addSourceISimple(1000000, 404, 809, "Source", 1.0, 0.0);
        circuit.addSourceUSimple(1000001, 297, 702, "cct", 0, 0);
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        double f = 1e7;
        /*/
         UniformMagneticField srcField = new UniformMagneticField(0.0, 0.0, 1.0); //1.0);
         sol.setSources(new SourceField[]{srcField});
         double f = 1e3;
         //*/
        sol.setFrequency(f);
        sol.setNormalization(normalize);
        //
        //
        //
        sol.setSolver(Solver.LU);
        ArrayList<Quantity> sortie = sol.solve();
        faces = sol.getFaces();
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
        FormulationLoopFaceBrutBeta.exportRes(sortie, "d:/tmp/C/" + (normalize ? "norm" : "") + "_" + f);
        return sol.getResMap().get(f);
    }

    private static double getError(Matrix ref, Matrix v) {
        double[] t = new double[ref.getRowCount() * ref.getColumnCount()];
        ref.get(t);
        double[] t1 = new double[ref.getRowCount() * ref.getColumnCount()];
        v.get(t1);
        return VectorCheckTools.getError(t1, t, false);
    }
}
