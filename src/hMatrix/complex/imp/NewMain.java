/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.complex.imp;

import static formulations.formulationT.cefc2016.Test_TIm.conduc;
import static formulations.formulationT.cefc2016.Test_TIm.ep;
import static formulations.formulationT.cefc2016.Test_TIm.f;
import static formulations.formulationT.cefc2016.Test_TIm.mesh;
import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.FormulationT.FormulationTJ_alpha;
import g2elab.mipse.material.AbstractMaterial;
import g2elab.mipse.material.conductor.LinearConductor;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementFactory;
import g2elab.mipse.meshCore.functionSpace.FunctionSpace;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.mipseCore.blockAssembly.complex.MatrixComplexPartReImConvert;
import g2elab.mipse.mipseCore.blockAssembly.real.RealMatrixBlocFmult;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.BinaryTree.BlockClusterTree;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.mipseCore.storage.StorageHmatrixComplex;
import g2elab.mipse.mipseCore.storage.StorageSparse;
import g2elab.mipse.numericalTools.matrix.complex.MatrixComplexPartReIm;
import g2elab.mipse.numericalTools.matrix.complex.bloc.ComplexMatrixAssembly;
import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author jsiau
 */
public class NewMain {

    /**
     * TEST with a simple H-matrix
     *
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        ElementFactory ef = new ElementFactory(ElementFactory.STEP, ElementFactory.STEP, ElementFactory.STEP);
        System.out.println("Test avec trou !");
        File file = new File("");
        String path;
        path = file.getAbsolutePath();
        path = path + mesh;
        System.out.println("Fichier maillage : " + path);
        ImportFlux IF = new ImportFlux(path);
        RegionsSet regSet = new RegionsSet(new Region[]{IF.getRegions()[0]});
        ((SurfaceRegion) regSet.getRegion(0)).setThickness(ep);

        // Matériaux
        AbstractMaterial[] materials = new AbstractMaterial[1];
        materials[0] = new LinearConductor(conduc);
        //
        //
        ////////////////////////////////////////////////////////////////////////
        FormulationTJ_alpha sol = new FormulationTJ_alpha();
        sol.DEBUGG_EXPORTMATRICES = false;
        sol.DEBUGG = false;
        sol.setPrecond(1);
        //
        //
        sol.setMeshGroup(regSet);
        sol.setMaterials(materials);
        sol.setCapacitiveComputation(false);
        sol.setFrequency(f);
        BlocElectriqueBasique circuit = new BlocElectriqueBasique();
        // tor
        circuit.addSourceISimple(1000000, 960, 1679, "source", 1.0, 0.0);
        circuit.finSaisie();
        sol.setElectricalCircuit(circuit);
        sol.setMeshStuff();
        sol.setCompression(Compression.HCA);
        StorageHmatrixComplex h = extractHmatComplex(new MatrixComplexPartReImConvert((MatrixComplexPartReIm) sol.integration().getBlock(0, 0)),sol.getRotAlpha());
        
        sol.setCompression(Compression.No);
        ComplexMatrixAssembly m = sol.integration();
        

    }

    public static StorageHmatrixComplex extractHmatComplex(MatrixComplexPartReImConvert riz, FunctionSpace rotAlpha) {
        // The real part is always sparse
        SparseMatrixRowReal r = ((StorageSparse) riz.getRealPart()).getSparseMat();
        //
        BlockClusterTree bct = new BlockClusterTree(rotAlpha, r, 3, 30, 2.0);
        StorageHmatrixComplex h;
        StorageHmatrix hreal;
        if (riz.getImagPart() instanceof StorageHmatrix) {
            hreal = (StorageHmatrix) riz.getImagPart();
            h = new StorageHmatrixComplex(new StorageHmatrix(bct), hreal, null);
        } else if (riz.getImagPart() instanceof RealMatrixBlocFmult) {
            RealMatrixBlocFmult tmp = (RealMatrixBlocFmult) riz.getImagPart();
            hreal = (StorageHmatrix) tmp.getMat();
            h = new StorageHmatrixComplex(1.0, new StorageHmatrix(bct), tmp.getScaleFactor(), hreal, null);
        } else {
            throw new UnsupportedOperationException();
        }
        return h;
    }

}
