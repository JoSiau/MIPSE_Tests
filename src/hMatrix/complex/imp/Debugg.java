/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.complex.imp;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.contraints.ExternalDriveFaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceRealDirichlet;
import g2elab.mipse.meshCore.elements.ElementLinSetHomogene;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.meshCore.region.LineRegion;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.IntegrationCorrectionStrategySource;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.BinaryTree.BlockClusterTree;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.Preconditionners.HmatrixLUDecompositionComplex;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.TruncationControlComplex;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.mipseCore.storage.StorageHmatrixComplex;
import g2elab.mipse.mipseCore.storage.StorageSparse;
import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
import g2elab.mipse.numericalTools.vector.full.VectorFullComplex;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.RowVector;

import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author jsiau
 */
public class Debugg {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        double f = 1e6;
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(Debugg.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir = meshDir + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/LOOPANTENNA_SURF3.DEC";
        ImportFlux mesh = new ImportFlux(meshDir);
        mesh.getRegions();

        SurfaceRegion conductor = (SurfaceRegion) mesh.getRegion(0);
        LineRegion posFlux = (LineRegion) mesh.getRegion(1);
        LineRegion negFlux = (LineRegion) mesh.getRegion(2);

        LineRegion borderNull = (LineRegion) conductor.generateBorder();

        borderNull = new LineRegion((ElementLinSetHomogene) borderNull.getComplement((LineRegion) posFlux, 1));
        borderNull = new LineRegion((ElementLinSetHomogene) borderNull.getComplement((LineRegion) negFlux, 1));
        // When the conductor is a surface, we have to descativavte the border.
        FaceDeg1 FDcond = new FaceDeg1((ElementSetHomogene) conductor.getElementSet(),
                new FaceConstraint[]{new FaceRealDirichlet(borderNull, 0.0),
                        new ExternalDriveFaceConstraint(posFlux), new ExternalDriveFaceConstraint(negFlux)});

        double sigma = 1 / 1.68e-8, ep = 20e-6;

        System.out.println("===== Integration: [R] =====");
        long t = System.currentTimeMillis();
        RowVector v = new RowVector(FDcond.getNbElement());
        v.setAllElements(1 / (sigma * ep));
        // Create the quanty 1/sigma 
        Cell support = new Cell(FDcond.getElementSet());
        RealScalarCellQuantity q = new RealScalarCellQuantity(support, v);
        // Assemblage de [R]
        FiniteElementFormulation EF = new FiniteElementFormulation(FDcond);
        EF.assembly(9, q);
        SparseMatrixRowReal R = ((StorageSparse) EF.getStore()).getSparseMat();
        System.out.println("Assembly Time: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");
        //
        // L 
        //
        System.out.println("===== Integration: [L] =====");
        IntegrationCorrectionStrategySource strategy;
        // Choose the integration strategy.
        System.out.println("Strategie d'integration: ");
        strategy = new SelfElementFixedGauss(9, new InCreasedPGSourceNumber(25));
        v.setAllElements(4.0 * Math.PI * 1E-7);
        // Create the quanty
        q = new RealScalarCellQuantity(support, v);
        // Assembly the matrix [L] as an H-matrix or Full-matrix
        GalerkinIntegralFormulation IV;
        t = System.currentTimeMillis();

        double epsHmat = 1e-4, eta = 2.0;
        int nmin = 30, kmax = 50, order = 3;
        boolean recomp = true;

        IV = new GalerkinIntegralFormulationHCA(FDcond, FDcond, new MultGvect(), strategy, 9, 9,
                epsHmat, kmax, nmin, order, eta, recomp);
        IV.assembly(q);
        StorageHmatrix Lh = (StorageHmatrix) IV.getStore();
        System.out.println("Assembly Time: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");

        BlockClusterTree bct = new BlockClusterTree(FDcond, R, 3, nmin, eta);
        StorageHmatrix hR = new StorageHmatrix(bct);
        double omega = 2 * Math.PI * f;
        //
        // hc = R + j w L
        StorageHmatrixComplex hc = new StorageHmatrixComplex(1.0, hR, omega, Lh, null);
        StorageHmatrixComplex h2 = hc.copy(true);
        TruncationControlComplex tol = new TruncationControlComplex("rel", 3e-1);
        hc.Agglomerate(tol);
        /////////////////////////// >TMP
        double b[] = new double[FDcond.getActiveDofCount() * 2];
        Arrays.fill(b, 1.12);
        //
        HmatrixLUDecompositionComplex hlu = new HmatrixLUDecompositionComplex(hc, tol);
        double x[] = hlu.precond(b, new double[b.length]);
        //
        double b2[] = h2.product(x, new double[b.length]);
        // Compute the error
        VectorFullComplex v1 = new VectorFullComplex(b);
        VectorFullComplex v2 = new VectorFullComplex(b2);
        System.out.println("b2 = " + Arrays.toString(b2));
        v2.sub(v1, v2);
        System.out.println("Erreur hlu = " + v2.norm2() / v1.norm2());
        System.out.println("");


        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
