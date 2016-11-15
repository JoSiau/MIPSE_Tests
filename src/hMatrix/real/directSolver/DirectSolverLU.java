/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.directSolver;

import g2elab.mipse.analytical.sourceFields.UniformElectricField;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmsh;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.gaussValues.VectorGaussPointsValues;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.Preconditionners.HmatrixLUDecomposition;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.LuDecomposition;
import got.matrix.Matrix;
import got.matrix.RowVector;

import java.io.IOException;

/**
 * TEST LA RESOLUTION DIRECTE DE LA DECOMPOSITION H-LU ET COMPARE LES RESULTATS
 * AVEC LE STOCKAGE PLEIN
 *
 * @author jsiau
 */
public class DirectSolverLU {

    public static final String[] sizes = {"SPHERE_572", "SPHERE_2104", "SPHERE_8262", "22k", "44k", "82k", "132k", "165k", "200k"};

    /**
     * @param args qq chose
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        // Load the mesh
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/hMatrix/mesh/sphere/";

        for (int i = 0; i < sizes.length; i++) {

            ImportFlux mesh = new ImportFlux(meshDir + sizes[i] + ".DEC");
            ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();

            int d = 3;
            int n = ES.getNbElement();
            int nbGauss = 7;
            /*
             SECOND MEMBRE
             */
            UniformElectricField e0 = new UniformElectricField(100, 0, 0);
            VectorGaussPointsValues e0gpv = e0.getFieldGPValue(ES, 1);
            Matrix t = e0gpv.projWithExplicitDof(new Cell(ES));

            ColumnVector secondMembre = ((RowVector) t).transpose();
//        ColumnVector secondMembre = new ColumnVector(n);
//        for (int i = 0; i  <n; i++) {
//            secondMembre.setElement(i, Math.random() * 100);
//        }
        /*
             HMATRIX
             */
            long deb = System.currentTimeMillis();
            double eps = 1e-5;
            int order = 4;
            TruncationControl tol = new TruncationControl("rel", eps);
            int kmax = 50, nmin = 30;
            Cell C = new Cell(ES);
            /*
             GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss,
             eps, kmax, nmin);
             /*/
            GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss, nbGauss,
                    eps, kmax, nmin, order);
            //*/
            f.assembly();
            StorageHmatrix H0 = (StorageHmatrix) f.getStore();
            H0.Agglomerate(tol);

            ColumnVector resH = new ColumnVector(n);
            long t0 = System.currentTimeMillis();
            HmatrixLUDecomposition hLu = new HmatrixLUDecomposition(H0, tol);
            long t1 = System.currentTimeMillis();
            System.out.println("Time dec. = " + (t1 - t0) / 1e3 + " sec.");
            //////////////////////////////
            t0 = System.currentTimeMillis();
            hLu.solve(resH, secondMembre);
            t1 = System.currentTimeMillis();
            System.out.println("Time solve. = " + (t1 - t0) / 1e3 + " sec.");
            ///////////////////////////////
            long fin = System.currentTimeMillis();
            System.out.println("Total time Hmatrix storage = " + (fin - deb) * 1e-3);
//        System.out.println("resH: \n"+resH.transpose());

            RealScalarCellQuantity qH = new RealScalarCellQuantity(C, resH.transpose());
            qH.exportGmsh(new ExportGmsh(ES, "d:/tmp/qH.msh"), "qH");
        }

//        /*
//         FULL STORAGE
//         */
//        deb = System.currentTimeMillis();
//        GalerkinIntegralFormulation IFCell = new GalerkinIntegralFormulationFull(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss);
//        IFCell.assembly();
//        Matrix Ai = ((StorageFull) IFCell.getStore()).getMatrix();
//        LuDecomposition MD = new LuDecomposition(Ai);
//        ColumnVector resEx = new ColumnVector(n);
//        MD.solve(secondMembre, resEx);
//        fin = System.currentTimeMillis();
//        System.out.println("Total time Full Storage = " + (fin - deb) * 1e-3);
////        System.out.println("resEx: \n"+resEx.transpose());
//        RealScalarCellQuantity q = new RealScalarCellQuantity(C, resEx.transpose());
//        q.exportGmsh(new ExportGmsh(ES, "d:/tmp/q.msh"), "q");
//
//        /*
//         CHECK ERROR
//         */
//        resH.sub(resEx);
//        System.err.println("Error relative = " + resH.norm() / resEx.norm());
    }
}
