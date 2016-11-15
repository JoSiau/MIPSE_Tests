/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Arithmetics;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.BinaryTree.BlockClusterTree;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.StorageBlock.Admissible.RkMatrix;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.clustering.Cluster;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import got.matrix.ColumnVector;
import got.matrix.Matrix;
import got.matrix.SVDDecomposition;

import java.io.IOException;

/**
 * Test l'addition des matrices de bas rang (RkMatrix)
 *
 * @author jsiau
 */
public class RoundedAddition {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {

        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
        ImportFlux mesh = new ImportFlux(meshDir + "/sphere/SPHERE_1884.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();

        int d = 3;
        int n = ES.getNbElement();
        int nbGauss = 7;

        /*
         HMATRIX
         */
        double eps = 1e-2;
        TruncationControl tol = new TruncationControl("rel", eps);
        int kmax = 50, nmin = 30;
        Cell C = new Cell(ES);
        GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss,
                eps, kmax, nmin);

        BlockClusterTree bct = new BlockClusterTree(f.getFunctionSpaceSource(), 3, nmin, 2.0);
        f.setIndex(bct.getDof2IdxTrg(), bct.getDof2IdxSrc());

        RkMatrix R1 = f.getCompressedBlock(new Cluster(3 * n / 4, n / 4, 0, 3), new Cluster(0, n / 4, 0, 3));
        System.out.println("R1 rank= " + R1.getRank());

        Matrix M1 = f.getFullBlock(new Cluster(3 * n / 4, n / 4, 0, 3), new Cluster(0, n / 4, 0, 3)).getFullMatrix();
        double Mnorm = M1.norm();

        Matrix tmp = R1.getFullMatrix();
        tmp.sub(M1);
        System.out.println("Erreur relative de départ= " + tmp.norm() / Mnorm);
        M1.scale(2);
        RkMatrix R2 = R1.RoundedAdd(R1, tol);
        System.out.println("R2 rank= " + R2.getRank());
        tmp = R2.getFullMatrix();
        tmp.sub(M1);
        System.out.println("Erreur relative apres addition= " + tmp.norm() / (M1.norm()));

        ////////////////////////////////////////////////////////////////////////
        Matrix M2 = f.getFullBlock(new Cluster(n / 4, n / 4, 0, 3), new Cluster(0, n / 4, 0, 3)).getFullMatrix();
        long deb = System.nanoTime();
        SVDDecomposition Svd = new SVDDecomposition(M2);
        ColumnVector Sv = new ColumnVector(Svd.getColumnCount());
        Svd.getS(Sv);
        Matrix Uv = new Matrix(Svd.getRowCount(), Svd.getColumnCount());
        Svd.getU(Uv);
        Matrix Vv = new Matrix(Svd.getColumnCount(), Svd.getColumnCount());
        Svd.getV(Vv);

        RkMatrix R3 = RkMatrix.Truncation(Uv, Sv, Vv, tol);
        System.out.println("R3 rank= " + R3.getRank());

        R2 = R3.RoundedAdd(R1, tol);
        System.out.println("R1 + R3 rank= " + R2.getRank());
        long fin = System.nanoTime();
        System.out.println("T1= " + (fin - deb) * 1e-9);
        ////////////////////////////////////////////////////////////////////////

        deb = System.nanoTime();
        M2.add(R1.getFullMatrix());
        Svd.decompose(M2);
        Sv = new ColumnVector(Svd.getColumnCount());
        Svd.getS(Sv);
        Uv = new Matrix(Svd.getRowCount(), Svd.getColumnCount());
        Svd.getU(Uv);
        Vv = new Matrix(Svd.getColumnCount(), Svd.getColumnCount());
        Svd.getV(Vv);

        RkMatrix R4 = RkMatrix.Truncation(Uv, Sv, Vv, tol);
        System.out.println("R4 rank= " + R4.getRank());

        fin = System.nanoTime();
        System.out.println("T2= " + (fin - deb) * 1e-9);
    }

}
