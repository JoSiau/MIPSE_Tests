/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.complex.Arithmetics;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.elements.Hexaedre;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.KernelInterface;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.ComplexOperation;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.StorageBlock.Admissible.RkMatrixComplex;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.TruncationControlComplex;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.binaryTree.BlockClusterComplex;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.decompositions.matlab.DecompositionMatLab;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulation;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulationHmatrix;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.mipseCore.storage.StorageHmatrixComplex;
import g2elab.mipse.mipseCore.storage.StorageSparse;
import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D;
import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author jsiau
 */
public class TestRkMatrixOperations {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(2);
        System.out.println("TestComplexAssembly.java");
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/";

        deb = System.nanoTime();
        ImportFlux mesh1 = new ImportFlux(meshDir + "HCAFARMESH.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        fin = System.nanoTime();
        System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);

        int d = 3;
        int n = ES.getNbElement();
        System.out.println(" n = " + n + "\t d = " + d);
        /*
         FaceDeg1 C = new FaceDeg1(ES);
         KernelInterface k = new MultGvect();
         /*/
        Cell C = new Cell(ES);
        KernelInterface k = new MultG();
        //*/
        int nbPGSrc, nbPGTrg;
        if (ES.getElements(0) instanceof Hexaedre) {
            nbPGSrc = 8;
            nbPGTrg = 27;
        } else {
            nbPGSrc = 4;
            nbPGTrg = 7;
        }

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        FiniteElementFormulation IF = new FiniteElementFormulation(C, C);
        IF.assembly(nbPGSrc);
        SparseMatrixRowReal S = ((StorageSparse) IF.getStore()).getSparseMat();

        GalerkinIntegralFormulationFull full = new GalerkinIntegralFormulationFull(C, C, k, new SelfElementFixedGauss(nbPGSrc, new InCreasedPGSourceNumber(nbPGTrg)), nbPGSrc);
        full.assembly();
        Matrix m = ((StorageFull) full.getStore()).getMatrix();

        double mat[][] = new double[C.getActiveDofCount()][C.getActiveDofCount() * 2];
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < C.getActiveDofCount(); j++) {
                /*
                 mat[i][2 * j] = m.getElement(i, j);
                 mat[i][2 * j + 1] = S.getQuick(i, j);
                 /*/
                mat[i][2 * j] = S.getQuick(i, j);
                mat[i][2 * j + 1] = m.getElement(i, j);
                //*/
            }
        }
        Basic2D M = new Basic2D(mat);
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        int kmax = 50, nmin = 4;
        double eps = 1e-4, eta = 20.0;
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, C, k, new SelfElementFixedGauss(nbPGSrc, new InCreasedPGSourceNumber(nbPGTrg)), nbPGSrc, nbPGSrc,
                eps, kmax, nmin, 5, eta);
        //*
        f.assembly();
        StorageHmatrix hIV = (StorageHmatrix) f.getStore();
        hIV.CheckError();
        //
        FiniteElementFormulationHmatrix Ifh = new FiniteElementFormulationHmatrix(C, nmin, eta);
        Ifh.assembly(1);
        StorageHmatrix hFE = (StorageHmatrix) Ifh.getStore();
        //
        // Complex assembly                
        StorageHmatrixComplex hc = new StorageHmatrixComplex(hFE, hIV, null);
//        hc.printOnJFrame();
        System.out.println("index = " + Arrays.toString(hc.getDof2Idx()));

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ColumnVector rhsv, r_exact, r_hm;
        double vr[] = new double[2 * C.getActiveDofCount()];
        Arrays.fill(vr, 1.0);
        rhsv = new ColumnVector(vr);

        r_exact = new ColumnVector(M.product(vr, new double[M.getRows() * 2]));
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur assemblage = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        ////////////////////////////////////////////////////////////////////////
        // EXTRACTION !
        ////////////////////////////////////////////////////////////////////////
        BlockClusterComplex bc = hc.getRoot().getSon(1);
        RkMatrixComplex rk = bc.getRM();
        Basic2D subm = rk.getFullMatrix();
        int index[] = hc.getDof2Idx();
        double mm[][] = new double[bc.getRowSize()][2 * bc.getColSize()];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < bc.getColSize(); j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
        Basic2D fm = new Basic2D(mm);
        Basic2D s = Basic2D.sub(subm, fm);
        System.out.println("Erreur bloc = " + s.norm2() / fm.norm2());

        System.out.println("*******************************");
        System.out.println("Inner-product form operations :");
        System.out.println("*******************************");
        ////////////////////////////////////////////////////////////////////////
        // ADDITION !
        ////////////////////////////////////////////////////////////////////////
        Basic2D add = Basic2D.add(fm, fm);

        RkMatrixComplex radd = rk.RoundedAdd(rk, new TruncationControlComplex("rel", 1e-4));

        s = Basic2D.sub(radd.getFullMatrix(), add);
        System.out.println("Erreur add UV = " + s.norm2() / add.norm2());

        ////////////////////////////////////////////////////////////////////////
        // PRODUIT !
        ////////////////////////////////////////////////////////////////////////
        Basic2D mul = Basic2D.multiply(fm, fm);

        RkMatrixComplex rmul = rk.prod(rk);

        s = Basic2D.sub(rmul.getFullMatrix(), mul);
        System.out.println("Erreur mul rk-rk = " + s.norm2() / mul.norm2());

        ////////////////////////////////////////////////////////////////////////
        Basic2D mul2 = rk.prod(fm);

        s = Basic2D.sub(mul2.getFullMatrix(), mul);
        System.out.println("Erreur mul rk-full = " + s.norm2() / mul.norm2());

        ////////////////////////////////////////////////////////////////////////
        // orthonormal outer-product form
        ////////////////////////////////////////////////////////////////////////
        System.out.println("*******************************");
        System.out.println("Outer-product form operations :");
        System.out.println("*******************************");

        rk.SVD(new TruncationControlComplex("rel", 1e-4));
        radd = rk.RoundedAdd(rk, new TruncationControlComplex("rel", 1e-4));
        s = Basic2D.sub(radd.getFullMatrix(), add);
        System.out.println("Erreur add USV = " + s.norm2() / add.norm2());

        ////////////////////////////////////////////////////////////////////////
        // PRODUIT !
        ////////////////////////////////////////////////////////////////////////
        RkMatrixComplex rmulb = rk.prod(rk);

        s = Basic2D.sub(rmulb.getFullMatrix(), mul);
        System.out.println("Erreur mul rk-rk = " + s.norm2() / mul.norm2());

        ////////////////////////////////////////////////////////////////////////
        Basic2D mul2b = rk.prod(fm);

        s = Basic2D.sub(mul2b.getFullMatrix(), mul);
        System.out.println("Erreur mul rk-full = " + s.norm2() / mul.norm2());

        DecompositionMatLab.exit();
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
