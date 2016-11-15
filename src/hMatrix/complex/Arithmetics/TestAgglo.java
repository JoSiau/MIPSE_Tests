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
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.functionSpace.FunctionSpace;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.KernelInterface;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.ComplexOperation;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.TruncationControlComplex;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.binaryTree.BlockClusterComplex;
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
public class TestAgglo {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(2);
        System.out.println("TestAgglo.java");
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/hMatrix/mesh/";

        deb = System.nanoTime();
        /*
         ImportGmshMesh mesh1 = new ImportGmshMesh(meshDir + "plaque/plaque_4548.msh");
         mesh1.meshSummary();
         ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
         /*/
        ImportFlux mesh1 = new ImportFlux(meshDir + "sphere/SPHERE_572.DEC");
//        ImportFlux mesh1 = new ImportFlux(meshDir + "HCAFARMESH.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        //*/
        fin = System.nanoTime();
        System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);

        int d = 3;
        int n = ES.getNbElement();
        System.out.println(" n = " + n + "\t d = " + d);
//        /*
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
        int kmax = 50, nmin = 8;
        double eps = 1e-4, eta = 2.0;
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, C, k, new SelfElementFixedGauss(nbPGSrc, new InCreasedPGSourceNumber(nbPGTrg)), nbPGSrc, nbPGSrc,
                eps, kmax, nmin, 5, eta);
        //*
        f.assembly();
        StorageHmatrix hIV = (StorageHmatrix) f.getStore();
        hIV.CheckError();
        //
        FiniteElementFormulationHmatrix Ifh = new FiniteElementFormulationHmatrix(C, nmin, eta);
        Ifh.assembly(nbPGSrc);
        StorageHmatrix hFE = (StorageHmatrix) Ifh.getStore();
        //
        // Complex assembly                
        StorageHmatrixComplex hc = new StorageHmatrixComplex(hFE, hIV, null);
        hc.printOnJFrame();

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        /*
         System.out.println("Avant agglo : ");
         test0(hc, mat);
        
         hc.Agglomerate(new TruncationControlComplex("rel", 1e-2));
         hc.printOnJFrame("Agglo 1e-2");
         test0(hc, mat);
         /*/
//        testSingle(hc, M, C);
        testAll(hc, M, C);
        //*/
    }

    public static void testSingle(StorageHmatrixComplex hc, Basic2D M, FunctionSpace C) {
        ColumnVector rhsv, r_exact, x_exact, r_hm, x_hm;
        double vr[] = new double[2 * C.getActiveDofCount()];
        Arrays.fill(vr, 1.0);
        rhsv = new ColumnVector(vr);
        // Vecteur de reference
        r_exact = new ColumnVector(M.product(vr, new double[M.getRows() * 2]));
        //
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur assemblage = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        hc.Agglomerate(new TruncationControlComplex("rel", 1e-4));
        hc.printOnJFrame("1e-4");
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur agglo 1e-4 = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));  
    }

    public static void test0(StorageHmatrixComplex hc, double mat[][]) {
        BlockClusterComplex bc = hc.getRoot().getSon(1);
        Basic2D subm = bc.getM();
        double sm[][] = subm.getArray();
        int index[] = hc.getDof2Idx();
        double mm[][] = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));
        Basic2D mmm = new Basic2D(mm);
        Basic2D s = Basic2D.sub(subm, mmm);
        System.out.println("Erreur = " + s.norm2() / mmm.norm2());
    }

    public static void testAll(StorageHmatrixComplex hc, Basic2D M, FunctionSpace C) {
        ColumnVector rhsv, r_exact, x_exact, r_hm, x_hm;
        double vr[] = new double[2 * C.getActiveDofCount()];
        Arrays.fill(vr, 1.0);
        rhsv = new ColumnVector(vr);
        // Vecteur de reference
        r_exact = new ColumnVector(M.product(vr, new double[M.getRows() * 2]));
        //
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur assemblage = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        hc.Agglomerate(new TruncationControlComplex("rel", 1e-4));
        hc.printOnJFrame("1e-4");
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur agglo 1e-4 = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        hc.Agglomerate(new TruncationControlComplex("rel", 1e-3));
        hc.printOnJFrame("1e-3");
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur agglo 1e-3 = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        hc.Agglomerate(new TruncationControlComplex("rel", 1e-2));
        hc.printOnJFrame("1e-2");
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur agglo 1e-2 = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        hc.Agglomerate(new TruncationControlComplex("rel", 1e-1));
        hc.printOnJFrame("1e-1");
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur agglo 1e-1 = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        hc.Agglomerate(new TruncationControlComplex("rel", 3e-1));
        hc.printOnJFrame("3e-1");
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur agglo 3e-1 = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
