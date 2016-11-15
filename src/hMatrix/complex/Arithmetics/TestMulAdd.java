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
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.Graphics.printHmatrix;
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
public class TestMulAdd {

    public static void main(String[] args) throws IOException {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(2);
        System.out.println("TestMulAdd.java");
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/hMatrix/mesh/";

        deb = System.nanoTime();
        /*
         ImportGmshMesh mesh1 = new ImportGmshMesh(meshDir + "plaque/plaque_4548.msh");
         mesh1.meshSummary();
         ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
         /*/
        ImportFlux mesh1 = new ImportFlux(meshDir + "sphere/SPHERE_1884.DEC");
//        ImportFlux mesh1 = new ImportFlux(meshDir + "HCAFARMESH.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        //*/
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
        int kmax = 50, nmin = 50;
        double eps = 1e-4, eta = 2.0;
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
        hc.printOnJFrame("Complex");
        System.out.println("index = " + Arrays.toString(hc.getDof2Idx()));

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ColumnVector rhsv, r_exact, x_exact, r_hm, x_hm;
        double vr[] = new double[2 * C.getActiveDofCount()];
        Arrays.fill(vr, 1.0);
        rhsv = new ColumnVector(vr);

        r_exact = new ColumnVector(M.product(vr, new double[M.getRows() * 2]));
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur assemblage = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        StorageHmatrixComplex hc2 = hc.copy(true);
        hc2.mulAdd(hc, hc, new TruncationControlComplex("rel", 1e-4));
        hc2.printOnJFrame("MulAdded");

        M = Basic2D.add(M.copy(), Basic2D.multiply(M.copy(), M.copy()));

        r_exact = new ColumnVector(M.product(vr, new double[M.getRows() * 2]));
        r_hm = hc2.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur muladd = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void mainSIMPLE(String[] args) throws IOException {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(2);
        System.out.println("TestComplexAssembly.java");
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/";

        deb = System.nanoTime();
        /*
         ImportGmshMesh mesh1 = new ImportGmshMesh(meshDir + "plaque/plaque_4548.msh");
         mesh1.meshSummary();
         ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
         /*/
//        ImportFlux mesh1 = new ImportFlux(meshDir + "sphere/SPHERE_104.DEC");
        ImportFlux mesh1 = new ImportFlux(meshDir + "HCAFARMESH.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        //*/
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
        hc.printOnJFrame();
        System.out.println("index = " + Arrays.toString(hc.getDof2Idx()));

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ColumnVector rhsv, r_exact, x_exact, r_hm, x_hm;
        double vr[] = new double[2 * C.getActiveDofCount()];
        Arrays.fill(vr, 1.0);
        rhsv = new ColumnVector(vr);

        r_exact = new ColumnVector(M.product(vr, new double[M.getRows() * 2]));
        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur assemblage = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        StorageHmatrixComplex hc2 = hc.copy(true);
        hc2.mulAdd(hc, hc, new TruncationControlComplex("rel", 1e-4));
        hc2.printOnJFrame();

        M = Basic2D.add(M.copy(), Basic2D.multiply(M.copy(), M.copy()));

        r_exact = new ColumnVector(M.product(vr, new double[M.getRows() * 2]));
        r_hm = hc2.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur muladd = " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        checkError(hc2, M.getArray());

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        System.out.println();
        System.out.println();
        System.out.println();
        //
        //
        // On va tester le bloc en haut a droite.

        //
        BlockClusterComplex bc = hc.getRoot().getSon(1);
        RkMatrixComplex rk = bc.getRM();

        // 1ere etape !
        rk.addProd(hc.getRoot().getSon(0), bc, new TruncationControlComplex("rel", 1e-4));
//        mat = M.getArray();
        Basic2D m1 = getSubMatrix(bc, hc, mat);
        Basic2D m0 = getSubMatrix(hc.getRoot().getSon(0), hc, mat);
        Basic2D map = Basic2D.add(m1, Basic2D.multiply(m0, m1));

        Basic2D s = Basic2D.sub(rk.getFullMatrix(), map);
        System.out.println("Erreur add prod 0 = " + s.norm2() / map.norm2());

        // 2eme etape !
        rk.addProd(bc, hc.getRoot().getSon(3), new TruncationControlComplex("rel", 1e-4));
        m0 = getSubMatrix(bc, hc, mat);
        m1 = getSubMatrix(hc.getRoot().getSon(3), hc, mat);
        map = Basic2D.add(map, Basic2D.multiply(m0, m1));

        s = Basic2D.sub(rk.getFullMatrix(), map);
        System.out.println("Erreur add prod 1 = " + s.norm2() / map.norm2());


        printHmatrix.closeAll();
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    protected static Basic2D getSubMatrix(BlockClusterComplex bc, StorageHmatrixComplex hc, double mat[][]) {
        int index[] = hc.getDof2Idx();
        double mm[][] = new double[bc.getRowSize()][2 * bc.getColSize()];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < bc.getColSize(); j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
        return new Basic2D(mm);
    }

    protected static void checkError(StorageHmatrixComplex hc, double mat[][]) {
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        BlockClusterComplex bc = hc.getRoot().getSon(0);
//        System.out.println(bc.getRM());
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
        System.out.println("Erreur sub bloc 0 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(0).getSon(0);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));
        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
        System.out.println("Erreur sub bloc 0 -> 0 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(0).getSon(1);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));

        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
//        System.out.println("s=" + s);
        System.out.println("Erreur sub bloc 0 -> 1 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(0).getSon(2);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));

        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
        System.out.println("Erreur sub bloc 0 -> 2 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(0).getSon(3);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));

        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
        System.out.println("Erreur sub bloc 0 -> 3 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(1);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));

        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
        System.out.println("Erreur sub bloc 1 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(2);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));

        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
        System.out.println("Erreur sub bloc 2 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(3);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));

        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
        System.out.println("Erreur sub bloc 3 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(3).getSon(0);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));
        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
        System.out.println("Erreur sub bloc 3 -> 0 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(3).getSon(1);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));

        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
//        System.out.println("s=" + s);
        System.out.println("Erreur sub bloc 3 -> 1 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(3).getSon(2);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));

        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
        System.out.println("Erreur sub bloc 3 -> 2 = " + s.norm2() / mmm.norm2());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        // cas hca far mesh
        bc = hc.getRoot().getSon(3).getSon(3);
//        System.out.println(bc.getRM());
        subm = bc.getM();
        sm = subm.getArray();
        index = hc.getDof2Idx();
        mm = new double[sm.length][sm[0].length];
        for (int i = 0; i < mm.length; i++) {
            for (int j = 0; j < sm[0].length / 2; j++) {
                mm[i][2 * j] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()]];
                mm[i][2 * j + 1] = mat[index[i + bc.getRowStart()]][2 * index[j + bc.getColStart()] + 1];
            }
        }
//        System.out.println("M=" + Basic2D.toString(mm));
//        System.out.println("Aprox=" + Basic2D.toString(sm));

        mmm = new Basic2D(mm);
        s = Basic2D.sub(subm, mmm);
        System.out.println("Erreur sub bloc 3 -> 3 = " + s.norm2() / mmm.norm2());

    }

}
