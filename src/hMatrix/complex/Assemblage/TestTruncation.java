/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.complex.Assemblage;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.elements.Hexaedre;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.KernelInterface;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.ComplexOperation;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.StorageBlock.Admissible.RkMatrixComplex;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.TruncationControlComplex;
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
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class TestTruncation {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(2);
        System.out.println("testTruncation.java");
        double deb, fin;
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(TestTruncation.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/hMatrix/mesh/";
        deb = System.nanoTime();
        ImportFlux mesh1 = new ImportFlux(meshDir + "HCAFARMESH.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        fin = System.nanoTime();
        System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);
//        /*
         FaceDeg1 C = new FaceDeg1(ES);
         KernelInterface k = new MultGvect();
         /*/
        Cell C = new Cell(ES);
        KernelInterface k = new MultG();
        //*/
        int nbPG = 4, nbPGp = 7;
        if (ES.getElements(0) instanceof Hexaedre) {
            nbPG = 8;
            nbPGp = 27;
        }   ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        FiniteElementFormulation IF = new FiniteElementFormulation(C, C);
        IF.assembly(nbPG);
        SparseMatrixRowReal S = ((StorageSparse) IF.getStore()).getSparseMat();
        GalerkinIntegralFormulationFull full = new GalerkinIntegralFormulationFull(C, C, k, new SelfElementFixedGauss(nbPG, new InCreasedPGSourceNumber(nbPGp)), nbPG);
        full.assembly();
        Matrix m = ((StorageFull) full.getStore()).getMatrix();
        double a = 0, b = 10, c = 10, d = 2;
        double mat[][] = new double[C.getActiveDofCount()][C.getActiveDofCount() * 2];
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < C.getActiveDofCount(); j++) {
                /*
                 mat[i][2 * j] = m.getElement(i, j);
                 mat[i][2 * j + 1] = S.getQuick(i, j);
                 /*/
                mat[i][2 * j] = a * S.getQuick(i, j) + c * m.getElement(i, j);
                mat[i][2 * j + 1] = b * S.getQuick(i, j) + d * m.getElement(i, j);
                //*/
            }
        }
        Basic2D M = new Basic2D(mat);
        int kmax = 50, nmin = 10;
        double eps = 1e-4, eta = 2.0;
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, C, k, new SelfElementFixedGauss(nbPG, new InCreasedPGSourceNumber(nbPGp)), nbPG, nbPG,
                eps, kmax, nmin, 5, eta);
        //*
        f.assembly();
        StorageHmatrix hIV = (StorageHmatrix) f.getStore();
        hIV.CheckError();
        hIV.printOnJFrame("IV Reel");
//        System.out.println(hIV.getRoot().getSon(1).getRM());
        //        System.out.println("S:\n" + S.getFullMatrix());
        FiniteElementFormulationHmatrix Ifh = new FiniteElementFormulationHmatrix(C, nmin, eta);
        Ifh.assembly(nbPG);
        StorageHmatrix hFE = (StorageHmatrix) Ifh.getStore();
        //        BlockCluster b = hFE.getRoot().getBC(1, 57);
//        System.out.println("b= "+b);
//        System.out.println("bFE=\n"+b.getRM().getFullMatrix());
        hFE.printOnJFrame("FE");
        //
        // Complex assembly                
//        StorageHmatrixComplex hc = new StorageHmatrixComplex(1.0, hFE, 50.0, hIV, null);
        StorageHmatrixComplex hc = new StorageHmatrixComplex(a, b, hFE, c, d, hIV, null);
        //        System.out.println(hc.getRoot().getSon(1).getRM());
        System.out.println("dof2idx = " + Arrays.toString(hc.getDof2Idx()));
        hc.printOnJFrame("RL Complex");
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        double vr[] = new double[2 * C.getActiveDofCount()];
        Arrays.fill(vr, 1.0);
        ColumnVector v = new ColumnVector(vr);
        ColumnVector rb = new ColumnVector(M.product(vr, new double[M.getRows() * 2]));
        ColumnVector rh = hc.prod(v);
        RkMatrixComplex rk = hc.getRoot().getSon(1).getRM().copy();
        rh.sub(rb);
        System.out.println("Erreur Assemblage Complexe = " + ComplexOperation.norm2(rh) / ComplexOperation.norm2(rb));
        hc.Agglomerate(new TruncationControlComplex("rel", 1e-4));
        RkMatrixComplex rtk = hc.getRoot().getSon(1).getRM().copy();
        hc.printOnJFrame("RL agglo'd");
        rh = hc.prod(v);
        rh.sub(rb);
        System.out.println("Erreur after agglo(1e-1) complexe = " + ComplexOperation.norm2(rh) / ComplexOperation.norm2(rb));
        //
        //
        //
        //
        //
        //
//        System.out.println("r.S = " + rk.getS());
        System.out.println("r2.S = " + rtk.getS());

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
