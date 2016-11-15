/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.complex.imp;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.elements.Hexaedre;
import g2elab.mipse.meshCore.elements.QuadrangleDroit;
import g2elab.mipse.meshCore.elements.TriangleDroit;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.KernelInterface;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.Preconditionners.HmatrixLUDecompositionComplex;
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
import g2elab.mipse.numericalTools.vector.full.VectorFull;
import g2elab.mipse.numericalTools.vector.full.VectorFullComplex;
import got.matrix.Matrix;
import java.io.IOException;

/**
 *
 * @author jsiau
 */
public class NewMain1 {

    /**
     * TRY TO UNDERSTAND SOMETHING !
     *
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/hMatrix/mesh/";
        ImportFlux mesh = new ImportFlux(meshDir + "sphere/SPHERE_572.DEC");
//        ImportFlux mesh = new ImportFlux(meshDir + "4CUBES_4x4.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();
//        ES = new ElementSetHomogene(ES.getBorder(0).getElements());
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
        } else if (ES.getElements(0) instanceof QuadrangleDroit) {
            nbPGSrc = 4;
            nbPGTrg = 9;
        } else if (ES.getElements(0) instanceof TriangleDroit) {
            nbPGSrc = 4;
            nbPGTrg = 9;
        } else {
            throw new InternalError("Fait les PG des " + ES.getElements(0));
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
                mat[i][2 * j] = S.getQuick(i, j);
                mat[i][2 * j + 1] = m.getElement(i, j);
            }
        }
        Basic2D M = new Basic2D(mat);
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        int kmax = 50, nmin = 16;
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
        hc.printOnJFrame("ass");
        ////////////////////////
        BlockClusterComplex bc = hc.getRoot().getSon(new int[]{0, 2, 0, 0});
        System.out.println("bc = " + bc.toString());
        bc.printOnJFrame("lll",hc.getRows());
        //
        //
//        for (double er : new double[]{1e-1 / 3, 1e-1 / 2.5, 1e-1 / 2, 1e-1 / 1.5}) {
//            StorageHmatrixComplex ha = NewMain1.agglo(hc.copy(true), er);
//            ha.printOnJFrame("agglo " + er);
//            ////////////////////////
//            bc = hc.getRoot().getSon(new int[]{0, 2});
//            System.out.println("bc = " + bc.toString());
//        }
        ////////////////////////////////////////////////////////////////////////
        Basic2D mh = hc.reOrder_Full2HM(M);
        //
        Basic2D minih = mh.submatrix(bc.getRowStart(), bc.getColStart(), bc.getRowSize(), bc.getColSize());
        
        for (double er : new double[]{1e-1 / 3, 1e-1 / 2.5, 1e-1 / 2, 1e-1 / 1.5}) {
            System.out.println("er = "+er);
            BlockClusterComplex bb = bc.copy(true);
            bb.agglomerate(new TruncationControlComplex("rel", er));
            bb.printOnJFrame("agglo " + er, hc.getRows());
            Basic2D tmp = bb.getM();
            if(bb.isAdm()){
                System.out.println("Rank = "+bb.getRM().getRank());
            }
            tmp.sub(minih);
            System.out.println("Erreur agglo = "+tmp.norm2()/minih.norm2());
        }
        //
        VectorFullComplex b = VectorFullComplex.generateRandom(C.getActiveDofCount());
        VectorFullComplex xh = new VectorFullComplex(hc.product(b.getValues(false), new double[2 * C.getActiveDofCount()]));
        VectorFullComplex xm = new VectorFullComplex(M.product(b.getValues(false), new double[2 * C.getActiveDofCount()]));
        VectorFull diff = xh.sub(xm, null);
        System.out.println("Error assembly = " + diff.norm2() / xm.norm2());
        
        bc = hc.getRoot().getSon(new int[]{0, 2, 0});//, 0, 3});
        bc.printOnJFrame("tt", hc.getRows());

    }

    public static StorageHmatrixComplex agglo(StorageHmatrixComplex h, double eps) {
//        StorageHmatrixComplex hc = h.copy(true);
        TruncationControlComplex tol = new TruncationControlComplex("rel", eps);
        h.agglomerate(tol);
        return h;
    }

}
