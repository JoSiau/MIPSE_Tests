/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.KernelInterface;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.Preconditionners.HmatrixLUDecomposition;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulation;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulationHmatrix;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.mipseCore.storage.StorageSparse;
import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author jsiau
 */
public class TestFiniteElementAssembly {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(1);
        System.out.println("TestFiniteElementAssembly.java");
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/";

        deb = System.nanoTime();
        /*
         ImportGmshMesh mesh1 = new ImportGmshMesh(meshDir + "plaque/plaque_4548.msh");
         mesh1.meshSummary();
         ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
         /*/
        ImportFlux mesh1 = new ImportFlux(meshDir + "sphere/SPHERE_2104.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        //*/
        fin = System.nanoTime();
        System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);

        int d = 3;
        int n = ES.getNbElement();
        System.out.println(" n = " + n + "\t d = " + d);
        //*
        FaceDeg1 C = new FaceDeg1(ES);
        KernelInterface k = new MultGvect();
        /*/
         Cell C = new Cell(ES);
         KernelInterface k = new MultG();
         //*/

        FiniteElementFormulation IF = new FiniteElementFormulation(C, C);
        long t = System.currentTimeMillis();
        IF.assembly(1);
        long t2 = System.currentTimeMillis();
        System.out.println("Time to assemble FE maatrix = " + (t2 - t) + " msec.");
        int kmax = 50, nmin = 30;
        double eps = 1e-4;
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, C, k, new SelfElementFixedGauss(4, new InCreasedPGSourceNumber(7)), 4, 4,
                eps, kmax, nmin, 4, 2.0);
        //*
        f.assembly();
        StorageHmatrix H = (StorageHmatrix) f.getStore();
//        System.out.println("f.IV= " + H.getRoot().getFullMatrix());
        H.printOnJFrame("H");

        SparseMatrixRowReal S = ((StorageSparse) IF.getStore()).getSparseMat();
        System.out.println("S.Mem = " + S.getMemoryUsed());
//        System.out.println("S:\n" + S.getFullMatrix());

        FiniteElementFormulation Ifh = new FiniteElementFormulationHmatrix(C, nmin, 2.0);
        Ifh.assembly(1);
        StorageHmatrix h = (StorageHmatrix) Ifh.getStore();

        System.out.println(Arrays.toString(h.getDof2Idx()));

//        System.out.println("f.FE= " + h.getRoot().getFullMatrix());
        h.printOnJFrame("FE");
        h.addH(H.copy(true), new TruncationControl("rel", eps));
//        System.out.println("f.FE= " + h.getRoot().getFullMatrix());
//        h.printOnJFrame("FE + H");

//        h.Agglomerate(new TruncationControl("rel", 1e-4));
//        h.printOnJFrame();
//        System.out.println("h: " + h.getApprox());
        double x[] = new double[C.getActiveDofCount()];
        Arrays.fill(x, 1.0);

        double r[] = S.product(x, new double[C.getActiveDofCount()]);
        r = H.product(x, r);
//        System.out.println("r=" + Arrays.toString(r));

        double rf[] = h.product(x, new double[C.getActiveDofCount()]);
//        r = H.product(x, r);
//        System.out.println("rf=" + Arrays.toString(rf));

        ColumnVector rv = new ColumnVector(r);

        ColumnVector rvf = new ColumnVector(rf);

        rvf.sub(rv);
        System.out.println("Erreur= " + rvf.norm() / rv.norm());

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        TruncationControl tolPrecond = new TruncationControl("rel", 3e-1);
        h.Agglomerate(tolPrecond);
        h.printOnJFrame("h agglo");
        HmatrixLUDecomposition lu = new HmatrixLUDecomposition(h, tolPrecond);

        GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
