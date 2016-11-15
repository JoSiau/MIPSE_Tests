/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Arithmetics;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.KernelInterface;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.Preconditionners.HmatrixLUDecomposition;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.StorageBlock.Admissible.RkMatrix;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import g2elab.mipse.numericalTools.matrix.exportMatlab.MatlabExporter;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.LuDecomposition;
import got.matrix.Matrix;
import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author jsiau
 */
public class tmp {

    private static int nbPGSrc;
    private static int nbPGTrg;

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(2);
        System.out.println("TestAgglo.java");
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/";

        deb = System.nanoTime();
        ImportFlux mesh1 = new ImportFlux(meshDir + "HCAFARMESH_2.DEC");
        ElementSetHomogene ES = new ElementSetHomogene(mesh1.getRegion(0).getElementSet().computeBorder(0).getElements());
        fin = System.nanoTime();
//        ES.plotElementSet("d:/tmp/tmp.vtk");
//        Exec.openFile("d:/tmp/tmp.vtk");
        System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);

        int d = 3;
        int n = ES.getNbElement();
        System.out.println(" n = " + n + "\t d = " + d);
        Cell C = new Cell(ES);
        KernelInterface k = new MultG();
        nbPGSrc = 4;
        nbPGTrg = 9;

        ////////////////////////////////////////////////////////////////////////
        int kmax = 1, nmin = 24;
        double eps = 1e-4, eta = 0.0;
        StorageHmatrix hIV = comp(C, kmax, nmin, eps, eta, k);
        hIV.printOnJFrame("1");
        System.out.println("dof2idx = " + Arrays.toString(hIV.getDof2Idx()));
        Matrix m = hIV.getRoot().getSon(1).getFullMatrix();

        ////////////////////////////////////////////////////////////////////////
        StorageHmatrix hIV2 = comp(C, 1, nmin, eps, 2.0, k);
        System.out.println("dof2idx'' = " + Arrays.toString(hIV2.getDof2Idx()));
        RkMatrix r = hIV2.getRoot().getSon(1).getRM();
        hIV2.printOnJFrame("2");
        Matrix m2 = r.getFullMatrix();
        m2.sub(m);
        System.out.println("error2 = " + m2.norm() / m.norm());

        ////////////////////////////////////////////////////////////////////////
        StorageHmatrix hIV3 = comp(C, 2, nmin, eps, 2.0, k);
        System.out.println("dof2idx'' = " + Arrays.toString(hIV3.getDof2Idx()));
        RkMatrix r3 = hIV3.getRoot().getSon(1).getRM();
        hIV3.printOnJFrame("3");
        Matrix m3 = r3.getFullMatrix();
        m3.sub(m);
        System.out.println("error3 = " + m3.norm() / m.norm());

        ////////////////////////////////////////////////////////////////////////
        StorageHmatrix hIV4 = comp(C, 3, nmin, eps, 2.0, k);
        System.out.println("dof2idx'' = " + Arrays.toString(hIV4.getDof2Idx()));
        RkMatrix r4 = hIV4.getRoot().getSon(1).getRM();
        hIV4.printOnJFrame("3");
        Matrix m4 = r4.getFullMatrix();
        m4.sub(m);
        System.out.println("error4 = " + m4.norm() / m.norm());
        ////////////////////////////////////////////////////////////////////////
        StorageHmatrix hIV5 = comp(C, 10, nmin, 1e-6, 2.0, k);
        System.out.println("dof2idx'' = " + Arrays.toString(hIV5.getDof2Idx()));
        RkMatrix r5 = hIV5.getRoot().getSon(1).getRM();
        hIV5.printOnJFrame("4");
        Matrix m5 = r5.getFullMatrix();
        m5.sub(m);
        System.out.println("error5 = " + m5.norm() / m.norm());
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        System.out.println("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
        System.out.println("$$$$$$$$$$$$$$$$$$  >LU< $$$$$$$$$$$$$$$$$$$$$$$$$");
        double epsLU = 3e-1;
        StorageHmatrix matLU = hIV5.copy(true);
        matLU.Agglomerate(new TruncationControl("rel", epsLU));
        matLU.printOnJFrame("matLU");
        HmatrixLUDecomposition hlu = new HmatrixLUDecomposition(matLU, new TruncationControl("rel", epsLU));
        // RHS
        double bv[] = new double[n];
        for (int i = 0; i < n; i++) {
            bv[i] = Math.random() * 100000;
        }
        ColumnVector b = new ColumnVector(n, bv);
        // Solution
        ColumnVector x = new ColumnVector(n);
        hlu.solve(x, b);
        //
        // Check error !
        ColumnVector b2 = hIV5.prod(x);
        b2.sub(b);
        System.err.println("Error lu Ax-b = " + b2.norm() / b.norm());
        //
        // FULL
        GalerkinIntegralFormulationFull f = new GalerkinIntegralFormulationFull(C, C, k, new SelfElementFixedGauss(nbPGSrc, new InCreasedPGSourceNumber(nbPGTrg)), nbPGSrc);
        f.assembly();
        Matrix F = ((StorageFull) f.getStore()).getMatrix();
        Matrix Fh = hIV5.reOrder_Full2HM(F);
        Matrix Hm = hIV5.getApprox();
        Hm.sub(Fh);
        System.err.println("Error Ass. = " + Hm.norm() / Fh.norm());
        // LU Full
        LuDecomposition lu = new LuDecomposition(F);
        ColumnVector x_ref = new ColumnVector(n);
        lu.solve(b, x_ref);
        x.sub(x_ref);
        System.err.println("Error lu = " + x.norm() / x_ref.norm());

        System.out.println("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
        System.out.println("$$$$$$$$$$$$$$$$$$  >PREC< $$$$$$$$$$$$$$$$$$$$$$$$$");
        FGMResReal solver = new FGMResReal(hIV5, null);
        solver.setInfoResolution(new double[]{1, -1e-8, 1, -1});
        double[] xv1 = solver.solve(new double[n], bv);
        //
        solver.setProduitPrecond(hIV5, hlu);
        double[] xv2 = solver.solve(new double[n], bv);

        ColumnVector x1 = new ColumnVector(xv1);
        ColumnVector x2 = new ColumnVector(xv2);
        x2.sub(x1);
        System.out.println("Error between solvers = " + x2.norm() / x1.norm());

        x1.sub(x_ref);
        System.out.println("Error it./directe = " + x1.norm() / x_ref.norm());
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    private static StorageHmatrix comp(Cell C, int kmax, int nmin, double eps, double eta, KernelInterface k) {
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, C, k, new SelfElementFixedGauss(nbPGSrc, new InCreasedPGSourceNumber(nbPGTrg)), nbPGSrc, nbPGSrc,
                eps, kmax, nmin, 5, eta);
        f.assembly();
        return (StorageHmatrix) f.getStore();
    }

}
