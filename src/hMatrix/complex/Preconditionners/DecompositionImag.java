/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.complex.Preconditionners;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.elements.Hexaedre;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.KernelInterface;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.ComplexOperation;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.Preconditionners.HmatrixLUDecompositionComplex;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.TruncationControlComplex;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.decompositions.matlab.DecompositionMatLab;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.mipseCore.storage.StorageHmatrixComplex;
import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D;
import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.SolverLUBasic2D;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.Arrays;

/**
 * Test la décomposition H-LU d'une matrice purement imaginaire !
 *
 * @author jsiau
 */
public class DecompositionImag {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(1);
        long t;
        System.out.println("DecompositionImag.java");
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/";
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
        GalerkinIntegralFormulationFull full = new GalerkinIntegralFormulationFull(C, C, k, new SelfElementFixedGauss(nbPGSrc, new InCreasedPGSourceNumber(nbPGTrg)), nbPGSrc);
        full.assembly();
        Matrix m = ((StorageFull) full.getStore()).getMatrix();

        double mat[][] = new double[C.getActiveDofCount()][C.getActiveDofCount() * 2];
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < C.getActiveDofCount(); j++) {
                mat[i][2 * j + 1] = m.getElement(i, j);
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
        // Complex assembly                
        StorageHmatrixComplex hc = new StorageHmatrixComplex(null, hIV, null);
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
        System.out.println("Erreur= " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

        TruncationControlComplex tolDec = new TruncationControlComplex("rel", 1e-4);
        hc.Agglomerate(tolDec);
        hc.printOnJFrame();

        r_hm = hc.prod(rhsv);
        r_hm.sub(r_exact);
        System.out.println("Erreur after agglo= " + ComplexOperation.norm2(r_hm) / ComplexOperation.norm2(r_exact));

//        System.out.println("hm =\n" + hc.getRoot().getFullMatrix());
//        System.out.println("mat =\n" + Basic2D.toString(mat));
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        HmatrixLUDecompositionComplex hlu = new HmatrixLUDecompositionComplex(hc, tolDec);
        x_hm = new ColumnVector(2 * n);
        t = System.currentTimeMillis();
        hlu.solve(x_hm, rhsv);
        System.out.println("Time to solve H-LU = " + (System.currentTimeMillis() - t) + " msec.");
        System.out.println("x_hm = " + x_hm.transpose());

        t = System.currentTimeMillis();
        SolverLUBasic2D lu = new SolverLUBasic2D(mat, 1e-20, false, false, false);
        System.out.println("Time to decompose LU = " + (System.currentTimeMillis() - t) + " msec.");
        t = System.currentTimeMillis();
        x_exact = new ColumnVector(lu.solveLU(vr));
        System.out.println("Time to solve LU = " + (System.currentTimeMillis() - t) + " msec.");
        System.out.println("x_exact = " + x_exact.transpose());

        x_hm.sub(x_exact);
        System.out.println("Erreur relative = " + ComplexOperation.norm2(x_hm) / ComplexOperation.norm2(x_exact));

//        StorageHmatrixComplex l = hlu.getL();
//        StorageHmatrixComplex u = hlu.getU();
        DecompositionMatLab.exit();
//        printHmatrix.closeAll();
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
