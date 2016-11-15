/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.complex.imp;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.elements.Hexaedre;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.KernelInterface;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.Graphics.printHmatrix;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.Preconditionners.HmatrixLUDecompositionComplex;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.TruncationControlComplex;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.complex.decompositions.matlab.DecompositionMatLab;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulation;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulationHmatrix;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.mipseCore.storage.StorageHmatrixComplex;
import g2elab.mipse.mipseCore.storage.StorageSparse;
import g2elab.mipse.numericalTools.iterativeSolver.complex.FGMResComplex;
import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
import g2elab.mipse.numericalTools.vector.full.VectorFullComplex;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author jsiau
 */
public class SolutionPrecond {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(1);
        long t;
        System.out.println("SolutionPrecond.java");
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
//        ImportFlux mesh1 = new ImportFlux(meshDir + "HCAFARMESH.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        //*/
        fin = System.nanoTime();
        System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);

        int d = 3;
        //*
        FaceDeg1 C = new FaceDeg1(ES);
        KernelInterface k = new MultGvect();
        /*/
         Cell C = new Cell(ES);
         KernelInterface k = new MultG();
         //*/
        int n = C.getActiveDofCount();
        int nbPGSrc, nbPGTrg;
        if (ES.getElements(0) instanceof Hexaedre) {
            nbPGSrc = 8;
            nbPGTrg = 27;
        } else {
            nbPGSrc = 4;
            nbPGTrg = 7;
        }
        TruncationControlComplex tolAss = new TruncationControlComplex("rel", 1e-4);
        TruncationControlComplex tolPrecond = new TruncationControlComplex("rel", 1e-1);

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        FiniteElementFormulation IF = new FiniteElementFormulation(C, C);
        IF.assembly(nbPGSrc);
        SparseMatrixRowReal S = ((StorageSparse) IF.getStore()).getSparseMat();

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
        /*
         BlockClusterTree bct = new BlockClusterTree(C, S, 3, nmin, eta);
         StorageHmatrix hFE = new StorageHmatrix(bct);
         /*/
        FiniteElementFormulationHmatrix Ifh = new FiniteElementFormulationHmatrix(C, nmin, eta);
        Ifh.assembly(nbPGSrc);
        StorageHmatrix hFE = (StorageHmatrix) Ifh.getStore();
        //*/
        //
        // Complex assembly                
        StorageHmatrixComplex hc = new StorageHmatrixComplex(hFE, hIV, null);
        hc.printOnJFrame();

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        hc.getRoot().checkMatricesSize();
        hc.Agglomerate(tolAss);
        hc.printOnJFrame("Hmatrix solution");
        hc.getRoot().checkMatricesSize();

        StorageHmatrixComplex hcP = hc.copy(true);
        hcP.Agglomerate(tolPrecond);
        hcP.getRoot().checkMatricesSize();
        hcP.printOnJFrame("Hmatrix precond.");
        HmatrixLUDecompositionComplex hlu = new HmatrixLUDecompositionComplex(hcP, tolPrecond);

        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        System.out.println("\nSolution with precond. !\n");
        FGMResComplex gmres = new FGMResComplex(hc, hlu);

        double x[] = new double[2 * n];
        double b[] = new double[2 * n];
        Arrays.fill(b, 1.1);
        gmres.setInfoResolution(new double[]{10, 1e-8, 1.0, -1});
        t = System.currentTimeMillis();
        x = gmres.solve(x, b);
        t = System.currentTimeMillis() - t;
        System.out.println("Time to solve = " + t + " mSec.");
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////
        System.out.println("\nSolution without precond. !\n");
        gmres = new FGMResComplex(hc, null);

        double xe[] = new double[2 * n];
        gmres.setInfoResolution(new double[]{10, 1e-8, 1.0, -1});
        t = System.currentTimeMillis();
        xe = gmres.solve(xe, b);
        t = System.currentTimeMillis() - t;
        System.out.println("Time to solve = " + t + " mSec.");

        VectorFullComplex v = new VectorFullComplex(x);
        VectorFullComplex ve = new VectorFullComplex(xe);

        v.sub(ve, v);
        System.out.println("Erreur solution = " + v.norm2() / ve.norm2());

        System.out.println("Reorganisation time = " + DecompositionMatLab.getReOrganisationTime() + " msec");

        DecompositionMatLab.exit();
        printHmatrix.closeAll();
        GestionnaireTaches.getGestionnaireTaches().stop();

    }

}
