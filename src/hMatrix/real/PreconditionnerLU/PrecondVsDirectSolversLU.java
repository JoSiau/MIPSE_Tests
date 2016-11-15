/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.PreconditionnerLU;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.Preconditionners.HmatrixLUDecomposition;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import got.matrix.ColumnVector;
import got.matrix.LuDecomposition;
import got.matrix.Matrix;

import java.io.IOException;

/**
 * @author jsiau
 */
public class PrecondVsDirectSolversLU {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        // Load the mesh
        String meshDir = new java.io.File(".").getCanonicalPath();
        //*
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
        ImportFlux mesh = new ImportFlux(meshDir + "/sphere/SPHERE_" + "3538.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();
        /*/
         meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere_gmsh";
         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/SPHERE_4881.msh");                    
         ElementSetHomogene ES = mesh1.createHomogeneESet();
         //*/

        int d = 3;
        int n = ES.getNbElement();
        int nbGauss = 7;
        /*
         SECOND MEMBRE
         */
        ColumnVector secondMembre = new ColumnVector(n);
        for (int i = 0; i < n; i++) {
            secondMembre.setElement(i, Math.random() * 100);
        }
        /*
         HMATRIX
         */
        double eps = 1e-5;
        TruncationControl tol = new TruncationControl("rel", eps);
        int kmax = 50, nmin = 30;
        Cell C = new Cell(ES);

        long deb = System.currentTimeMillis();
        GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss,
                eps, kmax, nmin);
        f.assembly();
        StorageHmatrix H0 = (StorageHmatrix) f.getStore();
        H0.Agglomerate(tol);

        ColumnVector resD = new ColumnVector(n);
        HmatrixLUDecomposition hLu = new HmatrixLUDecomposition(H0, tol);
        hLu.solve(resD, secondMembre);
        long fin = System.currentTimeMillis();
        double tpsD = (fin - deb) * 1e-3;
        System.err.println("Total time to solve direct = " + tpsD);


        deb = System.currentTimeMillis();
        double epsPrecond = 1e-1;
        TruncationControl tolPrecond = new TruncationControl("rel", epsPrecond);
        f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss,
                eps, kmax, nmin);
        f.assembly();
        H0 = (StorageHmatrix) f.getStore();
        // Preconditionneur
        StorageHmatrix H1 = H0.AgglomerateACopy(epsPrecond);
        hLu = new HmatrixLUDecomposition(H1, tolPrecond);
        // Solveur
        FGMResReal t1 = new FGMResReal(H0, hLu);
        double configRes[] = {100, eps, 1, -50};
        t1.setInfoResolution(configRes);
        double[] resPv = new double[n];
        double[] secMem = new double[n];
        secondMembre.get(secMem);
        t1.solve(resPv, secMem);
        ColumnVector resP = new ColumnVector(resPv);
        fin = System.currentTimeMillis();
        double tpsP = (fin - deb) * 1e-3;
        System.err.println("Total time to solve with preconditionner = " + tpsP);
        /*
         FULL STORAGE
         */
        System.out.println("Calcul du stockage plein");
        deb = System.currentTimeMillis();
        GalerkinIntegralFormulation IFCell = new GalerkinIntegralFormulationFull(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss);
        IFCell.assembly();
        Matrix Ai = ((StorageFull) IFCell.getStore()).getMatrix();
        LuDecomposition MD = new LuDecomposition(Ai);
        ColumnVector resEx = new ColumnVector(n);
        MD.solve(secondMembre, resEx);
        fin = System.currentTimeMillis();
        System.out.println("Total time Full Storage = " + (fin - deb) * 1e-3);

        /*
         CHECK ERROR
         */
        resD.sub(resEx);
        System.err.println("Relative error with direct solver = " + resD.norm() / resEx.norm());
        resP.sub(resEx);
        System.err.println("Relative error with preconditionned solver = " + resP.norm() / resEx.norm());
        System.err.println("Total time to solve direct = " + tpsD);
    }
}
