/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Assemblage;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMesh;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.IOException;

/**
 * TEST TO THE ASSEMBLY OF AN HMATRIX WITH ACA's COMPRESSION
 *
 * @author jsiau
 */
public class TestAssemblageACA {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        System.out.println("TestAssemblageACA.java");
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/plaque";

        deb = System.nanoTime();
        ImportGmshMesh mesh1 = new ImportGmshMesh(meshDir + "/plaque_4548.msh");
        mesh1.meshSummary();
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        fin = System.nanoTime();
        System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);

        int d = 3;
        int n = ES.getNbElement();
        System.out.println(" n = " + n + "\t d = " + d);

        Cell C = new Cell(ES);

        if (C.equals(C))
            System.out.println("1");

        if (C == C)
            System.out.println("1b");

        GalerkinIntegralFormulationFull IF = new GalerkinIntegralFormulationFull(C, C, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3);
        IF.assembly();
        Matrix M = ((StorageFull) IF.getStore()).getMatrix();

        ColumnVector v = new ColumnVector(n);
        for (int i = 0; i < n; i++) {
            v.setElement(i, Math.random() * 1000);
        }
//        v.setAllElements(1);
        ColumnVector ex = new ColumnVector(n);
        ex.mul(M, v);

        int kmax = 50, nmin = 30;
        double eps = 1e-8;
        GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3,
                eps, kmax, nmin, 2.0, false);
        //*
        f.assembly();
        StorageHmatrix H = (StorageHmatrix) f.getStore();
        /*/
        Hmatrix H = new Hmatrix(f);//, ProgressBarHmat.getAnimateur());
        H.Coarsen(new TruncationControl("rel",1e-5));
        //*/
        if (H.getRoot().isThereAnyNaN())
            System.err.println("WTF NAN");

        ColumnVector vH = H.prod(v);

        vH.sub(ex);
        System.out.println("Erreur relative= " + vH.norm() / ex.norm());
        
        
        /*
         CHECK THE ERROR BLOCK-WISE
         */
        Matrix M2 = new Matrix(n, n);
        int ind[] = H.getDof2Idx();
        for (int i = 0; i < M.getRowCount(); i++) {
            for (int j = 0; j < M.getColumnCount(); j++) {
                M2.setElement(i, j, M.getElement(ind[i], ind[j]));
            }
        }
        H.CheckError(M2, eps);

        H.printOnJFrame();

        H.Agglomerate(new TruncationControl("rel", 1e-5));
        //Affiche la Hmatrix
        H.printOnJFrame();
    }

}
