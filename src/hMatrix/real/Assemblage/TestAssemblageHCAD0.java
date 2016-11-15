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
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.IOException;

/**
 * @author jsiau
 */
public class TestAssemblageHCAD0 {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        System.out.println("TestAssemblageHCAD0.java");
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
        int kmax = 50, nmin = 30;
        int order = 8;
        double eps = Math.pow(10, -order - 1);
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, C, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3, 3,
                eps, kmax, nmin, order, 2.0, true);
        f.assembly();
        StorageHmatrix H = (StorageHmatrix) f.getStore();

//        H.Coarsen(new TruncationControl("rel",eps));
        GalerkinIntegralFormulationFull IF = new GalerkinIntegralFormulationFull(C, C, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3);
        IF.assembly();
        Matrix M = ((StorageFull) IF.getStore()).getMatrix();
        ColumnVector v = new ColumnVector(n);
        for (int i = 0; i < n; i++) {
            v.setElement(i, Math.random() * 1000);
        }
        v.setAllElements(1);
        ColumnVector ex = new ColumnVector(n);
        ex.mul(M, v);

        ColumnVector vH = new ColumnVector(n);
        vH = H.prod(v);

        vH.sub(ex);
        System.out.println("Erreur relative= " + vH.norm() / ex.norm());

        //Affiche la Hmatrix
        H.printOnJFrame();
    }
}
