/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package hMatrix.real.Arithmetics;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
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
 * Test the MulAdd() of Hmatrices
 *
 * @author jsiau
 */
public class MullAddHmat {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        String meshDir = new java.io.File(".").getCanonicalPath();
        //*
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
        ImportFlux mesh = new ImportFlux(meshDir + "/sphere/SPHERE_572.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();
        /*/
         meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/";
         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/plaque/plaque_576.msh");                    
         ElementSetHomogene ES = mesh1.createHomogeneESet();
         //*/
        int d = 3;
        int n = ES.getNbElement();
        int nbGauss = 7;
        
        /*
        HMATRIX
        */
        double eps = 1e-5;
        TruncationControl tol = new TruncationControl("rel", eps);
        int kmax = 50, nmin = 30;
        Cell C = new Cell(ES);
        GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss,
                eps, kmax, nmin, 2.0);
        f.assembly();
        StorageHmatrix H0 = (StorageHmatrix) f.getStore();
//        H0.Agglomerate(tol);


        long deb = System.currentTimeMillis();
        StorageHmatrix Hp = H0.copy(false);
        Hp.mulAdd(H0, H0, tol);
        long fin = System.currentTimeMillis();
        System.out.println("\nTime to compute the mulAdd= " + (fin - deb) * 1e-3);

        if (Hp.getRoot().isThereAnyNaN()) {
            System.out.println("WHAAAAAT !");
        }

        GalerkinIntegralFormulationFull full = new GalerkinIntegralFormulationFull(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss);
        full.assembly();
        Matrix Mp = ((StorageFull) full.getStore()).getMatrix();
        Matrix M = new Matrix(n, n);
        deb = System.currentTimeMillis();
        M.mulAdd(1, Mp, Mp);
        fin = System.currentTimeMillis();
        System.out.println("\nTime to compute the mulAdd FULL= " + (fin - deb) * 1e-3);

        ColumnVector v = new ColumnVector(n);
        for (int i = 0; i < n; i++) {
            v.setElement(i, Math.random() * 100);
        }

        ColumnVector vp = new ColumnVector(n);
        vp.mul(M, v);

        ColumnVector vh = Hp.prod(v);

        vh.sub(vp);
        System.out.println("\nErreur relative par PMV= " + vh.norm() / vp.norm());
        System.out.println("");
                
        
        
        /*
         CHECK THE ERROR BLOCK-WISE
         */
        Matrix M2 = new Matrix(n, n);
        int ind[] = H0.getDof2Idx();
        for (int i = 0; i < M.getRowCount(); i++) {
            for (int j = 0; j < M.getColumnCount(); j++) {
                M2.setElement(i, j, M.getElement(ind[i], ind[j]));
            }
        }
        Hp.CheckError(M2, eps);

//        Hp.printOnJFrame();
        H0.printOnJFrame();

    }

}
