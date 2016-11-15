/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Arithmetics;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMesh;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.IOException;

/**
 * Test l'agglomeration.
 *
 * @author jsiau
 */
public class TestAgglomeration {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        String meshDir = new java.io.File(".").getCanonicalPath();
//        /*
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
        ImportFlux mesh = new ImportFlux(meshDir + "/sphere/SPHERE_3538.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();
         /*/
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
        ImportGmshMesh mesh1 = new ImportGmshMesh(meshDir + "/plaque/plaque_4548.msh");
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        //*/

        int d = 3;
        int n = ES.getNbElement();
        System.out.println(" n = " + n + "\t d = " + d);

        Cell C = new Cell(ES);
        int kmax = 50, nmin = 30;
        int order = 4;
        double eps = Math.pow(10, -order - 1);
        /*
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, C, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3, 3,
                eps, kmax, nmin, order);
        /*/
        GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3,
                eps, kmax, nmin,2.0,false);
        //*/
        f.assembly();
        StorageHmatrix H = (StorageHmatrix) f.getStore();
        
        /*
        ERREUR DE DEPART
        */
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
        H.CheckError();
        System.out.println("Erreur relative de depart= " + vH.norm() / ex.norm());

        H.printOnJFrame();
        H.Agglomerate(new TruncationControl("rel", 1e-5));
        vH = H.prod(v);
        vH.sub(ex);
        H.CheckError();
        System.out.println("Erreur relative 1e-5= " + vH.norm() / ex.norm());
        H.printOnJFrame("1e-5");
        
        
        H.Agglomerate(new TruncationControl("rel", 1e-3));
        vH = H.prod(v);
        vH.sub(ex);
        H.CheckError();
        System.out.println("Erreur relative 1e-3= " + vH.norm() / ex.norm());
        H.printOnJFrame("1e-3");
        
        
        H.Agglomerate(new TruncationControl("rel", 1e-2));
        vH = H.prod(v);
        vH.sub(ex);
        H.CheckError();
        System.out.println("Erreur relative 1e-2= " + vH.norm() / ex.norm());
        H.printOnJFrame("1e-2");

    }

}
