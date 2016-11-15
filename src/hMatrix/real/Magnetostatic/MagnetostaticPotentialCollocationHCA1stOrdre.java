/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Magnetostatic;

import g2elab.mipse.analytical.singularity.SimpleTetraedreGradientPotential;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.GradNodalDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.NegDotDG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.HCA.HmatrixHCAMagnetoStatPotentialCollocDeg1;
import g2elab.mipse.mipseCore.numericalMethods.CollocationIntegralFormulation;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;

/**
 * A FILE TO TEST THE ASSEMBLAGE OF AN H-MATRIX HCA-P1 (BY COLLOCATION) FOR THE FORMULATION OF ANTO
 * <p/>
 * THIS FILE MAY BE OBSOLETE...
 *
 * @author siau
 */
public class MagnetostaticPotentialCollocationHCA1stOrdre {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        File fi = new File("");
        String path = fi.getAbsolutePath();
        String file = path + "./../formulations/src/g2elab/mipse/formulationInProgress/magnetostatic/potential/prob13_mesh_5718.DEC";
        ImportFlux mesh = new ImportFlux(file);
        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();

        int d = 3;
        int n = ES.getNbNoeud();

        System.out.println("nbNoeuds= " + ES.getNbNoeud());
        System.out.println("nbElmts= " + ES.getNbElement());
                
        
        /*
         * H-MATRIX
         */
        GradNodalDeg1 gradAlpha = new GradNodalDeg1(ES);
        int epsPow = 6;
        HmatrixHCAMagnetoStatPotentialCollocDeg1 f = new HmatrixHCAMagnetoStatPotentialCollocDeg1(gradAlpha, new NegDotDG(), new SimpleTetraedreGradientPotential(), 15,
                Math.pow(10, -epsPow), 50, 30, epsPow - 1, 2.0, true);
        double deb = System.nanoTime();
        StorageHmatrix H = new StorageHmatrix(f);
        double fin = System.nanoTime();
        System.out.println("Time to assembly the Hmatrix: " + (fin - deb) / 1e9);

//        Hmatrix HP = H.copy(true);        
//        deb = System.nanoTime();
//        TruncationControl tol =  new TruncationControl("rel",1e-1);
//        HP.Coarsen(tol);
//        HmatrixLUDecomposition Hlu = new HmatrixLUDecomposition(HP, tol);
//        fin = System.nanoTime();
//        System.out.println("Time to construct the Preconditionner: "+(fin-deb)/1e9);
        
        /*
         * STOCKAGE PLEIN
         */
        CollocationIntegralFormulation IntegralTerm = new CollocationIntegralFormulation(new GradNodalDeg1(ES), new NegDotDG(), new SelfElementFixedGauss(15, new AnalyticalCorrection()));
        IntegralTerm.matAssembly(new GradNodalDeg1(ES));
        Matrix M = ((StorageFull) IntegralTerm.getStore()).getMatrix();
        
        /*
         * Verification par PMV
         */
        ColumnVector v = new ColumnVector(n);
        for (int i = 0; i < n; i++)
            v.setElement(i, Math.random() * 100);

        ColumnVector exV = new ColumnVector(n);
        exV.mul(M, v);

        ColumnVector hV = H.prod(v);

        hV.sub(exV);
        System.out.println("eps= " + Math.pow(10, -epsPow) + "\t Erreur rel par PMV= " + hV.norm() / exV.norm());
    }


}
