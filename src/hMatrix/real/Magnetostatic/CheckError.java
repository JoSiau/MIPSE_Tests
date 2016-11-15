/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
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

/**
 * @author jsiau
 */
public class CheckError {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // Sphere Creuse
        System.out.println("hhhhh");
        String path = "D:";
        // SC_3846 SC_6955 SC_12938 SC_26420 SC_48154 SC_79380
        String file = "D:/Meshs/SphereCreuse/SC_";
        //*/

//        int nFile[] = new int[6];
//        nFile[0] = 3846;
//        nFile[1] = 6955;
//        nFile[2] = 12938;
//        nFile[3] = 26420;
//        nFile[4] = 48154;
//        nFile[5] = 79380;
//        boolean bool = true;
//
//        double tHmat[] = new double[6];
//        double tPlein[] = new double[6];
//        double nNodes[] = new double[6];
//        double nElmts[] = new double[6];
//        
//        int nF = nFile.length-3;
//        for (int k = 0; k < 1; k++) {
        int nFile[] = new int[3];
        nFile[0] = 3846;
        nFile[1] = 12938;
        nFile[2] = 48154;
        boolean bool = true;

        double tHmat[] = new double[3];
        double tPlein[] = new double[3];
        double nNodes[] = new double[3];
        double nElmts[] = new double[3];

        int nF = nFile.length;
        for (int k = 0; k < nF; k++) {
            // Time counters and time-variables-templates
            double beg = 0, end = 0;

            ImportFlux mesh = new ImportFlux(file + nFile[k] + ".DEC");

            ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();
            int d = 3;
            int n = ES.getNbNoeud();
            nNodes[k] = n;
            nElmts[k] = ES.getNbElement();

            System.out.println("nbNoeuds= " + ES.getNbNoeud() + "\t nbElmts= " + ES.getNbElement());

            /**
             *********************************************************************
             ***********************************************************************
             * H-MATRIX
             */
            beg = System.currentTimeMillis();
            int ordre = 4;
            double eps = 1e-5;
            HmatrixHCAMagnetoStatPotentialCollocDeg1 f = new HmatrixHCAMagnetoStatPotentialCollocDeg1(new GradNodalDeg1(ES), new NegDotDG(), new SimpleTetraedreGradientPotential(), 15,
                    eps, 50, 60, ordre, 2.0, true);
            double deb = System.nanoTime();
            StorageHmatrix H = new StorageHmatrix(f);
            double fin = System.nanoTime();
            System.out.println("Time to assembly the Hmatrix: " + (fin - deb) / 1e9);

//            H.Coarsen(new TruncationControl("rel", eps));
            end = System.currentTimeMillis();
            tHmat[k] = (end - beg) / 1e3;
            Matrix Mat = null;
            if (k < 2) {
                /*
                 * STOCKAGE PLEIN
                 */
                beg = System.currentTimeMillis();
                CollocationIntegralFormulation IntegralTerm = new CollocationIntegralFormulation(new GradNodalDeg1(ES), new NegDotDG(), new SelfElementFixedGauss(15, new AnalyticalCorrection()));
                IntegralTerm.matAssembly(new GradNodalDeg1(ES));
                Mat = ((StorageFull) IntegralTerm.getStore()).getMatrix();

                end = System.currentTimeMillis();
                tPlein[k] = (end - beg) / 1e3;
                /*
                 * Verification par PMV
                 */
                ColumnVector v = new ColumnVector(n);
                for (int i = 0; i < n; i++) {
                    v.setElement(i, Math.random() * 100);
                }

                ColumnVector exV = new ColumnVector(n);
                exV.mul(Mat, v);

                ColumnVector hV = H.prod(v);

                hV.sub(exV);
                System.out.println("esp= " + eps + "\t Erreur rel par PMV= " + hV.norm() / exV.norm());
            }

        }
    }
}
