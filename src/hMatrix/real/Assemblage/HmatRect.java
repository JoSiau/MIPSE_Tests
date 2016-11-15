/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package hMatrix.real.Assemblage;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.contraints.ExternalDriveFaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceRealDirichlet;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

/**
 * @author jsiau
 */
public class HmatRect {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.out.println("HmatRect.java");
        double deb, fin;
        deb = System.nanoTime();
        ImportFlux mesh1 = new ImportFlux("D:/Meshs/poutre/POUTRE.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        fin = System.nanoTime();
        System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);

        FaceDeg1 B = new FaceDeg1(ES, new FaceRealDirichlet(mesh1.getRegion(3), 0));

        int d = 3;
        int n = B.getActiveDofCount();//ES.getNbElement();
        System.out.println(" n = " + n + "\t d = " + d);

        FaceDeg1 C = new FaceDeg1(ES, new ExternalDriveFaceConstraint(mesh1.getRegion(3)));
        int kmax = 50, nmin = 30;
        int order = 3;
        double eps = Math.pow(10, -order - 1);
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, B, new MultGvect(), new SelfElementFixedGauss(15, new Cancel()), 15, 15,
                eps, kmax, nmin, order, 2.0, true);
        f.assembly();
        StorageHmatrix H = (StorageHmatrix) f.getStore();

//        H.Coarsen(new TruncationControl("rel",eps));
        GalerkinIntegralFormulationFull IF = new GalerkinIntegralFormulationFull(C, B, new MultGvect(), new SelfElementFixedGauss(15, new Cancel()), 15);
        IF.assembly();
        Matrix M = ((StorageFull) IF.getStore()).getMatrix();
        System.out.println("M.size= " + M.getRowCount() + " x " + M.getColumnCount());
        ColumnVector v = new ColumnVector(n);
        for (int i = 0; i < n; i++) {
            v.setElement(i, Math.random() * 1000);
        }
        v.setAllElements(1);
        ColumnVector ex = new ColumnVector(M.getRowCount());
        ex.mul(M, v);

        ColumnVector vH = H.prod(v);

        vH.sub(ex);
        System.out.println("Erreur relative= " + vH.norm() / ex.norm());

        H.printOnJFrame();
    }

}
