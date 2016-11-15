/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Magnetodynamic;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.functionSpace.CurlSEdgeDeg1;
import g2elab.mipse.meshCore.functionSpace.GradNodalDeg1;
import g2elab.mipse.meshCore.functionSpace.Hgrad;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.NegCrossDG;
import g2elab.mipse.mipseCore.integralIntegration.kernel.NegDotDG;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;

/**
 * @author jsiau
 */
public class TestDebugg {

    /**
     * SUB SPACES
     */
    NodalDeg1 alpha;
    Hgrad alphaN;
    GradNodalDeg1 gradalpha;
    Hgrad rotSalpha;
    int nDof;
    int nbGauss = 7;

    public TestDebugg(ElementSetHomogene mesh) {
        nDof = mesh.getNbNoeud();

        alpha = new NodalDeg1(mesh);
        alphaN = alpha.createProjOnNormal(1);

        gradalpha = new GradNodalDeg1(mesh);
        rotSalpha = new CurlSEdgeDeg1(gradalpha);

        nbGauss = 7;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        File f = new File("");
        String path = f.getAbsolutePath();
//        String file = path + "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/PLAQUE2000.DEC";
        String file = path + "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere/SPHERE_1884.DEC";
        ImportFlux ImF = new ImportFlux(file);
        Region reg = ImF.getRegion(0);
        ElementSurfSetHomogene mesh = (ElementSurfSetHomogene) reg.getElementSet();

        System.out.println("Nombre elements " + mesh.getNbElement());
        System.out.println("Nombre de noeuds " + mesh.getNbNoeud());

        TestDebugg t = new TestDebugg(mesh);

        /* Test si les stockages pleins des 2 formulations sont identiques.
         System.out.println("CALCUL DES STOCKAGES PLEINS");
         t.testFullStorage();
         /*/
        //*
        System.out.println("CALCUL DES STOCKAGES HMATRIX");
        t.testHmatStorage();
         /*/
        System.out.println("CALCULS DES H-MATRICES / FORMULATIONS");
        t.testHmats();
        //*/
    }

    /**
     * TEST LES STOCKAGES HMATRICES DES 2 FORMULATIONS
     */
    public void testHmats() {
        ColumnVector v = new ColumnVector(nDof);
        for (int i = 0; i < nDof; i++) {
            v.setElement(i, Math.random() * 1000);
        }
        // HMATRICES
        double eps = 1e-6;
        int kmax = 50, nmin = 30, order = 5;
        System.out.println("CALCUL DE LA HMATRIX 1");
        GalerkinIntegralFormulationHCA IF1 = new GalerkinIntegralFormulationHCA(alphaN, rotSalpha, new NegCrossDG(), new SelfElementFixedGauss(3, new Cancel()), nbGauss, nbGauss,
                eps, kmax, nmin, order, 2.0);
        IF1.assembly();
        StorageHmatrix H1 = (StorageHmatrix) IF1.getStore();
        int[] ind1 = H1.getDof2Idx();
        ColumnVector vH1 = H1.prod(v);

        System.out.println("CALCUL DE LA HMATRIX 2");
        GalerkinIntegralFormulationHCA IF2 = new GalerkinIntegralFormulationHCA(alpha, gradalpha, new NegDotDG(), new SelfElementFixedGauss(1, new Cancel()), nbGauss, nbGauss, eps, kmax, nmin, order, 2.0);
        IF2.assembly();
        StorageHmatrix H2 = (StorageHmatrix) IF2.getStore();
        int[] ind2 = H2.getDof2Idx();
        ColumnVector vH2 = H2.prod(v);

        vH1.sub(vH2);
        System.out.println("Erreur relative= " + vH1.norm() / vH2.norm());

    }


    /**
     * TEST DE LA NOUVELLE FORMULATION AVEC STOCKAGE PLEIN ET HMAT
     */
    public void testHmatStorage() {
        ColumnVector v = new ColumnVector(nDof);
        for (int i = 0; i < nDof; i++) {
            v.setElement(i, Math.random() * 1000);
        }
        // HMATRICES
        double eps = 1e-6;
        int kmax = 50, nmin = 30, order = 5;
        System.out.println("CALCUL DE LA HMATRIX");
        GalerkinIntegralFormulationHCA IF2 = new GalerkinIntegralFormulationHCA(alphaN, rotSalpha, new NegCrossDG(), new SelfElementFixedGauss(nbGauss, new Cancel()), nbGauss, nbGauss,
                eps, kmax, nmin, order, 2.0);
        IF2.assembly();
        StorageHmatrix H = (StorageHmatrix) IF2.getStore();
        ColumnVector vH = H.prod(v);
        // Check if the integral are wrong
        if (H.getRoot().isThereAnyInfinity()) {
            System.err.println("H HAS SOME INFINITY !!");
        }

        if (H.getRoot().isThereAnyNaN()) {
            System.err.println("H HAS SOME NaN !!");
        }

        // COMPUTE THE APPROXIMATION
        System.out.println("COMPUTE LAPPROXIMATION");
        Matrix Mh = H.getApprox();
        int in[] = H.getDof2Idx();
        Mh.permuteColumns(in, new int[nDof]);
        Mh.permuteRows(in, new int[nDof]);

        System.out.println("CALCUL INTEGRAL 1");
        GalerkinIntegralFormulation IF = new GalerkinIntegralFormulationFull(alphaN, rotSalpha, new NegCrossDG(), new SelfElementFixedGauss(3, new Cancel()), nbGauss);
        IF.assembly();
        Matrix M = ((StorageFull) IF.getStore()).getMatrix();
        ColumnVector vM = new ColumnVector(nDof);
        vM.mul(M, v);

        System.out.println("CALCUL DE LERREUR");
        Mh.sub(M);
        double Mhnorm = Mh.norm();
        System.out.println("ERREUR ABSOLUE= " + Mhnorm + "\t ERREUR RELATIVE= " + Mhnorm / M.norm());

        vH.sub(vM);
        System.out.println("erreur relative PMV= " + vH.norm() / vM.norm());
    }

    /**
     * TEST LES 2 FORMULATIONS EN PLEIN !
     */
    public void testFullStorage() {
        System.out.println("CALCUL INTEGRAL 1");
        GalerkinIntegralFormulation IF = new GalerkinIntegralFormulationFull(alphaN, rotSalpha, new NegCrossDG(), new SelfElementFixedGauss(3, new Cancel()), nbGauss);
        IF.assembly();
        Matrix M = ((StorageFull) IF.getStore()).getMatrix();

        System.out.println("CALCUL INTEGRAL 2");
        GalerkinIntegralFormulationFull IFanto = new GalerkinIntegralFormulationFull(alpha, gradalpha, new NegDotDG(), new SelfElementFixedGauss(nbGauss, new Cancel()), nbGauss);
        IFanto.assembly();
        Matrix Ma = ((StorageFull) IFanto.getStore()).getMatrix();

        System.out.println("CALCUL ERREUR");
        M.sub(Ma);
        double Mnorm = M.norm();
        System.out.println("Erreur absolue= " + Mnorm + "\t Erreur relative= " + (Mnorm / Ma.norm()));
        // Erreur absolue= 9.099984703929158E-17	 Erreur relative= 1.05053027798521E-16
    }

}
