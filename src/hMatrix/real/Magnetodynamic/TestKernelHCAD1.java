/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Magnetodynamic;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.functionSpace.CurlSEdgeDeg1;
import g2elab.mipse.meshCore.functionSpace.GradNodalDeg1;
import g2elab.mipse.meshCore.functionSpace.Hgrad;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.CrossDGmulAlpha;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;

/**
 * COMPARE L'HCA SUR LA FORMULATION AVEC LE STOCKAGE PLEIN
 *
 * @author jsiau
 */
public class TestKernelHCAD1 {

    /**
     * @param args the command line arguments
     * @throws java.lang.Exception
     */
    public static void main(String[] args) throws Exception {
        String className = "";
        if (args.length != 0) {
            className = args[0];
        }

        File f = new File("");
        String path = f.getAbsolutePath();
        /*
         //        String file = path + "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/PLAQUE2000.DEC";
//         String file = path + "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere/SPHERE_1884.DEC";
         String file = path + "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere/SPHERE_1884.DEC";
         ImportFlux ImF = new ImportFlux(file);
         Region reg = ImF.getRegions(0);
         ElementSurfSetHomogene mesh = (ElementSurfSetHomogene) reg.getElementSet();
         /*/
        ImportGmshMeshRegion ImF = new ImportGmshMeshRegion(0, "D:/Meshs/plaqueInclinee/plaqueIncline_2400e.msh");
//        ImportGmshMeshRegion ImF = new ImportGmshMeshRegion(0,"D:/Meshs/plaqueNonCoPlanaire/plaqueNCP_2800e.msh");
        ElementSurfSetHomogene mesh = (ElementSurfSetHomogene) ImF.createHomogeneESet();
        //*/

        System.out.println("Nombre elements " + mesh.getNbElement());
        System.out.println("Nombre de noeuds " + mesh.getNbNoeud());

//        InactiveBorderNode Dir0border = new InactiveBorderNode();
        NodalDeg1 alpha = new NodalDeg1(mesh);//, Dir0border);
        Hgrad alphaN = alpha.createProjOnNormal(1);

        GradNodalDeg1 gradalpha = new GradNodalDeg1(mesh);//, Dir0border);
        CurlSEdgeDeg1 rotSalpha = new CurlSEdgeDeg1(gradalpha);

        int n = mesh.getNbNoeud();
        ColumnVector v = new ColumnVector(n);
//        v.setAllElements(1);
//        v.setElement(0, 0);
//        v.setElement(1, 0);

        for (int i = 0; i < n; i++) {
            v.setElement(i, Math.random() * 10000);
        }

        int nbGauss = 7;
        /*
         COMPUTE THE FULL STORAGE
         */
        GalerkinIntegralFormulation IF = new GalerkinIntegralFormulationFull(alphaN, rotSalpha, new CrossDGmulAlpha(1), new SelfElementFixedGauss(3, new Cancel()), nbGauss);
        IF.assembly();
        Matrix M = ((StorageFull) IF.getStore()).getMatrix();
        // Sauvegarde la matrice pour eviter de la recalculer a chaque test !
//        Ecriture ec = new Ecriture("D:/M" + n + ".out");
//        double Mt[] = new double[n * n];
//        M.get(Mt);
//        ec.ecrire(Mt, ';');
//        ec.close();
        /*/
         Lecture lec = new Lecture("D:/M"+n+".out");
         double Mt[] = lec.lireLigneDoubleVecteur(';');
         Matrix M = new Matrix(n,n,Mt);
         //*/
        ColumnVector vF = new ColumnVector(n);
        vF.mul(M, v);
        /*
         COMPUTE THE HMATRIX WITH THE HCA
         */
        double eps = 1e-5;
        int kmax = 50, nmin = 30, order = 4;
        GalerkinIntegralFormulationHCA IF2 = new GalerkinIntegralFormulationHCA(alphaN, rotSalpha, new CrossDGmulAlpha(/*frequence * mu0 * epaisseur*/1), new SelfElementFixedGauss(3, new Cancel()), nbGauss, nbGauss,
                eps, kmax, nmin, order, 2.0, true);
        IF2.assembly();

        StorageHmatrix H = (StorageHmatrix) IF2.getStore();
        ColumnVector vH0 = new ColumnVector(n);
        vH0 = H.prod(v);

        ColumnVector vH = vH0.copy();
        vH.sub(vF);
        double errRel = vH.norm() / vF.norm();
        System.err.println("ERREUR RELATIVE = " + errRel);
        System.out.println("vH.norm= " + vH0.norm() + "\t vF.norm= " + vF.norm() + "\t Erreur absolue= " + vH.norm());

        // Check if the integral are wrong
        if (H.getRoot().isThereAnyInfinity()) {
            System.err.println("H HAS SOME INFINITY !!");
        }

        if (H.getRoot().isThereAnyNaN()) {
            System.err.println("H HAS SOME NaN !!");
        }

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

//        new printHmatrix(H);
    }

}
/*
 Sphere1884
 1e-3 , 2 , 2.3774528061524886E-5
 1e-4 , 3 , 4.182812750611399E-6
 1e-5 , 4 , 5.240563184846922E-7
 1e-6 , 5 , 1.6255123889614657E-7
 1e-7 , 6 , 1.1425800575640523E-7
 Plaque2000
 1e-3 , 2 , 6.805856043911309E-5
 1e-4 , 3 , 2.9028593733613067E-5
 1e-5 , 4 , 8.285774074106398E-7
 1e-6 , 5 , 5.758994593927439E-8
 1e-7 , 6 , 4.7732238531554634E-9
 */
