/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Magnetodynamic.Jfacet;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.contraints.FaceRealDirichlet;
import g2elab.mipse.meshCore.elements.ElementVolSetHomogene;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;

/**
 * Test l'assemblage avec des elements de facettes.
 *
 * @author jsiau
 */
public class Assemblage {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        System.out.println("Assemblage.java");
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(2);
        /*
        /**/
        String path = new File("./../FORMULATIONS/").getAbsolutePath();
        String nom = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/PLATE_1000_TED.DEC";
        /*/
//        String nom = "D:/Meshs/Problem7/PLATE_30k.DEC";
//        String nom = "D:/Meshs/Problem7/mesh_oringal/PLATE2.DEC";
        String nom = "D:/Meshs/Problem7/mesh_oringal/AidedMesh/PLATE2.DEC";
        //*/
        ImportFlux ImF = new ImportFlux(nom);
        Region reg = ImF.getRegion(0);
        ElementVolSetHomogene mesh = new ElementVolSetHomogene(reg.getElementSet().getElements());
//        ElementVolSetHomogene mesh = (ElementVolSetHomogene) reg.getElementSet();

        //* Test avec constraintes
        FaceRealDirichlet FluxNull = new FaceRealDirichlet(ImF.getRegion(1), 0);
        FaceDeg1 Wi = new FaceDeg1(mesh, FluxNull);
        /*/ // Test sans contraintes
        FaceDeg1 Wi = new FaceDeg1(mesh);//,constraints);
        //*/

        System.out.println("nbFacets= " + Wi.getFacesSet().getNbFaces() + "\t nbActiveDofCount= " + Wi.getActiveDofCount());

        /*
//        Test le decoupage de lespace
        BlockClusterTree bct = new BlockClusterTree(Wi, 3, 30, 2.0);
        bct.inspectLeaves();
        return;
        /*/
        int n = Wi.getActiveDofCount();
        ColumnVector v = new ColumnVector(n);
        for (int i = 0; i < n; i++) {
            v.setElement(i, Math.random() * 10000);
        }

        double eps = 1e-3;
        int kmax = 50, nmin = 50, order = 2, nbGauss = 15, nbGaussFar = nbGauss;
        GalerkinIntegralFormulationHCA IF2 = new GalerkinIntegralFormulationHCA(Wi, Wi, new MultGvect(), new SelfElementFixedGauss(nbGauss, new Cancel()), nbGauss, nbGaussFar,
                eps, kmax, nmin, order);
        IF2.assembly();


        //* Sauvegarde la matrice pour eviter de la recalculer a chaque test !
        double deb = System.nanoTime();
        GalerkinIntegralFormulationFull IF1 = new GalerkinIntegralFormulationFull(Wi, Wi, new MultGvect(), new SelfElementFixedGauss(nbGauss, new Cancel()), nbGauss);
        IF1.assembly();
        double fin = System.nanoTime();
        System.out.println("Total time to compute the Full matrix= " + (fin - deb) * 1e-9);
        Matrix M = ((StorageFull) IF1.getStore()).getMatrix();
        Ecriture ec = new Ecriture("D:/M" + n + ".out");
        double Mt[] = new double[n * n];
        M.get(Mt);
        ec.ecrire(Mt, ';');
        ec.close();
        /*/
         Lecture lec = new Lecture("D:/M"+n+".out");
         double Mt[] = lec.lireLigneDoubleVecteur(';');
         Matrix M = new Matrix(n,n,Mt);
         //*/
        ColumnVector vF = new ColumnVector(n);
        vF.mul(M, v);

        StorageHmatrix H = (StorageHmatrix) IF2.getStore();
        ColumnVector vH0 = H.prod(v);

        ColumnVector vH = vH0.copy();
        vH.sub(vF);
        double errRel = vH.norm() / vF.norm();
        System.err.println("ERREUR RELATIVE par pmv= " + errRel);
        
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

        ExportGmshHdiv exportTm = new ExportGmshHdiv(Wi, "D:/meshFacet.msh");
        exportTm.saveMesh();
    }

}
