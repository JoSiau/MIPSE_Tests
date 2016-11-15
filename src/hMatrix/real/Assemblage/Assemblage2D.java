/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Assemblage;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG2D;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author jsiau
 */
public class Assemblage2D {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        System.out.println("Assemblage2D.java");
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/plaque";
        /*
         deb = System.nanoTime();
         // plaque_13532
         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + "/plaque_4548.msh");
         ElementSetHomogene ES = mesh1.createHomogeneESet();
         fin = System.nanoTime();
         System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);
         /*/
//        ImportFlux mesh = new ImportFlux("D:/Meshs/SphereCreuse/SC_2151.DEC");
//        ElementSetHomogene ES = (ElementSetHomogene)mesh.getRegion(0).getElementSet();

        File fi = new File("");
        String path = fi.getAbsolutePath();
//        String file = path + "/src/formulations/g2elab/mipse/formulationInProgress/magnetostatic/potential/PLAQUE2000.DEC";
        String file = path + "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/PLAQUE2000.DEC";
        ImportFlux ImF = new ImportFlux(file);
        Region reg = ImF.getRegion(0);
        ElementSetHomogene ES = (ElementSetHomogene) reg.getElementSet();
        //*/
        int d = 3;
        int n = ES.getNbNoeud();
        System.out.println("n = " + n + "\t d = " + d);

        Cell SC = new Cell(ES);
        Cell PC = SC;
        n = SC.getActiveDofCount();

        int kmax = 50, nmin = 60;
        int order = 4, nbGauss = 7;
        double eps = Math.pow(10, -order - 1);
        int nbGaussFar = nbGauss;
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(PC, SC, new MultG2D(), new SelfElementFixedGauss(nbGauss, new InCreasedPGSourceNumber(16)), nbGauss, nbGaussFar,
                eps, kmax, nmin, order);
        f.assembly();
        StorageHmatrix H = (StorageHmatrix) f.getStore();

        GalerkinIntegralFormulationFull IF = new GalerkinIntegralFormulationFull(PC, SC, new MultG2D(), new SelfElementFixedGauss(nbGauss, new InCreasedPGSourceNumber(16)), nbGauss);
        IF.assembly();
        Matrix M = ((StorageFull) IF.getStore()).getMatrix();

        // SECOND MEMBRE
        ColumnVector v = new ColumnVector(n);
        for (int i = 0; i < n; i++) {
            v.setElement(i, Math.random() * 1000);
        }
        // LE VECTEUR EXACTE
        ColumnVector ex = new ColumnVector(n);
        ex.mul(M, v);
        // LE VECTEUR APPROCHE
        ColumnVector vH = H.prod(v);
        // CALCUL DE LERREUR
        vH.sub(ex);
        System.err.println("Erreur relative= " + vH.norm() / ex.norm());

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

        //Affiche la Hmatrice
        H.printOnJFrame();
    }
    
}
