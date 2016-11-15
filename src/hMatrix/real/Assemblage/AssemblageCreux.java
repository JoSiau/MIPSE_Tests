/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Assemblage;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMesh;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.mipseCore.integralIntegration.kernel.KernelInterface;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.BinaryTree.BlockClusterTree;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulation;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulationHmatrix;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.mipseCore.storage.StorageSparse;
import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;

import java.io.IOException;
import java.util.Arrays;

/**
 * @author jsiau
 */
public class AssemblageCreux {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(1);
        System.out.println("TestFiniteElementAssembly.java");
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/";

        deb = System.nanoTime();
        //*
        ImportGmshMesh mesh1 = new ImportGmshMesh(meshDir + "sphere/SPHERE_4893.msh");
        mesh1.meshSummary();
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
         /*/
        ImportFlux mesh1 = new ImportFlux(meshDir + "sphere/SPHERE_2104.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        //*/
        fin = System.nanoTime();
        System.out.println("Time to load Gmsh mesh: " + (fin - deb) / 1e9);

        //*
        FaceDeg1 C = new FaceDeg1(ES);
        KernelInterface k = new MultGvect();
        /*/
         Cell C = new Cell(ES);
         KernelInterface k = new MultG();
         //*/
        int d = 3;
        int n = C.getActiveDofCount();
        System.out.println(" n = " + n + "\t d = " + d);

        int nbGauss = 1;
        ////////////////////////////////////////////////////
        ////////////////////////////////////////////////////
        FiniteElementFormulation IF = new FiniteElementFormulation(C, C);
        long t = System.currentTimeMillis();
        IF.assembly(nbGauss);
        long t2 = System.currentTimeMillis();
        System.out.println("Time to assemble FE matrix = " + (t2 - t) + " msec.");
        SparseMatrixRowReal s = ((StorageSparse) IF.getStore()).getSparseMat();
        System.out.println("S.Mem = " + s.getMemoryUsed());
        ////////////////////////////////////////////////////
        BlockClusterTree bct = new BlockClusterTree(C, s, 3, 50, 2.0);
        StorageHmatrix hs1 = new StorageHmatrix(bct);
        hs1.printOnJFrame("1");
        ////////////////////////////////////////////////////
        FiniteElementFormulation Ifh = new FiniteElementFormulationHmatrix(C, 50, 2.0);
        Ifh.assembly(nbGauss);
        StorageHmatrix hs2 = (StorageHmatrix) Ifh.getStore();
        hs2.printOnJFrame("2");
        ////////////////////////////////////////////////////
        ////////////////////////////////////////////////////
        double x[] = new double[C.getActiveDofCount()];
        Arrays.fill(x, 1.0);
        double r_ex[] = s.product(x, new double[n]);
        //
        double r1[] = hs1.product(x, new double[n]);
        //
        double r2[] = hs2.product(x, new double[n]);
        //
        ColumnVector v_ex = new ColumnVector(r_ex);
        ColumnVector v1 = new ColumnVector(r1);
        ColumnVector v2 = new ColumnVector(r2);
        //
        v1.sub(v_ex);
        v2.sub(v_ex);
        double norm = v_ex.norm();
        System.out.println("Erreur v1 = " + v1.norm() / norm);
        System.out.println("Erreur v2 = " + v2.norm() / norm);

        v1.sub(v2);
        System.out.println("Diff 1-2 = " + v1.norm());

    }

}
