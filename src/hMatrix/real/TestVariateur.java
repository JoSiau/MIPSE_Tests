/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshCell;
import g2elab.mipse.meshCore.elements.Element;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.region.Region;
import got.matrix.Matrix;

import java.io.IOException;

/**
 * @author Siau
 */
public class TestVariateur {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

//        ATV71_47931 ; ATV71_80144 ; ATV71_105963
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/VariateurDEC";
        // ImportFlux mesh = new ImportFlux(meshDir+"/ATV71_47931.DEC"); 
        // ImportFlux mesh = new ImportFlux(meshDir+"/ATV71_80144.DEC"); 
        ImportFlux mesh = new ImportFlux(meshDir + "/ATV71_105963.DEC");

        Region[] R = mesh.getRegions();
        int nR = R.length;
        int n[] = new int[nR];
        int nTotal = 0;
        for (int i = 0; i < nR; i++) {
            n[i] = R[i].getElementSet().getNbElement();
            nTotal += n[i];
        }
        Element Elt[] = new Element[nTotal];

        int curN = 0;
        for (int i = 0; i < nR; i++) {
            System.arraycopy(R[i].getElementSet().getElements(), 0, Elt, curN, n[i]);
            curN += n[i];
        }

        Matrix X = new Matrix(nTotal, 3);
        Element E;
        for (int i = 0; i < nTotal; i++)
            X.setRow(i, Elt[i].getCentroid());

        System.out.println("nTotal= " + nTotal);

//        int kmax=50, nmin=30, nbGauss = 3;
//        HmatrixACAGalerkinIntegral f = new HmatrixACAGalerkinIntegral(Elt,new MultG(),new SimpleTriangleConstantChargePotential(),nbGauss,nbGauss);
//        Hmatrix H0 = new Hmatrix(X,f,1e-3,kmax,nTotal,nmin);


        for (int i = 0; i < 27; i++) {
            ExportGmshCell exportVP = new ExportGmshCell(new Cell((ElementSetHomogene) mesh.getRegion(i).getElementSet()), "D:/Var" + i + ".msh");
        }
    }
}
