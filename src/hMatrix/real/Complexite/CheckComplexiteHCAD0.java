/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Complexite;

import g2elab.mipse.meshCore.IO.gmsh.ExportGmsh;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.Preconditionners.HmatrixLUDecomposition;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import got.matrix.RowVector;

import java.io.IOException;

/**
 * @author jsiau
 */
public class CheckComplexiteHCAD0 {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/Sphere_gmsh/sphere_";

        int nFile[] = new int[5];
        nFile[0] = 6481;
        nFile[1] = 12529;
        nFile[2] = 34449;
        nFile[3] = 62713;
        nFile[4] = 132241;
//        nFile[5] = 301791;

        int nF = nFile.length;

        double tHmat[] = new double[nF];
        double nNodes[] = new double[nF];
        double nElmts[] = new double[nF];
        double tpsResHm[] = new double[nF];
        double stor[] = new double[nF];

        for (int k = 0; k < nF; k++) {
            System.out.println("=== HCA ===");
            // Time counters and time-variables-templates
            double beg = 0, end = 0;

            ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + nFile[k] + ".msh");
            ElementSetHomogene ES = mesh1.createHomogeneESet();

            int d = 3;
            int n = ES.getNbElement();
            nElmts[k] = ES.getNbElement();

            System.out.println("nbNoeuds= " + ES.getNbNoeud() + "\t nbElmts= " + ES.getNbElement());
            /*
             * DEFINIE LA FONCTION (LE NOYAU) QUE L'ON DESIRE
             */
            Cell alpha = new Cell(ES);
            deb = System.nanoTime();
            GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(alpha, alpha, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3, 3,
                    1e-4, 50, 30, 3, 2.0, false);
            System.out.println("Calculs de quadratures = " + (System.nanoTime() - deb) * 1e-9);
            f.assembly();
            StorageHmatrix Mh = (StorageHmatrix) f.getStore();
//            Mh.Agglomerate(new TruncationControl("rel", 1e-4));
            fin = System.nanoTime();
            tHmat[k] = (fin - deb) / 1e9;
            System.err.println("Time to assembly the Hmatrix = " + tHmat[k]);
            stor[k] = Mh.getStorageSizePerDofs();
            //*
            deb = System.nanoTime();
            TruncationControl tolPrecond = new TruncationControl("rel", 1e-1);
            StorageHmatrix precondH = Mh.copy(true);
            precondH.Agglomerate(tolPrecond);
            HmatrixLUDecomposition hLu = new HmatrixLUDecomposition(precondH, tolPrecond);
            fin = System.nanoTime();
            System.err.println("Time to compute the Hmatrix Preconditionner = " + (fin - deb) / 1e9);
            stor[k] += (hLu.getMemoryUsed() / 1024.0) / (double) n;
             /*/
            PrecondIdentityReal hLu = new PrecondIdentityReal();
            //*/

            RowVector valueNd = new RowVector(n);
            valueNd.setAllElements(1);
            RealScalarCellQuantity V0p1 = new RealScalarCellQuantity(alpha, valueNd);
            RowVector bg1 = (RowVector) V0p1.projWithExplicitDof(alpha, 3);
            double[] secondMembre = new double[n];
            bg1.get(secondMembre);

            deb = System.nanoTime();
            FGMResReal solvMH = new FGMResReal(Mh, hLu);
            double configRes[] = {100, 1e-8, 1, -1};
            solvMH.setInfoResolution(configRes);
            double x[] = solvMH.solve(new double[n], secondMembre);
            fin = System.nanoTime();
            tpsResHm[k] = (fin - deb) * 1e-9;
            System.err.println("Time to solve the Hmatrix = " + tpsResHm[k]);


            ExportGmsh expgmsh1 = new ExportGmsh(ES, "D:/HCA_phi.msh");
            RealScalarCellQuantity quant = new RealScalarCellQuantity(alpha, new RowVector(x));
            quant.exportGmsh(expgmsh1, "charge");

        }
        System.out.println("\n\n");
        for (int i = 0; i < nF; i++) {
            System.out.println(nElmts[i] + " , " + tHmat[i] + " , " + tpsResHm[i] + " , " + stor[i]);
        }
    }

}

/*
 (6198.0 , 20.500353)
 (12138.0 , 42.137424)
 (33806.0 , 136.41139)
 (61842.0 , 287.13654)
 (130974.0 , 662.65658)


 Time to assembly the Hmatrix = 38.035200382
 Time to assembly the Hmatrix = 74.893098273
 Time to assembly the Hmatrix = 237.61820245
 Time to assembly the Hmatrix = 540.547375859
 Time to assembly the Hmatrix = 1192.967740139
 */
