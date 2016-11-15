/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Complexite;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.octree.repartitionElements.RepartitionElemNbMaxElem;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFMM;
import g2elab.mipse.mipseCore.storage.StorageFMMGalerkine;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import g2elab.mipse.numericalTools.preconditioner.PrecondIdentityReal;
import got.matrix.RowVector;

import java.io.IOException;

/**
 * @author jsiau
 */
public class CheckComplexiteFMM {

    public static void main(String[] args) throws IOException {
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().setNbTaches(1);
        double deb, fin;
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/Sphere_gmsh/sphere_";

        int nFile[] = new int[6];
        nFile[0] = 6481;
        nFile[1] = 12529;
        nFile[2] = 34449;
        nFile[3] = 62713;
        nFile[4] = 132241;
        nFile[5] = 301791;

        int nF = nFile.length;

        double tFMM[] = new double[nF];
        double nNodes[] = new double[nF];
        double nElmts[] = new double[nF];
        double tpsResFMM[] = new double[nF];
        double stor[] = new double[nF];

        for (int k = 0; k < nF - 1; k++) {
            ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + nFile[k] + ".msh");
            ElementSetHomogene ES = mesh1.createHomogeneESet();

            int d = 3;
            int n = ES.getNbElement();
            nElmts[k] = ES.getNbElement();

            System.out.println("nbNoeuds= " + ES.getNbNoeud() + "\t nbElmts= " + ES.getNbElement());
            /*
             * DEFINIE LA FONCTION (LE NOYAU) QUE L'ON DESIRE
             */
            Cell C = new Cell(ES);
            deb = System.nanoTime();
            GalerkinIntegralFormulationFMM f = new GalerkinIntegralFormulationFMM(C, C, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3, 3,
                    new RepartitionElemNbMaxElem(30), 1.0, 3, 1, 1);
            f.assembly();
            StorageFMMGalerkine FMM = (StorageFMMGalerkine) f.getStore();
            fin = System.nanoTime();
            tFMM[k] = (fin - deb) / 1e9;
            System.err.println("Time to compute the FMM = " + tFMM[k]);

            stor[k] = (FMM.getMemoryUsed() / 1024.0) / (double) n;
            System.out.println("stor= " + stor[k]);
            /*
             deb = System.nanoTime();
             TruncationControl tolPrecond = new TruncationControl("rel", 3e-1);
             StorageHmatrix precondH = Mh.copy(true);
             precondH.Agglomerate(tolPrecond);
             HmatrixLUDecomposition hLu = new HmatrixLUDecomposition(precondH, tolPrecond);
             fin = System.nanoTime();
             System.err.println("Time to compute the Hmatrix Preconditionner = " + (fin - deb) / 1e9);
             storHm[k] += hLu.getMemoryUsed();
             /*/
            PrecondIdentityReal hLu = new PrecondIdentityReal();
            //*/

            RowVector valueNd = new RowVector(n);
            valueNd.setAllElements(1);
            RealScalarCellQuantity V0p1 = new RealScalarCellQuantity(C, valueNd);
            RowVector bg1 = (RowVector) V0p1.projWithExplicitDof(C, 3);
            double[] secondMembre = new double[n];
            bg1.get(secondMembre);

            deb = System.nanoTime();
            FGMResReal solvMH = new FGMResReal(FMM, hLu);
            double configRes[] = {100, 1e-8, 1, -1};
            solvMH.setInfoResolution(configRes);
            double x[] = solvMH.solve(new double[n], secondMembre);
            fin = System.nanoTime();
            tpsResFMM[k] = (fin - deb) * 1e-9;
            System.err.println("Time to solve the FMM = " + tpsResFMM[k]);
        }
        System.out.println("\n\n");
        for (int i = 0; i < nF; i++) {
            System.out.println(nElmts[i] + " , " + tFMM[i] + " , " + tpsResFMM[i] + " , " + stor[i]);
        }
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();

    }

}
/*
7.318401196
19.1451293550
100.027476147
266.332834117

// RepartitionNiveauConstant(7)
45.095827268
128.972560626
579.88962538
1052.770634758

6198.0 , 3101.0 , 42.118583157
12138.0 , 6071.0 , 116.333366104
33806.0 , 16905.0 , 531.767729789
61842.0 , 30923.0 , 1007.325898467
130974.0 , 65489.0 , 1675.139043498
*/
