/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real;

import g2elab.mipse.meshCore.IO.gmsh.ExportGmsh;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.matrixCompression.octree.repartitionElements.RepartitionElemNbMaxElem;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFMM;
import g2elab.mipse.mipseCore.storage.Storage;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import g2elab.mipse.numericalTools.preconditioner.PrecondIdentityReal;
import got.matrix.RowVector;

/**
 * @author jsiau
 */
public class FmmVsHmAuto {

    int n;

    public static void main(String[] args) {
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().setNbTaches(1);
        double deb, fin;
        // Sphere Creuse
        // SC_3846 SC_6955 SC_12938 SC_26420 SC_48154 SC_79380
        String file = "D:/Meshs/Sphere_gmsh/sphere_";

        int nFile[] = new int[5];
//        nFile[0] = 31696;
        nFile[0] = 6481;
        nFile[1] = 12529;
        nFile[2] = 34449;
        nFile[3] = 62713;
        nFile[4] = 132241;
//        nFile[5] = 301791;

        int nF = nFile.length;

        double nNodes[] = new double[nF];
        double nElmts[] = new double[nF];
        // total time
        double tpsTotalFMM[] = new double[nF];
        double tpsTotalHm[] = new double[nF];
        // Time to assemble
        double tpsAssFMM[] = new double[nF];
        double tpsAssHm[] = new double[nF];
        // Time to compute the resolution
        double tpsResFMM[] = new double[nF];
        double tpsResHm[] = new double[nF];
        // The storage
        double storFMM[] = new double[nF];
        double storHm[] = new double[nF];
        // Redundun data
        Cell C = null;
        double x[];
        double secondMembre[] = null;
        GalerkinIntegralFormulation f;
        for (int k = 0; k <= nF - 1; k++) {
            ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, file + nFile[k] + ".msh");
            ElementSetHomogene ES = mesh1.createHomogeneESetWithRenumbering();

            int d = 3;
            int n = ES.getNbElement();
            nElmts[k] = ES.getNbElement();

            secondMembre = new double[n];
            for (int i = 0; i < n; i++) {
                secondMembre[i] = 1;
            }
            System.out.println("nbElmts= " + ES.getNbElement());
            /*
             * DEFINIE LA FONCTION (LE NOYAU) QUE L'ON DESIRE
             */
            C = new Cell(ES);
            System.out.println("RESOLUTION FMM");
            /*$$
            
             * FMM
             */
            double eps = 1e-2;
            double debFMM = System.nanoTime();
            deb = System.nanoTime();
            f = new GalerkinIntegralFormulationFMM(C, C, new MultG(), new SelfElementFixedGauss(4, new AnalyticalCorrection()), 4, 4,
                    new RepartitionElemNbMaxElem(30), 1.0, 3, 1, 1);
//             f = new GalerkinIntegralFormulationHCA(C, C, new MultG(), new SelfElementFixedGauss(4),new AnalyticalCorrection(), 4, 4,
//                    eps, 50, 30, Math.abs((int) Math.floor(Math.log10(eps))), 2.0, true);
//            f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(4,new AnalyticalCorrection()),
//                     4, eps, 50, 30, 2.0, false);
            // comparer adaptatif et niveau constant
            // verifier la complexite du nivo constant
            // comparer les options memoire 
            f.assembly();
            fin = System.nanoTime();
            tpsAssFMM[k] = (fin - deb) / 1e9;
            System.err.println("Time to compute the FMM = " + tpsAssFMM[k]);
            Storage Mf = f.getStore();
            storFMM[k] = Mf.getMemoryUsed();
            System.out.println("Storage FMM = " + storFMM[k]);

            /*
             * RESOLUTION
             */
            // THE PRECOND IS THE IDENTITY
            PrecondIdentityReal Id = new PrecondIdentityReal();
            deb = System.nanoTime();
            FGMResReal solvFMM = new FGMResReal(Mf, Id);
            /* Parametres du solveurFGMRes
             1) Nombre d'iterations maximal du solveur
             2) Critère d'arrêt (Norme du résidu)
             3) Affichage de l'évolution du solveur
             0 -> Pas d'affichage
             1 -> Affichage
             4) Taille du sous-espace de Krylov (plus petit que le nombre d'inconnues du système)
             */
            double configRes[] = {100, 1e-8, 1, -1};
            solvFMM.setInfoResolution(configRes);
            x = solvFMM.solve(new double[n], secondMembre);
            fin = System.nanoTime();
            tpsResFMM[k] = (fin - deb) / 1e9;
            System.err.println("Time to solve the FMM = " + tpsResFMM[k]);
            tpsTotalFMM[k] = (fin - debFMM) * 1e-9;
            System.out.println("Total time FMM = " + tpsTotalFMM[k]);

            eps = 1e-2;// = 1%
            ExportGmsh expgmsh1 = new ExportGmsh(ES, "D:/ComparaisonMHFMM2.msh");
            RealScalarCellQuantity quant = new RealScalarCellQuantity(C, new RowVector(x));

            quant.exportGmsh(expgmsh1, "chargefmm");
            System.out.println("RESOLUTION HMATRIX");
            /*
             * HCA 
             */
            double debHM = System.nanoTime();
            deb = System.nanoTime();
//            f = new GalerkinIntegralFormulationHCA(C, C, new MultG(), new SelfElementFixedGauss(4),new AnalyticalCorrection(), 4, 4,
//                    eps, 50, 30, Math.abs((int) Math.floor(Math.log10(eps))), 2.0, true);
            f = new GalerkinIntegralFormulationACA(C, C, new MultG(), new SelfElementFixedGauss(4, new AnalyticalCorrection()),
                    4, eps, 50, 30);
            f.assembly();
            StorageHmatrix Mh = (StorageHmatrix) f.getStore();
            Mh.Agglomerate(new TruncationControl("rel", eps));
            fin = System.nanoTime();
            tpsAssHm[k] = (fin - deb) / 1e9;
            System.err.println("Time to compute the Hmatrix = " + tpsAssHm[k]);
            storHm[k] = Mh.getMemoryUsed();
            /*
             * Preconditionneur
             */
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
        /*
             * RESOLUTION
             */
            deb = System.nanoTime();
            FGMResReal solvMH = new FGMResReal(Mh, hLu);
            solvMH.setInfoResolution(configRes);
            x = solvMH.solve(new double[n], secondMembre);
            fin = System.nanoTime();
            tpsResHm[k] = (fin - deb) * 1e-9;
            System.err.println("Time to solve the Hmatrix = " + tpsResHm[k]);
            tpsTotalHm[k] = (fin - debHM) * 1e-9;
            RealScalarCellQuantity quant2 = new RealScalarCellQuantity(C, new RowVector(x));

            quant2.exportGmsh(expgmsh1, "chargeHmat");
            System.out.println("Total  \tFMM\tHmatrix");
            System.out.format("Time   \t\t%.3f\t%.3f", tpsTotalFMM[k], tpsTotalHm[k]);
            System.out.println("\n Storage FMM= " + storFMM[k] + "\t Storage HMatrix= " + storHm[k]);
        }
        System.out.println("\n\n");
        System.out.println("nElmt \t tpsAssFMM \t tpsResFMM \t storFMM \t tpsAssHm \t tpsResHm \t storHm");
        for (int i = 0; i < nF; i++) {
            System.out.println(nElmts[i] + " , " + tpsAssFMM[i] + " , " + tpsResFMM[i] + " , " + storFMM[i] + " , " + tpsAssHm[i] + " , " + tpsResHm[i] + " , " + storHm[i]);
        }
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();

    }

}
