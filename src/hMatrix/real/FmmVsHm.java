/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real;

import g2elab.mipse.meshCore.IO.gmsh.ExportGmshCell;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.matrixCompression.octree.repartitionElements.RepartitionElemNbMaxElem;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFMM;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFMMGalerkine;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import g2elab.mipse.numericalTools.preconditioner.PrecondIdentityReal;
import g2elab.mipse.tools.files.Ecriture;
import got.matrix.ColumnVector;
import got.matrix.LuDecomposition;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author jsiau
 */
public class FmmVsHm {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        double deb, fin, debFull = 0, finFull = 0, debFMM = 0, finFMM = 0, debHM = 0, finHM = 0;
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().setNbTaches(1);
        String type;
        //*
        type = "sphere_";
        // 6481, 12529, 34449, 62713, 132241, 301791;

        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, "D:/Meshs/Sphere_gmsh/" + type + "62713.msh");
        ElementSetHomogene ES = mesh1.createHomogeneESet();
        /*/
         type = "plaque";
         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, "D:/Meshs/plaque/"+type+"_13532.msh");
         ElementSetHomogene ES = mesh1.createHomogeneESet();
         //        type = "cube";  
         //        ImportFlux mesh1 = new ImportFlux( "D:/Meshs/Cube/"+type+"_18928.DEC");
         //        ElementSetHomogene ES = (ElementSetHomogene) mesh1.getRegions(0).getElementSet();
         //*/
        Cell C = new Cell(ES);

        int n = ES.getNbElement();
        // Check if the full matrix can be store 
        boolean computeFull = false;

        System.out.println("nbNoeuds= " + ES.getNbNoeud() + "\t nbElmts= " + ES.getNbElement());

        /*
         * SECOND MEMBRE
         */
        double[] secondMembre = new double[n];
        for (int i = 0; i < n; i++) {
            secondMembre[i] = 1;
        }
        ColumnVector resFull = new ColumnVector(n);
        ExportGmshCell exportPhi = null;
        if (computeFull) {
            System.out.println("RESOLUTION PLEINE");
            /*
             * FULL STORAGE
             */
            //*
            debFull = System.nanoTime();
            deb = System.nanoTime();
            GalerkinIntegralFormulationFull IF = new GalerkinIntegralFormulationFull(C, C, new MultG(), new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3);
            IF.assembly();
            fin = System.nanoTime();
            System.err.println("Time to compute the FullStorage = " + (fin - deb) / 1e9);

            deb = System.nanoTime();
            LuDecomposition MD = new LuDecomposition(((StorageFull) IF.getStore()).getMatrix());
            MD.solve(new ColumnVector(secondMembre), resFull);
            fin = System.nanoTime();
            System.err.println("Time to solve the FullStorage = " + (fin - deb) / 1e9);
            finFull = System.nanoTime();

            double resT[] = new double[n];
            resFull.get(resT);
            try {
                Ecriture fic = new Ecriture("D:/ChargeFull_" + type + "" + n + ".out");
                fic.ecrire(resT, ';');
                fic.close();
            } catch (IOException ex) {
                Logger.getLogger(FmmVsHm.class.getName()).log(Level.SEVERE, null, ex);
            }
            /*/
             double resT[] = new double[n];
             try {
             Lecture fic = new Lecture("D:/ChargeFull_"+type+""+n+".out");
             resT = fic.lireLigneDoubleVecteur(';');
             fic.close();
             } catch (IOException ex) {
             Logger.getLogger(FmmVsHm.class.getName()).log(Level.SEVERE, null, ex);
             }
        
             ColumnVector res = new ColumnVector(resT);
             //*/

            exportPhi = new ExportGmshCell(C, "D:/V.msh");
            exportPhi.addQuantity(new RealScalarCellQuantity(C, resFull.transpose()), "charge");
        }
        System.out.println("RESOLUTION FMM");
        /*
         * FMM
         */
        debFMM = System.nanoTime();
        deb = System.nanoTime();
        GalerkinIntegralFormulationFMM f = new GalerkinIntegralFormulationFMM(C, C, new MultG(),
                new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3, 3,
                new RepartitionElemNbMaxElem(30), 1.0, 2, 1, 0);
        // comparer adaptatif et niveau constant
        // verifier la complexite du nivo constant
        // comparer les options memoire 

        f.assembly();
        fin = System.nanoTime();
        System.err.println("Time to compute the FMM = " + (fin - deb) / 1e9);
        StorageFMMGalerkine Mf = (StorageFMMGalerkine) f.getStore();
        double storFMM = Mf.getMemoryUsed();
        System.out.println("Storage FMM = " + storFMM);

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
        double configRes[] = {100, 1e-8, 1, -50};
        solvFMM.setInfoResolution(configRes);
        double xFMM[] = new double[n];
        xFMM = solvFMM.solve(xFMM, secondMembre);
        fin = System.nanoTime();
        System.err.println("Time to solve the FMM = " + (fin - deb) / 1e9);
        ColumnVector resFMM = new ColumnVector(xFMM);
        finFMM = System.nanoTime();
        System.out.println("Total time FMM = " + (finFMM - debFMM) * 1e-9);

        double eps;
        if (computeFull) {
            ColumnVector tmp = resFMM.copy();
            tmp.sub(resFull);
            eps = tmp.norm() / resFull.norm();
            System.out.println("Erreur relative FMM = " + eps);
        } else {
            eps = 1e-2;// = 1%
        }
        exportPhi = new ExportGmshCell(C, "D:/V_FMM.msh");
        exportPhi.addQuantity(new RealScalarCellQuantity(C, resFMM.transpose()), "charge FMM");

        System.out.println("RESOLUTION HMATRIX");
        /*
         * HCA 
         */
        debHM = System.nanoTime();
        deb = System.nanoTime();
        GalerkinIntegralFormulationHCA f1 = new GalerkinIntegralFormulationHCA(C, C, new MultG(),
                new SelfElementFixedGauss(3, new AnalyticalCorrection()), 3, 3,
                eps, 50, 30, Math.abs((int) Math.floor(Math.log10(eps))), 2.0, false);
        f1.assembly();
        StorageHmatrix Mh = (StorageHmatrix) f1.getStore();
        Mh.Agglomerate(new TruncationControl("rel", eps));
        fin = System.nanoTime();
        System.err.println("Time to compute the Hmatrix = " + (fin - deb) / 1e9);
        double storHm = Mh.getMemoryUsed();
        /*
         * Preconditionneur
         */
        /*
         deb = System.nanoTime();
         Hmatrix precondH = Mh.copy(false);
         TruncationControl tolPrecond = new TruncationControl("rel",1e-1);
         precondH.Coarsen(tolPrecond);
         HmatrixLUDecomposition hLu = new HmatrixLUDecomposition(precondH,tolPrecond);
         fin = System.nanoTime();
         System.err.println("Time to compute the Hmatrix Preconditionner = " + (fin - deb) / 1e9);
         /*/
        PrecondIdentityReal hLu = new PrecondIdentityReal();
        //*/
        /*
         * RESOLUTION
         */
        deb = System.nanoTime();
        FGMResReal solvMH = new FGMResReal(Mh, hLu);
        solvMH.setInfoResolution(configRes);
        double xHM[] = new double[n];
        xHM = solvMH.solve(xHM, secondMembre);
        fin = System.nanoTime();

        System.err.println("Time to solve the Hmatrix = " + (fin - deb) / 1e9);
        ColumnVector resH = new ColumnVector(xHM);
        finHM = System.nanoTime();

        exportPhi = new ExportGmshCell(C, "D:/V_HM.msh");
        exportPhi.addQuantity(new RealScalarCellQuantity(C, resH.transpose()), "charge HM");

        if (computeFull) {
            resH.sub(resFull);
            System.out.println("Erreur relative Hmatrix = " + resH.norm() / resFull.norm());
        }

        System.out.println("Total  FULL\tFMM\tHmatrix");
        System.out.format("Time   %.3f\t%.3f\t%.3f", (finFull - debFull) * 1e-9, (finFMM - debFMM) * 1e-9, (finHM - debHM) * 1e-9);
        System.out.println("\n Storage FMM= " + storFMM + "\t Storage HMatrix= " + storHm);

        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

}
/* SANS RECOMP.
 ------------------------------6198
 Time   493,498	12,487	12,009
 Storage FMM= 5.034802E7	 Storage HMatrix= 1.532076E7

 ------------------------------12138
 Total  FULL	FMM	Hmatrix
 Time   3604,740	18,781	23,421
 Storage FMM= 8.41128E7	 Storage HMatrix= 3.4779152E7
 ------------------------------33806
 Time   0,000	77,015	78,236
 Storage FMM= 2.6868176E8	 Storage HMatrix= 1.06720568E8
 */
