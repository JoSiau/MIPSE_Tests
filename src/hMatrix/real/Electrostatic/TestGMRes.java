/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Electrostatic;

// La méthode de résolution

import g2elab.mipse.meshCore.IO.gmsh.ExportGmshCell;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMesh;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import g2elab.mipse.numericalTools.preconditioner.PrecondIdentityReal;
import got.matrix.Matrix;
import got.matrix.RowVector;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * @author siau
 */
public class TestGMRes {

    /**
     * @param args
     */
    public static void main(String[] args) {
        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().setNbTaches(2); // Mono proc
        double tPlein = 0, tHmat = 0;
        double deb, fin;
        int nbGauss = 7;

        String file = null;
        try {
            file = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(TestGMRes.class.getName()).log(Level.SEVERE, null, ex);
        }
        file += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
        /*
         // ImportFlux F = new ImportFlux(file+"/sphere/SPHERE_3538.DEC");
         //        ImportFlux F = new ImportFlux(file+"/PLAQUE2000.DEC");
         ImportFlux F = new ImportFlux(file+"/sphere/SPHERE_3538.DEC");
         ElementSurfSetHomogene mesh = (ElementSurfSetHomogene) F.getRegion(0).getElementSet();
         /*/
        System.out.println("Loading GMSH mesh");
        ImportGmshMesh mesh1 = new ImportGmshMesh(file + "/plaque/plaque_8694.msh");
        mesh1.meshSummary();
        ElementSetHomogene mesh = (ElementSetHomogene) mesh1.getRegion(0).getElementSet();
        System.out.println("GMSH mesh loaded");
        //*/

        // Assemblage
        Cell langda0 = new Cell(mesh);
        System.out.println("Calculating IF");
        deb = System.nanoTime();
        GalerkinIntegralFormulation IF = new GalerkinIntegralFormulationFull(langda0, langda0, new MultG(),
                new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss);
        IF.assembly();
        fin = System.nanoTime();
        System.out.println("IF done");
        tPlein = (fin - deb) / 1e9;
        System.out.println("Size store= " + langda0.getActiveDofCount() + " , " + langda0.getActiveDofCount());
        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

        // Second membre
        RowVector valueEl = new RowVector(mesh.getNbElement());
        valueEl.setAllElements(1);
        RealScalarCellQuantity V0 = new RealScalarCellQuantity(langda0, valueEl);
        // Ici on calcule (int(langda*V0))
        RowVector B = (RowVector) V0.projWithExplicitDof(langda0, nbGauss);
        double[] b = new double[B.getColumnCount()];
        B.get(b);

        /*
         * THE PRECOND IS THE IDENTITY
         */
        PrecondIdentityReal Id = new PrecondIdentityReal();
        FGMResReal solverPlein = new FGMResReal(IF.getStore(), Id);

        /* Parametres du solveurFGMRes
         1) Nombre d'iterations maximal du solveur
         2) Critère d'arrêt (Norme du résidu)
         3) Affichage de l'évolution du solveur
         0 -> Pas d'affichage
         1 -> Affichage
         4) Taille du sous-espace de Krylov (plus petit que le nombre d'inconnues du système)
         */
        //******************* Plutot 1e-6 (au minimum) ******************
        double configRes[] = {100, 1e-8, 1, -50};
        solverPlein.setInfoResolution(configRes);

        deb = System.nanoTime();
        double resPlein[] = solverPlein.solve(new double[langda0.getActiveDofCount()], b);
        fin = System.nanoTime();
        tPlein += (fin - deb) / 1e9;
        System.out.println("----- Time to compute FGMRES Plein: " + (fin - deb) / 1e9);
        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

        // Construction H-matrix
        int d = 3;//ES.getElements(0).getDim();
        int n = mesh.getNbElement();
        // System.out.println(" n = "+n+"\t d = "+d);

        Matrix X = new Matrix(n, d);
        for (int i = 0; i < n; i++) {
            X.setRow(i, mesh.getElements(i).getCentroid());
        }

        /*
         * ASSEMBLAGE DE LA HMATRIX
         */
        int kmax = 100;
        int nmin = kmax;
        double eps = 1e-5;
        deb = System.nanoTime();
        GalerkinIntegralFormulationHCA HGI = new GalerkinIntegralFormulationHCA(langda0, langda0, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss,
                nbGauss, eps, kmax, nmin, 2);
        HGI.assembly();
        StorageHmatrix Hm = (StorageHmatrix) HGI.getStore();
        fin = System.nanoTime();
        System.out.println("----- Time to Assembly the Hmatrix: " + (fin - deb) / 1e9);
        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        tHmat = (fin - deb) / 1e9;
        /*
         * EFFECTUE UNE AGGLOMERATION A PRECISION EGALE
         */
        deb = System.nanoTime();
        Hm.Agglomerate(new TruncationControl(eps));
        fin = System.nanoTime();
        System.out.println("----- Time to Agglomerate the Hmatrix: " + (fin - deb) / 1e9);
        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        tHmat += (fin - deb) / 1e9;
        /*
         * GMRES
         */
        FGMResReal solverHmatrix = new FGMResReal(Hm, Id);
        solverHmatrix.setInfoResolution(configRes);
        /*
         * RESOLUTION
         */
        deb = System.nanoTime();
        double resHmatrix[] = solverHmatrix.solve(new double[langda0.getActiveDofCount()], b);
        fin = System.nanoTime();
        tHmat += (fin - deb) / 1e9;
        System.out.println("----- Time to compute FGMRES Hmatrix: " + (fin - deb) / 1e9);
        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");

        RowVector resgP = new RowVector(resPlein);
        RealScalarCellQuantity qgP = new RealScalarCellQuantity(langda0, resgP);
//        System.out.println(qgP.getMatrixValues());
        System.out.println("Charge totale : " + qgP.integ(nbGauss));

        RowVector resgHM = new RowVector(resHmatrix);
        RealScalarCellQuantity qgHM = new RealScalarCellQuantity(langda0, resgHM);
//        System.out.println(qgHM.getMatrixValues());
        System.out.println("Charge totale : " + qgHM.integ(nbGauss));

        // EXPORT DES RESULTAT !
        ExportGmshCell exportVP = new ExportGmshCell(langda0, "Y:/chargeP.msh");
        exportVP.addQuantity(qgP, "charge");
        ExportGmshCell exportVHM = new ExportGmshCell(langda0, "Y:/chargeHM.msh");
        exportVHM.addQuantity(qgHM, "charge");
        //*/

        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("n= " + n);
        System.out.println("Time Plein= " + tPlein);
        System.out.println("Time Hmat= " + tHmat);

        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
