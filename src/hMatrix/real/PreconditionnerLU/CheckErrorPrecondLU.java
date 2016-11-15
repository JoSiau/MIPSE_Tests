/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.PreconditionnerLU;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.Preconditionners.HmatrixLUDecomposition;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationACA;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import g2elab.mipse.numericalTools.preconditioner.PrecondIdentityReal;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.RowVector;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * CHECK LES ERREURS DU PRECONDITIONNEUR
 *
 * @author siau
 */
public class CheckErrorPrecondLU {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {

        /////////////////////////////////////////////////////////////////////////
        ////////////////////////// HMATRIX CONSTRUCTION /////////////////////////
        /////////////////////////////////////////////////////////////////////////
        String meshDir = new java.io.File(".").getCanonicalPath();
        /*
        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
        ImportFlux mesh = new ImportFlux(meshDir+"/sphere/SPHERE_572.DEC");       
        ElementSet ES  =  mesh.getRegions(0).getElementSet();      
        /*/
//        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere_gmsh";
//        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/SPHERE_24091.msh");  
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
//        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/plaque/plaque_13532.msh"); 
        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(1, meshDir + "/plaque/plaque_24228.msh");
//        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/SPHERE_132241.msh");                    
        ElementSetHomogene ES = mesh1.createHomogeneESet();
        //*/

        /////////////////////////////////////////////////////////////////////////
        ///////////////////////// INITIALISATION DONNEES ////////////////////////
        /////////////////////////////////////////////////////////////////////////
        double deb, fin;
        int d = 3;
        int n = ES.getNbElement();
        int nbGauss = 7;
        double epsAssem = 1e-5;
        double epsFGMRes = 1e-10;
        double epsPrecond[] = new double[14];
        epsPrecond[0] = 3e-3;
        epsPrecond[1] = 1e-2;
        epsPrecond[2] = 3e-2;
        epsPrecond[3] = 8e-2;
        epsPrecond[4] = 8.5e-2;
        epsPrecond[5] = 9e-2;
        epsPrecond[6] = 9.5e-2;
        epsPrecond[7] = 1e-1;
        epsPrecond[8] = 1.2e-1;
        epsPrecond[9] = 1.4e-1;
        epsPrecond[10] = 1.6e-1;
        epsPrecond[11] = 1.8e-1;
        epsPrecond[12] = 2e-1;
        epsPrecond[13] = 3e-1;

        Cell Ces = new Cell((ElementSetHomogene) ES);

        GestionnaireTaches.getGestionnaireTaches().setNbTaches(7);
        
       /*
        * SECOND MEMBRE
        */
        RowVector valueEl = new RowVector(ES.getNbElement());
        /*
        valueEl.setAllElements(1);
        /*/
        for (int i = 0; i < valueEl.getColumnCount(); i++)
            valueEl.setElement(i, Math.random());
        //*/ 
        RealScalarCellQuantity V0 = new RealScalarCellQuantity(Ces, valueEl);
        // Ici on calcule (int(langda*V0))
        RowVector Bb = (RowVector) V0.projWithExplicitDof(Ces, nbGauss);
        double[] b = new double[Bb.getColumnCount()];
        Bb.get(b);
        /////////////////////////////////////////////////////////////////////////
        ////////////////////////// STOCKAGE PLEIN /////////////////////////
        /////////////////////////////////////////////////////////////////////////
        System.out.println("Calcul stockage Plein");
        GalerkinIntegralFormulation IFCell = new GalerkinIntegralFormulationFull(Ces, Ces, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss);
        IFCell.assembly();


        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ FGMRes AVEC PRECOND");
        FGMResReal t1 = new FGMResReal(IFCell.getStore().getMatrixPrecond(null), new PrecondIdentityReal());
        double configRes[] = {100, epsFGMRes, 1, -50};
        t1.setInfoResolution(configRes);
        deb = System.nanoTime();
        double r1[] = t1.solve(new double[Ces.getActiveDofCount()], b);
        fin = System.nanoTime();
        System.out.println("----- Time to compute FGMRES avec Precond: " + (fin - deb) * 1e-9);
        ColumnVector Ex = new ColumnVector(r1);

        IFCell = null;


        /////////////////////////////////////////////////////////////////////////
        ////////////////////////// HMATRIX CONSTRUCTION /////////////////////////
        /////////////////////////////////////////////////////////////////////////
        //*                
        int kmax = 50, nmin = 30;
        GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(Ces, Ces, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss,
                epsAssem, kmax, nmin);
        f.assembly();
        // Hmatrix(Matrix P, HmatrixCompatible HC, double epsilon, int k, int dof, int leafSize)
        StorageHmatrix H0 = (StorageHmatrix) f.getStore();

        deb = System.nanoTime();
        H0.Agglomerate(new TruncationControl("rel", epsAssem));
        fin = System.nanoTime();
        System.out.println("Time to Agglomerate H0= " + (fin - deb) * 1e-9);
        System.out.println("Csp= " + H0.getCsp());

//        Matrix A = f.getSubMatrix(0, 0, n, n);


        ///////////////////////////////////////////////////////////////////////// 
        ////////////////////////////   AGGLOMERATION   //////////////////////////
        ///////////////////////////////////////////////////////////////////////// 

        int epsLength = epsPrecond.length;
        ColumnVector r[] = new ColumnVector[epsLength];
        for (int i = 0; i < r.length; i++)
            r[i] = new ColumnVector(n);


        int csp[] = new int[epsLength];

        double TpsSolve[] = new double[epsLength];


        for (int k = 0; k < epsLength; k++) {
            TruncationControl tol = new TruncationControl("rel", epsPrecond[k]);
            System.err.println("Go to Coarsening: eps= " + epsPrecond[k]);

            StorageHmatrix Htmp = H0.copy(true);

            deb = System.nanoTime();
            Htmp.Agglomerate(tol);
            fin = System.nanoTime();
            System.out.println("Coarsening's over: " + (fin - deb) * 1e-9);
//            Htmp.getError("rel");

            csp[k] = Htmp.getCsp();

            System.out.println("Go to LU");
            deb = System.nanoTime();
            HmatrixLUDecomposition hLu = new HmatrixLUDecomposition(Htmp, tol);
            fin = System.nanoTime();
            System.out.println("LU's over: " + (fin - deb) * 1e-9);


            System.out.println("Résolution Directe");
            deb = System.nanoTime();
            hLu.solve(r[k], Bb.transpose());
            fin = System.nanoTime();
            TpsSolve[k] = (fin - deb) * 1e-9;
            System.out.println("Time to solve= " + TpsSolve[k]);

            System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++");
//            fprint("H"+k+".out",Htmp.getBCT().getTransfertStruct(false));
        }
        PrintWriter fic = new PrintWriter(new BufferedWriter(new FileWriter("Tps" + n + ".out")));

        System.out.println("n= " + n);
        double normEx = Ex.norm();
        for (int i = 0; i < epsLength; i++) {
            r[i].sub(Ex);
            double normi = r[i].norm();
            System.out.format("Erreur de " + epsPrecond[i] + "\t absolue= %1.4e \t relative= %1.4e \t TpsSolve= %1.4e \t Csp= %d %n", normi, normi / normEx, TpsSolve[i], csp[i]);
            fic.println(i + " " + epsPrecond[i] + "  " + normi / normEx + " " + csp[i]);
        }
        fic.close();

        GestionnaireTaches.getGestionnaireTaches().stop();
    }


}
