/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.PreconditionnerLU;

import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.Element;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.Preconditionners.HmatrixLUDecomposition;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TransfertStruct;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import g2elab.mipse.numericalTools.preconditioner.PrecondIdentityReal;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;
import got.matrix.RowVector;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * CHECK LE COMPORTEMENT DU PRECONDITIONNEUR (nbr Iterations)
 *
 * @author siau
 */
public class CheckIterationPrecondLU {

    /**
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        System.out.println("CheckIterationPrecondLU");
        ///////////////////////////////////////////////////////////////////////// 
        ////////////////////////// HMATRIX CONSTRUCTION /////////////////////////
        ///////////////////////////////////////////////////////////////////////// 
        String meshDir = new java.io.File(".").getCanonicalPath();
        /*
         meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
         //        ImportFlux mesh = new ImportFlux(meshDir+"/sphere/SPHERE_572.DEC");   
         ImportFlux mesh = new ImportFlux(meshDir+"/sphere_flux/25838.DEC");   
         ElementSet ES  =  mesh.getRegions(0).getElementSet();      
         /*/
//        meshDir += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere_gmsh"; 
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
//        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(1,meshDir+"/sphere/sphere_3197.msh"); 
        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(1, meshDir + "/plaque/plaque_24228.msh");
//        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(1, meshDir + "/plaque/plaque_4548.msh");
//        ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(1,meshDir+"/plaque/plaque_95006.msh");                    
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
        double epsFGMRes = -1e-8;
        double epsPrecond[] = new double[14];
        epsPrecond[0] = 1e-3;
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

        ///////////////////////////////////////////////////////////////////////// 
        ////////////////////////// HMATRIX CONSTRUCTION /////////////////////////
        ///////////////////////////////////////////////////////////////////////// 
        //*        
        Matrix X = new Matrix(n, d);
        Element E;
        for (int i = 0; i < X.getRowCount(); i++) {
            // Load the next element
            E = ES.getElements(i);
            // Get the gravity center and Copy it in X !
            X.setRow(i, E.getCentroid());
        }
        // f1 f = new f1(X);  
        int kmax = 30, nmin = 30;
        Cell C = new Cell(ES);
        GalerkinIntegralFormulation f = new GalerkinIntegralFormulationHCA(C, C, new MultG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()), nbGauss, nbGauss,
                epsAssem, kmax, nmin, 4);
        f.assembly();
        // Hmatrix(Matrix P, HmatrixCompatible HC, double epsilon, int k, int dof, int leafSize)
        StorageHmatrix H0 = (StorageHmatrix) f.getStore();
        //fprint("Ha.out",H0.getTransfertStruct(false));
        deb = System.nanoTime();
        H0.Agglomerate(new TruncationControl("rel", epsAssem));
        fin = System.nanoTime();
        System.out.println("Time to Agglomerate H0= " + (fin - deb) * 1e-9);
        System.out.println("Csp= " + H0.getCsp());

        //fprint("H.out",H0.getTransfertStruct(false));
//        Matrix A = f.getSubMatrix(0, 0, n, n);
        Cell Ces = new Cell((ElementSetHomogene) ES);
//        GalerkinIntegralFormulation IFCell = new GalerkinIntegralFormulation(Ces, new MultG(),new SimpleTriangleConstantChargePotential());
//        IFCell.matAssembly(Ces, nbGauss, nbGauss);

        /*/ 
         Cell Ces = new Cell((ElementSetHomogene)ES);
         GalerkinIntegralFormulation IFCell = new GalerkinIntegralFormulation(Ces, new MultG(),new SimpleTriangleConstantChargePotential());
         IFCell.matAssembly(Ces, nbGauss, nbGauss);
         Matrix A = ((StorageFull) IFCell.getStore()).getMatrix();
         Matrix Ai = A;
         int nlvl=2;
         double eps=1e-20;
         TruncationControl tol = new TruncationControl(eps);   
        
         //        BlockClusterTree Bc0 = new BlockClusterTree(A.copy(),new TruncationControl(1),nlvl);
         //        Hmatrix H0 = new Hmatrix(Bc0 ,n);
        
         BlockClusterTree Bc0 = new BlockClusterTree(n,nlvl,d);
         HmatrixGalerkinIntegral f = new HmatrixGalerkinIntegral(ES.getElements(),new MultG(),new SimpleTriangleConstantChargePotential(),nbGauss,nbGauss);
         Hmatrix H0 = new Hmatrix(Bc0, f, n, eps, 1);
         //*/
        double tpsLU[] = new double[epsPrecond.length];
        double tpsAgglo[] = new double[epsPrecond.length];
        double tpsFGMResPr[] = new double[epsPrecond.length];

        double mem[] = new double[epsPrecond.length];
        int nit[] = new int[epsPrecond.length];
        int csp[] = new int[epsPrecond.length];


        /*
         * SECOND MEMBRE
         */
        RowVector valueEl = new RowVector(ES.getNbElement());
        /*
         valueEl.setAllElements(1);
         /*/
        for (int i = 0; i < valueEl.getColumnCount(); i++) {
            valueEl.setElement(i, Math.random());
        }
        //*/ 
        RealScalarCellQuantity V0 = new RealScalarCellQuantity(Ces, valueEl);
        // Ici on calcule (int(langda*V0))
        RowVector Bb = (RowVector) V0.projWithExplicitDof(Ces, nbGauss);
        double[] b = new double[Bb.getColumnCount()];
        Bb.get(b);

        GestionnaireTaches.getGestionnaireTaches().setNbTaches(4);
        ///////////////////////////////////////////////////////////////////////// 
        ////////////////////////////   AGGLOMERATION   //////////////////////////
        ///////////////////////////////////////////////////////////////////////// 

        double[][] r = new double[epsPrecond.length][];
        for (int k = 0; k < epsPrecond.length; k++) {
            TruncationControl tol = new TruncationControl("rel", epsPrecond[k]);
            System.err.println("Go to Coarsening: eps= " + epsPrecond[k]);

            StorageHmatrix Htmp = H0.copy(true);

            deb = System.nanoTime();
            Htmp.Agglomerate(tol);
            fin = System.nanoTime();
            tpsAgglo[k] = (fin - deb) * 1e-9;
            System.out.println("Coarsening's over: " + tpsAgglo[k]);
            csp[k] = Htmp.getCsp();
            System.out.println("Csp= " + csp[k]);
//        Htmp.getError("rel");

            //((double)this.mem*8/1024)/(double)this.ndof
            mem[k] = ((double) Htmp.getMemoryUsed() * 8 / 1024) / n;

            //fprint("H"+k+".out",Htmp.getBCT().getTransfertStruct(false));
            ///////////////////////////////////////////////////////////////////////// 
            /////////////////////// CHECK THE H-MATRICES ERROR //////////////////////
            ///////////////////////////////////////////////////////////////////////// 
//        System.out.println("Checking the H-matrices error");
//        
//        TransfertStruct Tr[] = H0.getBCT().getTransfertStruct(true);
//        System.out.println("Tr.length= "+Tr.length);
//        fprint("transfert0.out",Tr);
//        
//        Matrix exact = new Matrix(n,n);
//        for(int i=0;i<Tr.length;i++)
//            exact.setSubmatrix(Tr[i].getRowStart(), Tr[i].getColStart(), Tr[i].getM());
//        
//        Matrix ex = exact.copy();
//        ex.sub(A);
//        System.out.println("Error de la hmatrix initiale exacte= "+ex.norm()+" \t relative= "+ex.norm()/A.norm());     
//        Tr = Htmp.getBCT().getTransfertStruct(true);
//        System.out.println("Tr.length= "+Tr.length);
//        fprint("transfert0.out",Tr);
//        
//        exact = new Matrix(n,n);
//        for(int i=0;i<Tr.length;i++)
//            exact.setSubmatrix(Tr[i].getRowStart(), Tr[i].getColStart(), Tr[i].getM());
//        
//        ex = exact.copy();
//        ex.sub(A);
//        System.out.println("Error de la hmatrix agglomerée exacte= "+ex.norm()+" \t relative= "+ex.norm()/A.norm());     
            ///////////////////////////////////////////////////////////////////////// 
            //////////////////////////// LU DECOMPOSITION ///////////////////////////
            ///////////////////////////////////////////////////////////////////////// 
            System.out.println("Go to LU");
            deb = System.nanoTime();
            HmatrixLUDecomposition hLu = new HmatrixLUDecomposition(Htmp, tol);
            fin = System.nanoTime();
            tpsLU[k] = (fin - deb) * 1e-9;
            System.out.println("LU's over: " + tpsLU[k]);
            System.out.println("Memory used for templates = " + (hLu.getMemoryUsed() / (double) 1024 / 1024) + " Mo");

            /////////////////////////////////////////////////////////////////////////
            System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ FGMRes AVEC PRECOND");

            FGMResReal t1 = new FGMResReal(H0, hLu);

            double configRes[] = {100, epsFGMRes, 1, -50};
            t1.setInfoResolution(configRes);

            deb = System.nanoTime();
            r[k] = t1.solve(new double[Ces.getActiveDofCount()], b);
            fin = System.nanoTime();
            tpsFGMResPr[k] = (fin - deb) * 1e-9;
            System.out.println("----- Time to compute FGMRES avec Precond: " + tpsFGMResPr[k]);
            nit[k] = t1.getNumberOfSubIterations() + 1;

            System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++");
            System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        }

        PrecondIdentityReal Id = new PrecondIdentityReal();
        FGMResReal t2 = new FGMResReal(H0, Id);

        double configRes[] = {100, epsFGMRes, 1, -50};
        t2.setInfoResolution(configRes);

        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++ FGMRes SANS PRECOND");
        deb = System.nanoTime();
        double r_ref[] = t2.solve(new double[Ces.getActiveDofCount()], b);
        fin = System.nanoTime();
        double tpsFGMResSPr = (fin - deb) * 1e-9;
        System.out.println("----- Time to compute FGMRES sans Precond: " + tpsFGMResSPr);

        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();

        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println(" Nbr Elts = \t|      Agglomération   \t \t| Décomposition\t|     \t  FGMRes      \t\t| \t Temps  \t| Constante de \t|");
        System.out.format(" %05d \t\t|    Temps  \t|    Stockage  \t| \t LU \t| Temps \t| Nbr Iter° \t| \t Total  \t|  Sparsité \t| %n", n);
        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        for (int i = 0; i < epsPrecond.length; i++) {
            System.out.format(" eps=" + epsPrecond[i] + " \t| %.9f \t| \t %.2f \t| %.9f \t| %.9f \t| \t %d \t| \t %.9f  \t| \t %d \t| %n", tpsAgglo[i], mem[i], tpsLU[i], tpsFGMResPr[i], nit[i], (tpsAgglo[i] + tpsLU[i] + tpsFGMResPr[i]), csp[i]);
            System.out.println("-----------------------------------------------------------------------------------------------------------------------------------------");
        }
        System.out.println("-----------------------------------------------------------------------------------------------------------------------------------------");
        System.out.format(" Sans Precond. \t| \t 0 \t| \t 0 \t| \t 0 \t| %.9f \t| \t %d \t| \t %.9f  \t| \t %d \t| %n", tpsFGMResSPr, (t2.getNumberOfSubIterations() + 1), tpsFGMResSPr, H0.getCsp());
        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("LU H-matrix Complète");
        deb = System.nanoTime();
        HmatrixLUDecomposition hLu = new HmatrixLUDecomposition(H0, new TruncationControl("rel", epsAssem));
        fin = System.nanoTime();
        System.out.println("LU's assembly : " + (fin - deb) * 1e-9);
        deb = System.nanoTime();
        ColumnVector x = new ColumnVector(n);
        hLu.solve(x, Bb.transpose());
        fin = System.nanoTime();
        System.out.println("LU's solver direct : " + (fin - deb) * 1e-9);

        System.out.println("");
        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("Check the errors !");
        ColumnVector ref = new ColumnVector(r_ref);
        double nr = ref.norm();

        for (int i = 0; i < epsPrecond.length; i++) {
            ColumnVector rc = new ColumnVector(r[i]);

            rc.sub(ref);
            System.out.println("Erreur a la precision preconditionneur " + epsPrecond[i] + " donne: " + (rc.norm() / nr));
        }
        x.sub(ref);
        System.out.println("Erreur de la complete est de : " + (x.norm() / ref.norm()));

        System.err.println("Bien entendu les résultats doivent être équivalents !");

    }

    /**
     * @param name
     * @param T
     */
    public static void fprint(String name, TransfertStruct T[]) {
        PrintWriter fic;
        try {
            fic = new PrintWriter(new BufferedWriter(new FileWriter(name)));
            for (int i = 0; i < T.length; i++) {
                fic.println(T[i]);
            }
            fic.close();
        } catch (IOException ex) {
        }
    }

    /**
     * @param name
     * @param A
     */
    public static void fprint(String name, Matrix A) {
        PrintWriter fic;
        try {
            fic = new PrintWriter(new BufferedWriter(new FileWriter(name)));
            fic.println(A);
            fic.close();
        } catch (IOException ex) {
        }
    }
}
