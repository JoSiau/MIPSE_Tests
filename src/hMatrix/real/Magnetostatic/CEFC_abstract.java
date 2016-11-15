/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Magnetostatic;

import g2elab.mipse.analytical.singularity.SimpleTetraedreGradientPotential;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHgrad;
import g2elab.mipse.meshCore.IO.paraview.ExportVtkHgrad;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.GradNodalDeg1;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.meshCore.gaussValues.GaussPointsPositions;
import g2elab.mipse.meshCore.gaussValues.GaussPointsValues;
import g2elab.mipse.meshCore.gaussValues.GaussPointsWeights;
import g2elab.mipse.meshCore.gaussValues.VectorGaussPointsValues;
import g2elab.mipse.meshCore.quantity.RealNodalQuantity;
import g2elab.mipse.meshCore.quantity.RealVectorCellQuantity;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.Dipole;
import g2elab.mipse.mipseCore.integralIntegration.kernel.NegDotDG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.HCA.HmatrixHCAMagnetoStatPotentialCollocDeg1;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.Preconditionners.HmatrixLUDecomposition;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.CollocationIntegralFormulation;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.LuDecomposition;
import got.matrix.Matrix;
import got.matrix.RowVector;

import java.io.IOException;

/**
 * @author jsiau
 */
public class CEFC_abstract {

    /**
     * COMPUTE THE ASSEMBLY TO PLOT THE ASYMPTOTIC BEHAVIOR
     *
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException {
        // Sphere Creuse
        System.out.println("CEFC_abstract.java");
        String path = "D:";
        // SC_3846 SC_6955 SC_12938 SC_26420 SC_48154 SC_79380
        double deb, fin;

        ImportFlux mesh = new ImportFlux("D:/Meshs/SphereCreuse/SC_" + "3846.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();
        int d = 3;
        int n = ES.getNbNoeud();
        int nbGauss = 4;

        System.out.println("nbNoeuds= " + ES.getNbNoeud() + "\t nbElmts= " + ES.getNbElement());

        /**
         * *****************************************************************
         * SECOND MEMBRE
         */
        GradNodalDeg1 gradAlpha0 = new GradNodalDeg1(ES);
        NodalDeg1 alpha = new NodalDeg1(ES);

        double a = 2.5;
        // Compute phi0(X) = -a*x with X = (x, y, z);
        RowVector Pot = new RowVector(n);
        for (int i = 0; i < n; i++) {
            Pot.setElement(i, -a * ES.getLocalNodes()[i].getCoord(0));
        }

        RealNodalQuantity phi0 = new RealNodalQuantity(alpha, Pot);

        double[] secondMembre = new double[n];
        phi0.getMatrixValues().submatrix(0, 0, 1, n).get(secondMembre);

        // Effectue une copie
        double secondMembreb[] = new double[n];
        System.arraycopy(secondMembre, 0, secondMembreb, 0, n);

        /**
         * *********************************************************************
         * RESOLUTION
         */
        double Qi = 999;
        double muR = Qi + 1;

        int ordre = 2;
        double eps = 1e-3;
        GradNodalDeg1 gradAlpha = new GradNodalDeg1(ES);

        HmatrixHCAMagnetoStatPotentialCollocDeg1 f = new HmatrixHCAMagnetoStatPotentialCollocDeg1(gradAlpha, new NegDotDG(), new SimpleTetraedreGradientPotential(), nbGauss,
                eps, 50, 60, ordre, 2.0, true);
        deb = System.nanoTime();
        StorageHmatrix H = new StorageHmatrix(f);
        fin = System.nanoTime();
        System.out.println("Time to assembly the Hmatrix: " + (fin - deb) / 1e9);

        H.scale(Qi);
        H.addIdentity();

        StorageHmatrix HP = H.copy(true);
        deb = System.nanoTime();
        TruncationControl tol = new TruncationControl("rel", 1e-1);
        HP.Agglomerate(tol);
        HmatrixLUDecomposition Hlu = new HmatrixLUDecomposition(HP, tol);
        fin = System.nanoTime();
        System.out.println("Time to construct the Preconditionner: " + (fin - deb) / 1e9);

        FGMResReal fgmres = new FGMResReal(H, Hlu);

        fgmres.setInfoResolution(new double[]{1000, eps, 1, -50});
        double[] x = new double[n];
        fgmres.solve(x, secondMembre);
        RowVector resH = new RowVector(x);
        RealNodalQuantity phiH = new RealNodalQuantity(alpha, resH);

        /*
         * Export
         */
        ExportGmshHgrad exportPhiH = new ExportGmshHgrad(gradAlpha, path + "/phiH_" + n + ".msh");
        exportPhiH.addQuantity(phiH, "phiH");

        ExportVtkHgrad expVtkPhiH = new ExportVtkHgrad(new NodalDeg1(ES), path + "/phiH_" + n + ".vtk");
        expVtkPhiH.addQuantity(phiH, "phi");

        /**
         * *********************************************************************
         * STOCKAGE PLEIN
         */
        long beg = System.currentTimeMillis();

        CollocationIntegralFormulation IntegralTerm = new CollocationIntegralFormulation(new GradNodalDeg1(ES), new NegDotDG(), new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection()));
        IntegralTerm.matAssembly(new GradNodalDeg1(ES));
        Matrix Mat = ((StorageFull) IntegralTerm.getStore()).getMatrix();
        Mat.scale(Qi);

        for (int i = 0; i < alpha.getActiveDofCount(); i++) {
            Mat.setElement(i, i, Mat.getElement(i, i) + 1);
        }

        System.err.println("alpha.getActiveDofCount= " + alpha.getActiveDofCount() + "\t M.size=" + Mat.getRowCount() + " ; " + Mat.getColumnCount() + "\t n= " + n);

        ColumnVector res = new ColumnVector(alpha.getActiveDofCount());
        LuDecomposition MD = new LuDecomposition(Mat);

        MD.solve(new ColumnVector(secondMembreb), res);
        RealNodalQuantity phi = new RealNodalQuantity(alpha, res.transpose());
        ExportGmshHgrad exportPhi = new ExportGmshHgrad(gradAlpha, path + "/phi_" + n + ".msh");
        exportPhi.addQuantity(phi, "phi");

        /*
         * Erreur totale (méthode+compression)
         */
        resH.sub(res.transpose());
        System.out.println("NbNodes= " + n + "\t Erreur rel= " + resH.norm() / res.norm());

        /**
         * *********************************************************************
         * Check analytique
         */
        double mu0 = 1.25663706e-6;// Perméabilité du vide

        // Compute H = - grad phi (magnetic field)
        RealVectorCellQuantity M = phi.grad();
        RealVectorCellQuantity MH = phiH.grad();

        // Compute M = (mu_r -1). H
        M.scale(-(muR - 1));
        MH.scale(-(muR - 1));

        // Compute GPV
        VectorGaussPointsValues MPG = M.computeGaussPointValues(nbGauss);
        VectorGaussPointsValues MPGH = MH.computeGaussPointValues(nbGauss);

        // Calculate the Field
//        ComputeFieldPointsWithDipoleKernel I = new ComputeFieldPointsWithDipoleKernel(MPG);
//        ComputeFieldPointsWithDipoleKernel IH = new ComputeFieldPointsWithDipoleKernel(MPGH);
        // at the point (0, 0, 0)
        ColumnVector Hred = FieldPointAir(MPG, new double[]{0, 0, 0});
        ColumnVector HredH = FieldPointAir(MPGH, new double[]{0, 0, 0});

        // H0 = -grad phi0
        ColumnVector H0 = new ColumnVector(new double[]{a, 0, 0});

        // H = H0 + Hred
        Hred.add(H0);
        HredH.add(H0);

        System.out.println("H: \n" + Hred.transpose());

        System.out.println("HH: \n" + HredH.transpose());

        ColumnVector resi = new ColumnVector(d);
        for (int i = 0; i < d; i++) {
            resi.setElement(i, Hred.getElement(i) - HredH.getElement(i));
        }

        System.out.println("Erreur relative des H= " + resi.norm() / Hred.norm());

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    /**
     * Calcul le champ H au point dans l'air sans correction annalytique Cette
     * methode est pour les points qui situent loins de region active
     *
     * @param MPG
     * @param cibleCood
     * @return
     */
    public static ColumnVector FieldPointAir(VectorGaussPointsValues MPG, double[] cibleCood) {
        int nbGaussSource = MPG.getNbGauss();
        ColumnVector fieldPoint = new ColumnVector(3);
        Dipole dipole = new Dipole();
        // Recuperation maillage
        ElementSetHomogene mesh = MPG.getElementSet();

        // Gauss point postion computation
        GaussPointsPositions GPP = new GaussPointsPositions(mesh, nbGaussSource);

        // Gauss points weight computation
        GaussPointsWeights GPWSource = new GaussPointsWeights(mesh, nbGaussSource);


        // Values de kernel aux points de Gauss
        GaussPointsValues kerGaussPoints = dipole.getKernelGPV(mesh, MPG, new ColumnVector(cibleCood), GPP);
        //Passe au repere global
        kerGaussPoints.scale(GPWSource);

        fieldPoint = kerGaussPoints.getInfluenceM();
//        //Integration des elements sur point cible
//        Matrix integraElementsPoint = kerGaussPoints.getMatrixInfluenceM();
//        
//        ColumnVector fieldPoint1 = new ColumnVector(3);
//        for (int i = 0; i < mesh.getNbElement(); i++) {
//            fieldPoint1.add(integraElementsPoint.column(i));
//        }


        return fieldPoint;
    }

}

/*
 nNode, nElmt , Tps Resolution(en sec.)
 3846 , 11885 , 165.829
 6955 , 23913 , 391.918
 12938 , 51445 , 1138.598
 26420 , 113486 , 3108.899
 48154 , 218385 , 7336.522
 375788 , 79380 , 14538.232
 */
