/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.condVol_dielecVol;

import g2elab.mipse.circuit.solverCircuitComplex.BlocCircuit;
import g2elab.mipse.circuit.solverCircuitComplex.SolveurMaillesIndependantes;
import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmsh;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.contraints.ExternalDriveFaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceConstraint;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.geomElements.GeomElementSurf;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.meshCore.region.RegionsSet;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageSparse;
import g2elab.mipse.numericalTools.iterativeSolver.ProductComplex;
import g2elab.mipse.numericalTools.matrix.MatriceIncidence;
import g2elab.mipse.numericalTools.matrix.complex.MatrixComplexPartReIm;
import g2elab.mipse.numericalTools.matrix.real.dense.gotMatrix.GOTMatrix;
import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
import g2elab.mipse.numericalTools.vector.sparse.SparseVectorComplex;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * PEEC generalized with capacitive's effects and consideration of dielectrics.
 * Najoute pas de noeud sur la frontiere.
 *
 * @author jsiau
 */
public class PEEC_RLMPC_WIRE implements BlocCircuit {

    /**
     * Stockage de l'espace fonctionnel
     */
    private final FaceDeg1 Omega;
    /**
     * Matrice des resistances
     */
    private SparseMatrixRowReal R[] = new SparseMatrixRowReal[2];
    /**
     * Matrice des capacités
     */
    private SparseMatrixRowReal C[] = new SparseMatrixRowReal[2];
    /**
     * Matrice d'induction
     */
    private Matrix L[] = new Matrix[2];
    /**
     * Matrice de capacité
     */
    private Matrix P[] = new Matrix[2];
    /**
     * Nombre de points de gauss source et cible pour l'integration de L
     */
    private int nbPtGCiblesInduc, nbPtGSourcesInduc;
    /**
     * Nombre de points de gauss source et cible pour l'integration de P
     */
    private int nbPtGCiblesCapa, nbPtGSourcesCapa;
    /**
     * Objet portant le solveur circuit
     */
    private final SolveurMaillesIndependantes solveurCircuit;
    /**
     * Source U
     */
    private double SourceU[] = null;
    /**
     * Source I
     */
    private double SourceI[] = null;
    /**
     * Matrice de passage
     */
    private MatriceIncidence M;
    /**
     * mu_0
     */
    private final double mu0 = 4 * Math.PI * 1e-7;
    /**
     * eps_0
     */
    private static final double eps0 = 8.85418782 * 1e-12;
    /**
     * sigmaStar (conductivite equivalente) = sigma + j * omega * epsilon
     */
    private final double sigmaStar[][];
    /**
     * The external border (conductor or/and dielectric) / air
     */
    private final Cell GammaExt;
    /**
     * The interface conductor / dielectric
     */
    private final Cell GammaInt;
    /**
     * Nb de face conductrice
     */
    private final int nDofCond;
    /**
     * Nb de face dielectrique
     */
    private final int nDofDielec;
    /**
     * Table des surfaces
     */
    private double[] s;
    /**
     *
     */
    CombLineaireComplexe cl = null;
    /**
     *
     */
    MatrixComplexPartReIm Pc = null;

    /**
     * Default constructor. Use all the functions "set...()" you must do before
     * solving anything ! May the Force be with you !
     *
     * @param omega FaceDeg1
     * @param conductivity sizes = [ 2 ][nbRegions], conductivity[0][?] = sigma
     * , conductivity[1][?] = epsilon
     * @param ExternalBorder Cell
     * @param InternalBorder Cell
     * @param nDof = {nbDofConductor , nbDielectric , ... }
     */
    public PEEC_RLMPC_WIRE(FaceDeg1 omega, double[][] conductivity, Cell ExternalBorder, Cell InternalBorder, int nDof[]) {
        Omega = omega;
        this.sigmaStar = conductivity;
        GammaExt = ExternalBorder;
        GammaInt = InternalBorder;

        this.nDofCond = nDof[0];
        this.nDofDielec = nDof[1];

        solveurCircuit = new SolveurMaillesIndependantes(2);
        solveurCircuit.setBloc(this, 0, 0);
    }

    /**
     * Resolution directe du probleme
     *
     * @param f Liste des frequence a analyser
     * @return Vecteur contenant les resultats (courants complexe pour chaque
     * frequence puis tension)
     */
    public double[][] resolutionIterative(double f) {
        this.integration(f);
        this.solveurCircuit.makeCircuit();
        this.solveurCircuit.setConfiguration(new double[]{1, 100, 1e-6, 1, -50}, new double[]{0, 50});
        return this.solveurCircuit.resolutionIterative(new double[]{f}, true);
    }

    public double[][] resolutionDirecte(double f) {
        this.integration(f);
        this.solveurCircuit.makeCircuit();
        return this.solveurCircuit.resolutionDirecte(new double[]{f}, true);
    }

    /**
     * Affectation du bloc circuit electrique pure (description des elements
     * autours de la partie PEEC : R, L , C , SU, SI
     *
     * @param circuitPur Description des elements de circuit exterieur a la
     * geometrie
     */
    public void setCircuitElectrique(BlocCircuit circuitPur) {
        this.solveurCircuit.setBloc(circuitPur, 1, 1);
    }

    /**
     * Integration de la matrice R
     *
     * @param f la frequence
     */
    protected void integrationR(double f) {
        // Compute 1 / sigmaStar = conjuguee(sigmaStar) / norm(sigmaStar)^2
        double invSigma[][] = new double[2][nDofCond + nDofDielec];
        // Conductor
        double norm = Math.hypot(sigmaStar[0][0], 2 * Math.PI * f * sigmaStar[1][0]);
        double Re = sigmaStar[0][0] / (norm * norm), Im = -2 * Math.PI * f * sigmaStar[1][0] / (norm * norm);
        for (int i = 0; i < nDofCond; i++) {
            // La partie reelle 
            invSigma[0][i] = Re;
            //La partie complexe
            invSigma[1][i] = Im;
        }
        // Dielectric
        norm = Math.hypot(sigmaStar[0][1], 2 * Math.PI * f * sigmaStar[1][1]);
        Re = sigmaStar[0][1] / (norm * norm);
        Im = - 2 * Math.PI * f * sigmaStar[1][1] / (norm * norm);
        for (int i = nDofCond; i < nDofCond + nDofDielec; i++) {
            // La partie reelle 
            invSigma[0][i] = Re;
            //La partie complexe
            invSigma[1][i] = Im;
        }

        //
        // Les calculs des matrices se fera en 2 parties: Reelle et Complexe.
        //
        /*
         Reelle
         */
        System.out.println("Partie Reelle");
        double tmp[] = new double[nDofCond + nDofDielec];
        System.arraycopy(invSigma[0], 0, tmp, 0, nDofCond + nDofDielec);
        Matrix supp = new Matrix(1, tmp.length, tmp);
        Cell supportConductivite = new Cell(this.Omega.getElementSet());
        RealScalarCellQuantity iSigma = new RealScalarCellQuantity(supportConductivite, supp);
        // Calcul matrice element finis (integral(Wi.Wj))
        FiniteElementFormulation FE = new FiniteElementFormulation(Omega);
        // Nombre de point de gauss pour integration element finis 4  minimum  (on prends les pt de gauss des sources)
        FE.assembly(this.nbPtGSourcesInduc, iSigma);
        StorageSparse SFE = (StorageSparse) FE.getStore();
        this.R[0] = SFE.getMatrixPrecond(null);
        /*
         Complexe
         */
        System.out.println("Partie Complexe");
        System.arraycopy(invSigma[1], 0, tmp, 0, nDofCond + nDofDielec);
        supp = new Matrix(1, tmp.length, tmp);
        iSigma = new RealScalarCellQuantity(supportConductivite, supp);
        // Calcul matrice element finis (integral(Wi.Wj))
        FE.assembly(this.nbPtGSourcesInduc, iSigma);
        SFE = (StorageSparse) FE.getStore();
        this.R[1] = SFE.getMatrixPrecond(null);
    }

    /**
     * Integration de la matrice C
     *
     * @param f la frequence
     */
    protected void integrationC(double f) {
        double omega = 2 * Math.PI * f;

        // Compute 1 / sigmaStar = conjuguee(sigmaStar) / norm(sigmaStar)^2
        double invSigma[][] = new double[2][nDofCond + nDofDielec];
        // Conductor
        double norm = Math.hypot(sigmaStar[0][0], omega * sigmaStar[1][0]);
        double Re = -sigmaStar[0][0] / (omega * norm * norm), Im = -sigmaStar[1][0] / (norm * norm);
        for (int i = 0; i < nDofCond; i++) {
            // La partie reelle 
            invSigma[0][i] = Re;
            //La partie complexe
            invSigma[1][i] = Im;
        }
        // Dielectric
        norm = Math.hypot(sigmaStar[0][1], omega * sigmaStar[1][1]);
        Re = -sigmaStar[0][1] / (omega * norm * norm);
        Im = -sigmaStar[1][1] / (norm * norm);
        for (int i = nDofCond; i < nDofCond + nDofDielec; i++) {
            // La partie reelle 
            invSigma[0][i] = Re;
            //La partie complexe
            invSigma[1][i] = Im;
        }

        //
        // Les calculs des matrices se fera en 2 parties: Reelle et Complexe.
        //
        /*
         Reelle
         */
        System.out.println("Partie Reelle");
        double tmp[] = new double[nDofCond + nDofDielec];
        System.arraycopy(invSigma[0], 0, tmp, 0, nDofCond + nDofDielec);
        Matrix supp = new Matrix(1, tmp.length, tmp);
        Cell supportConductivite = new Cell(this.Omega.getElementSet());
        RealScalarCellQuantity iSigma = new RealScalarCellQuantity(supportConductivite, supp);
        // Calcul matrice element finis (integral(Wi.Wj))
        FiniteElementFormulation FE = new FiniteElementFormulation(Omega);
        // Nombre de point de gauss pour integration element finis 4  minimum  (on prends les pt de gauss des sources)
        FE.assembly(this.nbPtGSourcesInduc, iSigma);
        StorageSparse SFE = (StorageSparse) FE.getStore();
        this.C[0] = SFE.getMatrixPrecond(null);
        /*
         Complexe
         */
        System.out.println("Partie Complexe");
        System.arraycopy(invSigma[1], 0, tmp, 0, nDofCond + nDofDielec);
        supp = new Matrix(1, tmp.length, tmp);
        iSigma = new RealScalarCellQuantity(supportConductivite, supp);
        // Calcul matrice element finis (integral(Wi.Wj))
        FE.assembly(this.nbPtGSourcesInduc, iSigma);
        SFE = (StorageSparse) FE.getStore();
        this.C[1] = SFE.getMatrixPrecond(null);
    }

    /**
     * Integration de la matrice L
     *
     * @param f la frequence
     */
    protected void integrationL(double f) {
        double omega = 2 * Math.PI * f;
        double cst[][] = new double[2][nDofCond + nDofDielec];
        // Compute mu_0 * ( sigmaStar - j omega epsilon_0 ) / (4 * pi * sigmaStar ) = a + j b
        /*
         Partie Conductrice
         */
        double norm = Math.hypot(sigmaStar[0][0], omega * sigmaStar[1][0]);
        // a = mu_0 * omega * sigma * omega * epsilon_0 / norm^2
        double Re = omega * mu0 * sigmaStar[0][0] * omega * eps0 / (norm * norm);
        // b = mu_0 * omega * ( sigma^2 - omega^2 * epsilon * ( epsilon_0 - epsilon ) ) / norm^2
        double Im = omega * mu0 * (sigmaStar[0][0] * sigmaStar[0][0] - omega * omega * sigmaStar[1][0] * (eps0 - sigmaStar[1][0])) / (norm * norm);
        for (int i = 0; i < nDofCond; i++) {
            // La partie reelle 
            cst[0][i] = Re;
            //La partie complexe
            cst[1][i] = Im;
        }
        /*
         Partie Dielectric
         */
        norm = Math.hypot(sigmaStar[0][1], omega * sigmaStar[1][1]);
        // a = sigma * omega * epsilon_0 / norm^2
        Re = omega * mu0 * sigmaStar[0][1] * omega * eps0 / (norm * norm);
        // b = ( sigma^2 - omega^2 * epsilon * ( epsilon_0 - epsilon ) ) / norm^2
        Im = omega * mu0 * (sigmaStar[0][1] * sigmaStar[0][1] - omega * omega * sigmaStar[1][1] * (eps0 - sigmaStar[1][1])) / (norm * norm);
        for (int i = nDofCond; i < nDofCond + nDofDielec; i++) {
            // La partie reelle 
            cst[0][i] = Re;
            //La partie complexe
            cst[1][i] = Im;
        }
        /**
         ***********************************************************************
         * ASSEMBLAGE DE LA MATRICE REELLE
         * **********************************************************************
         */
        System.out.println("Partie Reelle");
        double tmp[] = new double[nDofCond + nDofDielec];
        System.arraycopy(cst[0], 0, tmp, 0, nDofCond + nDofDielec);
        Matrix supp = new Matrix(1, tmp.length, tmp);
        Cell supportConductivite = new Cell(this.Omega.getElementSet());
        RealScalarCellQuantity iSigma = new RealScalarCellQuantity(supportConductivite, supp);
        // Calcul matrice integral (integral (Wi.wj/r)) avec 0 pour la singularité
        GalerkinIntegralFormulation IV = new GalerkinIntegralFormulationFull(Omega, Omega, new MultGvect(), new SelfElementFixedGauss(this.nbPtGSourcesInduc,new Cancel()),this.nbPtGCiblesInduc);
        IV.assembly(iSigma);
        this.L[0] = ((StorageFull) IV.getStore()).getMatrix();
        /**
         ***********************************************************************
         * ASSEMBLAGE DE LA MATRICE COMPLEXE
         * **********************************************************************
         */
        System.out.println("Partie Complexe");
        System.arraycopy(cst[0], 0, tmp, 0, nDofCond + nDofDielec);
        supp = new Matrix(1, tmp.length, tmp);
        supportConductivite = new Cell(this.Omega.getElementSet());
        iSigma = new RealScalarCellQuantity(supportConductivite, supp);
        // Calcul matrice integral (integral (Wi.wj/r)) avec 0 pour la singularité
        IV.assembly(iSigma);
        this.L[1] = ((StorageFull) IV.getStore()).getMatrix();
    }

    /**
     * Integration de la matrice P
     *
     * @param f la frequence
     */
    protected void integrationP(double f) {
        double omega = 2 * Math.PI * f;

        double normCond = Math.hypot(sigmaStar[0][0], omega * sigmaStar[1][0]);
        normCond *= normCond;
        double normDielec = Math.hypot(sigmaStar[0][1], omega * sigmaStar[1][1]);
        normDielec *= normDielec;
        double normAir = Math.abs(omega * eps0); // La valeur absolue est decorative !
        normAir *= normAir;
        // Bonne chance, ou fais moi confiance.
        // Conducteur - dielectrique
        double ReCondDielec = (omega * eps0 / (normCond * normDielec)) * (sigmaStar[0][0] * sigmaStar[0][0] * sigmaStar[0][1] - sigmaStar[0][0] * sigmaStar[0][1] * sigmaStar[0][1]
                + omega * omega * (sigmaStar[1][0] * sigmaStar[1][0] * sigmaStar[0][1] + sigmaStar[1][1] * sigmaStar[1][1] * sigmaStar[0][0]));
        double ImCondDielec = (omega * eps0 / (normCond * normDielec)) * omega * (sigmaStar[1][0] * sigmaStar[0][1] * sigmaStar[0][1] - sigmaStar[1][1] * sigmaStar[0][0] * sigmaStar[0][0]
                + omega * omega * sigmaStar[1][0] * sigmaStar[1][1] * sigmaStar[1][1] - omega * omega * sigmaStar[1][0] * sigmaStar[1][0] * sigmaStar[1][1]);
        // Dielectrique - air
        double ReDielecAir = (omega * eps0 / (normCond * normDielec)) * omega * omega * eps0 * sigmaStar[0][0];
        double ImDielecAir = (omega * eps0 / (normCond * normDielec)) * omega * (-eps0 * sigmaStar[0][1] * sigmaStar[0][1]
                + omega * omega * sigmaStar[1][1] * eps0 * eps0 - omega * omega * sigmaStar[1][1] * sigmaStar[1][1] * eps0);
        /**
         ***********************************************************************
         * ASSEMBLAGE DE LA MATRICE REELLE
         * **********************************************************************
         */
        P[0] = new Matrix(GammaExt.getActiveDofCount(), GammaExt.getActiveDofCount() + GammaInt.getActiveDofCount());
        GalerkinIntegralFormulationFull IFCellexterne = new GalerkinIntegralFormulationFull(this.GammaExt, this.GammaExt, new MultG(), new SelfElementFixedGauss(nbPtGSourcesCapa,new AnalyticalCorrection()),nbPtGCiblesCapa);
        IFCellexterne.assembly();
        P[0].setSubmatrix(0, GammaInt.getActiveDofCount(), ((StorageFull) IFCellexterne.getStore()).getMatrix());

        GalerkinIntegralFormulationFull IFCellinterne = new GalerkinIntegralFormulationFull(this.GammaExt, this.GammaInt, new MultG(), new SelfElementFixedGauss(nbPtGSourcesCapa,new AnalyticalCorrection()),nbPtGCiblesCapa);
        IFCellinterne.assembly();
        P[0].setSubmatrix(0, 0, ((StorageFull) IFCellinterne.getStore()).getMatrix());

        // Copie
        P[1] = P[0].copy();

        // La partie reelle et exterieur (dielec-air)
        P[0].submatrix(0, GammaInt.getActiveDofCount(), GammaExt.getActiveDofCount(), GammaExt.getActiveDofCount()).scale(ReCondDielec);
        // La partie reelle et interieur (cond-dielec)
        P[0].submatrix(0, 0, GammaExt.getActiveDofCount(), GammaInt.getActiveDofCount()).scale(ReDielecAir);

        // La partie imaginaire et exterieur (dielec-air)
        P[1].submatrix(0, GammaInt.getActiveDofCount(), GammaExt.getActiveDofCount(), GammaExt.getActiveDofCount()).scale(ImCondDielec);
        // La partie imaginaire et interieur (cond-dielec)
        P[1].submatrix(0, 0, GammaExt.getActiveDofCount(), GammaInt.getActiveDofCount()).scale(ImDielecAir);

        // Table des surfaces interrieures
        double[] si = new double[GammaInt.getActiveDofCount()];
        for (int i = 0; i < GammaInt.getActiveDofCount(); i++) {
            si[i] = ((GeomElementSurf) GammaInt.getElementSet().getElements(i).getGeom()).getSurface();
        }
        // Table des surfaces exterieures
        double[] se = new double[GammaExt.getActiveDofCount()];
        for (int i = 0; i < GammaExt.getActiveDofCount(); i++) {
            se[i] = ((GeomElementSurf) GammaExt.getElementSet().getElements(i).getGeom()).getSurface();
        }

        for (int i = 0; i < GammaExt.getActiveDofCount(); i++) {
            // Partie Int/Ext
            for (int j = 0; j < GammaInt.getActiveDofCount(); j++) {
                // La partie reelle et complexe subissent le même traitement
                P[0].setElement(i, j, P[0].getElement(i, j) / (si[j] * se[i]));
                P[1].setElement(i, j, P[1].getElement(i, j) / (si[j] * se[i]));
            }
            // Partie Ext/Ext            
            for (int j = GammaInt.getActiveDofCount(); j < GammaInt.getActiveDofCount() + GammaExt.getActiveDofCount(); j++) {
                // La partie reelle et complexe subissent le même traitement
                P[0].setElement(i, j, P[0].getElement(i, j) / (se[j - GammaInt.getActiveDofCount()] * se[i]));
                P[1].setElement(i, j, P[1].getElement(i, j) / (se[j - GammaInt.getActiveDofCount()] * se[i]));
            }
        }
    }

    /**
     * Integration of all the matrices
     *
     * @param f
     */
    public void integration(double f) {

        System.out.println("*** Integration ***");

        // ========================================================================================= //
        // ========================================================================================= //
        // ========================================================================================= //
        System.out.println("Calcul des resistances (matrice element finis)");
        long t0aFEM = System.currentTimeMillis();
        integrationR(f);
        long t1aFEM = System.currentTimeMillis();
        System.out.println(" Temps d'intégration matrice Eléments finis : " + (t1aFEM - t0aFEM) + " ms");
        System.out.println(" ");

        // ========================================================================================= //
        // ========================================================================================= //
        // ========================================================================================= //
        System.out.println("Calcul des inductances (matrice methode integrale)");
        long t0aIV = System.currentTimeMillis();
        integrationL(f);
        long t1aIV = System.currentTimeMillis();
        System.out.println(" Temps d'intégration matrice methode integrale: " + (t1aIV - t0aIV) + " ms");
        System.out.println(" ");

        // ========================================================================================= //
        // ========================================================================================= //
        // ========================================================================================= //
        System.out.println("Calcul des capacites (matrice methode integrale)");
        long t2aIV = System.currentTimeMillis();
        integrationP(f);
        long t3aIV = System.currentTimeMillis();
        System.out.println(" Temps d'intégration matrice methode integrale: " + (t3aIV - t2aIV) + " ms");
        System.out.println(" ");

        // ========================================================================================= //
        // ========================================================================================= //
        // ========================================================================================= //
        System.out.println("Calcul des capacites (matrice element finis)");
        t0aFEM = System.currentTimeMillis();
        integrationC(f);
        t1aFEM = System.currentTimeMillis();
        System.out.println(" Temps d'intégration matrice Eléments finis C : " + (t1aFEM - t0aFEM) + " ms");
        System.out.println(" ");

    }

    /**
     * Affectation du nombre de points de Gauss sur les sources et les cibles
     *
     * @param nbPtsGaussCibles Nombre de poinbts de Gauss sur les cibles
     * @param nbPtsGaussSources Nombre de points de Gauss sur les sources
     */
    public void setPtsDeGaussInductifs(int nbPtsGaussCibles, int nbPtsGaussSources) {
        this.nbPtGCiblesInduc = nbPtsGaussCibles;
        this.nbPtGSourcesInduc = nbPtsGaussSources;
    }

    /**
     * Affectation du nombre de points de Gauss sur les sources et les cibles
     *
     * @param nbPtsGaussCibles Nombre de poinbts de Gauss sur les cibles
     * @param nbPtsGaussSources Nombre de points de Gauss sur les sources
     */
    public void setPtsDeGaussCapacitifs(int nbPtsGaussCibles, int nbPtsGaussSources) {
        this.nbPtGCiblesCapa = nbPtsGaussCibles;
        this.nbPtGSourcesCapa = nbPtsGaussSources;
    }

    @Override
    public int[][] getTopologie() {
        if (this.Omega == null) {
            throw new UnsupportedOperationException("Aucun espace fonctionnel d'associe, employer la methode setMesh en amont");
        }
        // Construction de la matrice topologique
        int matTopo[][] = new int[this.Omega.getActiveDofCount()][4];

        // Indice de noeud libre
        int indNoeudInfinite = this.Omega.getElementSet().getNbElement() + 1;
        // Recuperation des liens entre face et elements
        int lienFaceElement[][] = this.Omega.getFacesSet().getPosElementFace();

        for (int i = 0; i < this.Omega.getActiveDofCount(); i++) {
            // Test si l'on est sur une face exterieure active
            if (lienFaceElement[i][0] == -1 || lienFaceElement[i][1] == -1) {
                if (Omega.getImplicitConstraintAdresses()[0] <= i && i < Omega.getImplicitConstraintAdresses()[1]) {
                    // Dans la borne 1
                    matTopo[i][0] = i; // Numero de branche = numero de face
                    matTopo[i][1] = lienFaceElement[i][0] == -1 ? lienFaceElement[i][1] : lienFaceElement[i][0]; // Numero du noeud interieur
                    matTopo[i][2] = indNoeudInfinite + 1; // Numero du noeud commun a la borne
                    matTopo[i][3] = 1; // Branche de type Source U
                } else if (Omega.getImplicitConstraintAdresses()[1] <= i && i < Omega.getImplicitConstraintAdresses()[2]) {
                    // Dans la borne 2
                    matTopo[i][0] = i; // Numero de branche = numero de face
                    matTopo[i][1] = lienFaceElement[i][0] == -1 ? lienFaceElement[i][1] : lienFaceElement[i][0]; // Numero du noeud interieur
                    matTopo[i][2] = indNoeudInfinite + 2; // Numero du noeud commun a la borne
                    matTopo[i][3] = 1; // Branche de type Source U
                } else {
                    // On ajoute la partie dans le materiau
                    matTopo[i][0] = i; // Numero de branche = numero de face
                    matTopo[i][1] = lienFaceElement[i][0] == -1 ? lienFaceElement[i][1] : lienFaceElement[i][0]; // Numero du noeud interieur
                    matTopo[i][2] = indNoeudInfinite; // Numero du noeud "infini"
                    matTopo[i][3] = -3; // Branche de type capacitive
                }
                System.out.println("Face exterieur d'indice " + i + " : N1 = " + matTopo[i][1] + " N2 = " + matTopo[i][2]);
            } else {
                // On est sur un dof interne
                matTopo[i][0] = i; // Numero de branche = numero de face
                matTopo[i][1] = lienFaceElement[i][0] < lienFaceElement[i][1] ? lienFaceElement[i][0] : lienFaceElement[i][1]; // Numero du noeud interieur de plus petite valeure (normal allant du plus petit numero d'element au plus grand)
                matTopo[i][2] = lienFaceElement[i][0] > lienFaceElement[i][1] ? lienFaceElement[i][0] : lienFaceElement[i][1]; // Numero du noeud interieur de plus grande valeure (normal allant du plus petit numero d'element au plus grand)
                matTopo[i][3] = -2; // Branche de type inductive
                System.out.println("Face interieure d'indice " + i + " : N1 = " + matTopo[i][1] + " N2 = " + matTopo[i][2]);
            }
        }

        return matTopo;
    }

    @Override
    public int getNbLignes() {
        return Omega.getActiveDofCount();
    }

    @Override
    public double[][] getZbFull(double[][] zb, int indDepLigne, int indFinLigne, int indDebCol, int indFinCol, double f) {
//        // On parcours toutes les lignes (faces donc branches) et on remplie la matrice pleine
//        int nL = 0;
//        double omega = Math.PI * 2 * f;
//
//        for (int i = indDepLigne; i <= indFinLigne; i++) {
//            for (int j = indDebCol; j <= indFinCol; j++) {
//                zb[i][2 * (j)] = this.R.getQuick(i, j);
//                zb[i][2 * (j) + 1] = L.getElement(i, j) * omega;
//            }
//        }
        return zb;
    }

    @Override
    public double[] getValeursSourceI(double f) {
        if (this.SourceI == null) {
            SourceI = new double[this.getNbLignes() * 2];
        }
        return SourceI;
    }

    @Override
    public double[] getValeursSourceU(double f) {
        if (this.SourceU == null) {
            SourceU = new double[2 * this.getNbLignes()];
        }
        return SourceU;
    }

    public void setSourceI(double[] I) {
        SourceI = I;
    }

    public void setSourceU(double[] U) {
        SourceU = U;

    }

    @Override
    public SparseVectorComplex getTermesExplicites(int indL, double f) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public SolveurMaillesIndependantes getSolveurCircuit() {
        return solveurCircuit;
    }

    @Override
    public double[] produit(double[] ib, double[] vb, double f) {
        if (cl == null && Pc == null) {
            long t = System.nanoTime();
            cl = new CombLineaireComplexe(new MatrixComplexPartReIm[]{new MatrixComplexPartReIm(R[0], R[1]),
                new MatrixComplexPartReIm(C[0], C[1]), new MatrixComplexPartReIm(new GOTMatrix(L[0]), new GOTMatrix(L[1]))});
            Pc = new MatrixComplexPartReIm(new GOTMatrix(P[0]), new GOTMatrix(P[1]));
            System.out.println("Time to compute the CL= " + (System.nanoTime() - t) * 1e-9 + " sec");
        }
        // Extrait le sous vecteur pour la multiplication avec P !
        double subIb[] = new double[2 * (GammaExt.getActiveDofCount() + GammaInt.getActiveDofCount())];
        System.arraycopy(ib, ib.length - subIb.length, subIb, 0, subIb.length);

        // Effectue le 'grand' produit
        vb = cl.product(ib, vb);

        // Vecteur reulstant de la multiplication avec P !
        double subVb[] = new double[2 * GammaExt.getActiveDofCount()];
        subVb = Pc.product(subIb, subVb);

        // Assemble les vecteurs
        int offSet = 2 * (Omega.getActiveDofCount() - GammaExt.getActiveDofCount());
        for (int i = 0; i < 2 * GammaExt.getActiveDofCount(); i++) {
            vb[offSet + i] += subVb[i];
        }

        return vb;
    }

    class CombLineaireComplexe implements ProductComplex {

        MatrixComplexPartReIm z[];

        public CombLineaireComplexe(MatrixComplexPartReIm S[]) {
            z = S;
        }

        @Override
        public double[] product(double[] x, double[] res) {
            if (res == null) {
                res = new double[2 * z[0].getRows()];
            }
            for (MatrixComplexPartReIm z1 : z) {
                double tmp[] = new double[2 * z[0].getRows()];
                tmp = z1.product(x, tmp);
                for (int j = 0; j < tmp.length; j++) {
                    res[j] += tmp[j];
                }
            }
            return res;
        }

    }

    /**
     * Print on a gmsh format !
     *
     * @param fileName
     */
    public void printCellBorder(String fileName) {
        ExportGmsh export = new ExportGmsh(GammaExt.getElementSet(), fileName);
    }

    public static void main(String[] args) throws IOException{                   String defaultPath = null;                  if (args.length != 0) {                          defaultPath = args[0]+args[1];                      new File(defaultPath).mkdirs();         }
        /*
         Nombre de regions importees : 6
         Region 0, Nom : CONDUCTOR, type : 3, 244 elements
         Region 1, Nom : DIELECTRIC, type : 3, 630 elements
         Region 2, Nom : BORNE1, type : 2, 4 elements
         Region 3, Nom : BORNE2, type : 2, 4 elements
         Region 4, Nom : COND_BORDER, type : 2, 164 elements
         Region 5, Nom : DIELEC_BORDER, type : 2, 288 elements
         */
//        String path = "D:\\Meshs\\Dielectric\\Wire\\WIRE.DEC";
        String path = "D:\\Meshs\\Dielectric\\Wire\\WIRE_SMALL.DEC";

        ImportFlux Imp = new ImportFlux(path);
        int nDof[] = new int[6];
        for (int i = 0; i < nDof.length; i++) {
            nDof[i] = Imp.getRegion(i).getElementSet().getNbElement();
        }
        RegionsSet Reg = new RegionsSet(Imp);
        // Organise les faces : Borne1 -> Borne2 -> GammaInt -> GammaExt
        FaceConstraint c[] = new FaceConstraint[]{new ExternalDriveFaceConstraint(Reg.getRegion(2)), new ExternalDriveFaceConstraint(Reg.getRegion(3)), new ExternalDriveFaceConstraint(Reg.getRegion(4)), new ExternalDriveFaceConstraint(Reg.getRegion(5))};

        FaceDeg1 Omega = new FaceDeg1((ElementSetHomogene) Reg.generateUnion(new int[]{0, 1}, "Omega").getElementSet(), c);

        double f = 2e5;
        /*
         CONSTRUCTION DU SUPPORT DES CONDUCTIVITES EQUIVALENTES
         */
        double Sigma[][] = new double[2][2];// Size = [2][nbRegions]
        // Conducteur
        double Re = 59.6 * 1e6, Im = eps0;
        Sigma[0][0] = Re;// Conductivite du cuivre
        Sigma[1][0] = Im;// Epsilon0
        // Dielectrique
        Re = 0;
        Im = 2.25;
        Sigma[0][1] = Re;// 0 pour le dielec
        Sigma[1][1] = Im;// Espilon

        // Compute the border
        Cell gammaInt = new Cell((ElementSetHomogene) Imp.getRegion(4).getElementSet());
        Cell gammaExt = new Cell((ElementSetHomogene) Imp.getRegion(5).getElementSet());
        // Create the solver
        PEEC_RLMPC_WIRE solP = new PEEC_RLMPC_WIRE(Omega, Sigma, gammaExt, gammaInt, nDof);
        solP.setPtsDeGaussInductifs(15, 15);
        solP.setPtsDeGaussCapacitifs(3, 3);

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        int add[] = Omega.getImplicitConstraintAdresses();
        int nbF[] = Omega.getImplicitConstraintCount();
        System.out.println("add= " + add[0] + " , " + add[1]);
        System.out.println("nbF= " + nbF[0] + " , " + nbF[1]);

        int nbBranches = Omega.getActiveDofCount() + 1;
        int C1 = Omega.getElementSet().getNbElement() + 2;
        int C2 = C1 + 1;

        circuitPur.addSourceISimple(nbBranches, C1, C2, "SourceI", 1, 0.0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);
        
        // Resolution
        double ib[][] = solP.resolutionIterative(f);
        Matrix res = new Matrix(2, Omega.getActiveDofCount());
        for (int i = 0; i < Omega.getActiveDofCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i] );
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1] );
        }

        // On reorganise les dof   
        RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
        RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));

        ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, (defaultPath != null ? defaultPath : "Y:")+"/Resultats/Dielectric/WireRE_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJreal.addQuantity(Jreal, "Jreal");

        ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, (defaultPath != null ? defaultPath : "Y:")+"/Resultats/Dielectric/WireIM_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJimag.addQuantity(Jimag, "Jimag");

        ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
        exportJreal = new ExportGmshHdiv(Omega, (defaultPath != null ? defaultPath : "Y:")+"/Resultats/Dielectric/WireJmod" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        exportJreal.addQuantityExportMod(J, "Jmod");

        File out_file = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources/Resultats/Dielectric")  +"/WireJmod" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        File out_file1 = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources/Resultats/Dielectric")  +"/WireRE_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
        File out_file2 = new File((defaultPath != null ? defaultPath : "D:/jsiau/_Backup_Sources/Resultats/Dielectric")  +"/WireIM_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");        
        try {      
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);
        } catch (IOException ex) {
            Logger.getLogger(PEEC_RLMPC_WIRE.class.getName()).log(Level.SEVERE, null, ex);
        }
       

        //GestionnaireTaches.getGestionnaireTaches().stop();
    }
}
