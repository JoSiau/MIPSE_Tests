/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC;

import g2elab.mipse.circuit.solverCircuitComplex.BlocCircuit;
import g2elab.mipse.circuit.solverCircuitComplex.SolveurMaillesIndependantes;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmsh;
import g2elab.mipse.meshCore.elements.*;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.functionSpace.Hdiv;
import g2elab.mipse.meshCore.geomElements.GeomAbstrait;
import g2elab.mipse.meshCore.geomElements.GeomElementSurf;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.meshCore.tree.faceTree.Face;
import g2elab.mipse.meshCore.tree.faceTree.FaceSet;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.NoCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.integrationstrategies.FixedGaussSource;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageSparse;
import g2elab.mipse.numericalTools.matrix.MatriceIncidence;
import g2elab.mipse.numericalTools.matrix.complex.MatrixComplexPartReIm;
import g2elab.mipse.numericalTools.matrix.real.dense.gotMatrix.GOTMatrix;
import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
import g2elab.mipse.numericalTools.vector.sparse.SparseVectorComplex;
import g2elab.mipse.tools.files.Ecriture;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class PEECwithCapa_Full implements BlocCircuit {//, ProductComplex, ComplexPrecond{

    /**
     * Stockage de l'espace fonctionnel
     */
    private FaceDeg1 faceDeg1;
    /**
     * Matrice des resistances
     */
    private SparseMatrixRowReal R;
    /**
     * Matrice d'induction
     */
    private Matrix L;
    /**
     * Matrice de capacité
     */
    private Matrix P;

    /**
     * Matrice imaginaire
     */
    private MatrixComplexPartReIm zb = new MatrixComplexPartReIm();
    /**
     * Nombre de points de gauss source et cible pour l'integration de L
     */
    private int nbPtGCiblesInduc, nbPtGSourcesInduc;
    /**
     * Nombre de points de gauss source et cible pour l'integration de P
     */
    private int nbPtGCiblesCapa, nbPtGSourcesCapa;
    /**
     * L'offset entre la matrice L et P
     */
    private int offSet;
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
    private final double eps0 = 8.85418782 * 1e-12;
    /**
     * the width of the mesh if it is a surf mesh
     */
    private double thickness = 1;
    /**
     * rho
     */
    private double ro = 1 / (3.526 * Math.pow(10, 7)); // Par defaut
    /**
     * The border
     */
    private Cell borderCell = null;
    /**
     * tmp
     */
    private int typeElement = -1;
    /**
     * number of branches
     */
    private int nbBranch = 0;
    /**
     * Control the formulation 1=R,2=R+LM,3=R+LM+P
     */
    private int formulation;
    /**
     * Set the verbosity
     */
    boolean verbose = true;
    private boolean FA_capa;
    private int nbPtsSourcesInducCorr = 0;

//    /**
//     * Default constructor. Use all the functions "set...()" you must do before
//     * solving anything ! May the Force be with you !
//     *
//     * @param FS the function space (must be an instance of Hdiv)
//     * @param indBorder the index of the constraint renumbering (implicit)
//     */
//    public PEECwithCapa_Full(FaceDeg1 FS, int indBorder) {
//        if (!(FS instanceof Hdiv)) {
//            throw new IllegalArgumentException("The function space must be an Hdiv");
//        }
//        faceDeg1 = FS;
//        borderCell = FS.generateCellSpaceFunction(indBorder);
//        solveurCircuit = new SolveurMaillesIndependantes(2);
//        solveurCircuit.setBloc(this, 0, 0);
//    }
    /**
     * Default constructor. Use all the functions "set...()" you must do before
     * solving anything ! May the Force be with you !
     *
     * @param FS the function space (must be an instance of Hdiv)
     */
    public PEECwithCapa_Full(FaceDeg1 FS) {
        if (!(FS instanceof Hdiv)) {
            throw new IllegalArgumentException("The function space must be an Hdiv");
        }
        faceDeg1 = FS;
        computeBorderCell();
        solveurCircuit = new SolveurMaillesIndependantes(2);
        solveurCircuit.setBloc(this, 0, 0);
    }

    /**
     * Default constructor. Use all the functions "set...()" you must do before
     * solving anything ! May the Force be with you !
     *
     * @param FS the function space (must be an instance of Hdiv)
     * @param formulation 0: Resistive ; 1: R+ Inductive ; 2: R+LM+Capcitive
     */
    public PEECwithCapa_Full(FaceDeg1 FS, int formulation) {
        if (!(FS instanceof Hdiv)) {
            throw new IllegalArgumentException("The function space must be an Hdiv");
        }
        if (formulation != 0 && formulation != 1 && formulation != 2) {
            throw new IllegalArgumentException("Wrong arugment ! formulation must be in {1, 2, 3}");
        }
        this.formulation = formulation;
        faceDeg1 = FS;
        // We need the border only for the Capacitive's effect
        if (formulation > 1) {
            computeBorderCell();
        }
        solveurCircuit = new SolveurMaillesIndependantes(2);
        solveurCircuit.setBloc(this, 0, 0);

        faceDeg1.getElementSet().plotElementSet("d:/faceDeg1.vtk");
    }

    /**
     * Resolution directe du probleme
     *
     * @param f Liste des frequence a analyser
     * @return Vecteur contenant les resultats (courants complexe pour chaque
     * frequence puis tension)
     */
    public double[][] resolutionIterative(double f) {

        this.integration();
        this.solveurCircuit.makeCircuit();
        this.solveurCircuit.setConfiguration(new double[]{1, 100, 1e-6, 1, -50}, new double[]{0, 50});
        //*
        return this.solveurCircuit.resolutionIterative(new double[]{f}, true);
        /*/
         double [][]tmp = this.solveurCircuit.resolutionIterative(new double[]{f}, true);
         g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D matRes = this.solveurCircuit.getFullMatrix();
         g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D.geneMatlab(matRes.getArray(), "A", "d:\\");
         return tmp;
         //*/
    }

    public double[][] resolutionDirecte(double f) {
        this.solveurCircuit.makeCircuit();
        // Check if the matrices are already computed
        if (R == null || L == null || P == null) {
            this.integration();
        }
        //this.currentMI.getCircuit().createFileGephi("Y:/CAPAnodes.csv", "Y:/CAPAbranches.csv", ",");
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
     * Set the thickness of the mesh (only used for surf mesh !)
     *
     * @param e the thickness
     */
    public void setThickness(double e) {
        this.thickness = e;
    }

    /**
     * Integration of all the matrices
     */
    public void integration() {

        System.out.println("*** Integration ***");
        System.out.println("Calcul des resistances (matrice element finis)");
        long t0aFEM = System.currentTimeMillis();
        // support pour les conductivites
        Matrix resVal = new Matrix(1, this.faceDeg1.getElementSet().getNbElement());

        resVal.setAllElements(this.ro);
        // Espace focntionnel constant par element
        Cell supportConductivite = new Cell(this.faceDeg1.getElementSet());
        RealScalarCellQuantity invSigma = new RealScalarCellQuantity(supportConductivite, resVal);

        // Calcul matrice element finis (integral(Wi.Wj))
        FiniteElementFormulation FE = new FiniteElementFormulation(faceDeg1);
        // Nombre de point de gauss pour integration element finis 4  minimum  (on prends les pt de gauss des sources)
        FE.assembly(this.nbPtGSourcesInduc, invSigma);
        StorageSparse SFE = (StorageSparse) FE.getStore();
        this.R = SFE.getMatrixPrecond(null);

        if (this.thickness != 1 && faceDeg1.getElementSet().isSurf()) { // Only for a surf mesh
            R.scale(1 / this.thickness);
        }

        long t1aFEM = System.currentTimeMillis();
        System.out.println(" Temps d'intégration matrice Eléments finis : " + (t1aFEM - t0aFEM) + " ms");
        System.out.println(" ");

        if (formulation > 0) {
            System.out.println("Calcul des inductances (matrice methode integrale)");
            long t0aIV = System.currentTimeMillis();
            // Calcul matrice integral (integral (Wi.wj/r)) avec 0 pour la singularité
            GalerkinIntegralFormulation IV;
            if (nbPtsSourcesInducCorr == 0) {
                IV = new GalerkinIntegralFormulationFull(faceDeg1, faceDeg1, new MultGvect(), new NoCorrection(new FixedGaussSource(this.nbPtGSourcesInduc)), this.nbPtGCiblesInduc);
            } else {
                IV = new GalerkinIntegralFormulationFull(faceDeg1, faceDeg1, new MultGvect(), new SelfElementFixedGauss(nbPtGSourcesInduc, new InCreasedPGSourceNumber(nbPtsSourcesInducCorr)), this.nbPtGCiblesInduc);
            }
            IV.assembly();
            this.L = ((StorageFull) IV.getStore()).getMatrix();
            L.scale(mu0);
            long t1aIV = System.currentTimeMillis();
            System.out.println(" Temps d'intégration matrice methode integrale: " + (t1aIV - t0aIV) + " ms");
            System.out.println(" ");
        }

        if (formulation > 1) {
            System.out.println("Calcul des capacites (matrice methode integrale)");
            long t2aIV = System.currentTimeMillis();
            // Compute the matrice [P]
            // L = intBorder 1/Si * 1/Sj ( Wi * 1/r * Wj) avec correction analytique pour influence element sur lui meme
            GalerkinIntegralFormulation IFCell;
            switch (this.typeElement) {
                case 1:
                    IFCell = new GalerkinIntegralFormulationFull(borderCell, borderCell, new MultG(), new SelfElementFixedGauss(nbPtGSourcesCapa, new AnalyticalCorrection()), nbPtGCiblesCapa);
                    break;
                case 2:
                    IFCell = new GalerkinIntegralFormulationFull(borderCell, borderCell, new MultG(), new SelfElementFixedGauss(nbPtGSourcesCapa, new AnalyticalCorrection()), nbPtGCiblesCapa);
                    break;
                default:
                    throw new RuntimeException(" Not implemented yet!");
            }
            IFCell.assembly();
            // On recupere la matrice (un peu laborieux mais bon....)
            P = ((StorageFull) IFCell.getStore()).getMatrix();

            // table des surface
            double[] S = new double[P.getRowCount()];
            for (int i = 0; i < P.getRowCount(); i++) {
                GeomAbstrait geom = borderCell.getElementSet().getElements()[i].getGeom();
                S[i] = ((GeomElementSurf) geom).getSurface();
            }

            borderCell.getElementSet().plotElementSet("d:/borderCell.vtk");

            double invEpsiZero = 1 / eps0;
            for (int i = 0; i < P.getRowCount(); i++) {
                for (int j = 0; j < P.getColumnCount(); j++) {
                    P.setElement(i, j, invEpsiZero * P.getElement(i, j) / S[i] / S[j]);
                }
            }
            long t3aIV = System.currentTimeMillis();
            System.out.println(" Temps d'intégration matrice methode integrale: " + (t3aIV - t2aIV) + " ms");
            System.out.println(" ");

            this.offSet = L.getRowCount() - P.getRowCount();
            System.err.println("offSet= " + offSet + " \t size2= " + P.getRowCount() + "\t size3= " + borderCell.getTotalDofCount());
        }

    }

    private void ecrireMat(Matrix r1, double s, String name) {
        System.out.println("Ecriture de la matrix (" + r1.getRowCount() + " , " + r1.getColumnCount() + ")");
        Ecriture e1 = null;
        try {
            e1 = new Ecriture(name);
            double tab[][] = new double[r1.getRowCount()][r1.getColumnCount()];
            for (int i = 0; i < r1.getRowCount(); i++) {
                for (int j = 0; j < r1.getColumnCount(); j++) {
                    tab[i][j] = r1.getElement(i, j) * s;
                }
            }
            e1.ecrire(tab, ' ');
            e1.close();
        } catch (IOException ex) {
            Logger.getLogger(LoopAntenna.class.getName()).log(Level.SEVERE, null, ex);
        }
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
     * @param nbPtsSourcesCorr
     */
    public void setPtsDeGaussInductifs(int nbPtsGaussCibles, int nbPtsGaussSources, int nbPtsSourcesCorr) {
        this.nbPtGCiblesInduc = nbPtsGaussCibles;
        this.nbPtGSourcesInduc = nbPtsGaussSources;
        this.nbPtsSourcesInducCorr = nbPtsSourcesCorr;
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

    /**
     * Affecte la valeur de \rho (constante devant [R]).
     *
     * @param rho a double
     */
    public void setRho(double rho) {
        this.ro = rho;
    }

    /**
     * Return the matrix L
     *
     * @return
     */
    public Matrix getL() {
        return this.L;
    }

    public FaceDeg1 getFunctionSpace() {
        return this.faceDeg1;
    }

    public Cell getBorderCell() {
        return borderCell;
    }

    public SparseMatrixRowReal getR() {
        return R;
    }

    @Override
    public int[][] getTopologie() {
        if (this.faceDeg1 == null) {
            throw new UnsupportedOperationException("Aucun espace fonctionnel d'associe, employer la methode setMesh en amont");
        }
        // Construction de la matrice topologique
        int matTopo[][] = new int[this.faceDeg1.getActiveDofCount() + this.getSizeC()][4];

        // Indice de noeud libre
        int indNoeudLibre = this.faceDeg1.getElementSet().getNbElement() + 1;
        // Recuperation des liens entre face et elements
        int lienFaceElement[][] = this.faceDeg1.getFacesSet().getPosElementFace();

        int type = -2;
        if (formulation == 0) { // Si la formulation est uniquement resistive        
            type = -1;
        }

        // On parcours toutes les faces et l'on construit la topologie
        // Le numero des face correspond au numero des branches
        // Le numero des elements correspond au numero des noeuds
        // Les normales allant du plus petit element vers le plus grand (en numero d'identification) les branche sont orientés de la meme facon afin d'avoir la coherence entre les signatures des matrices integrees et les convention du solveur circuit
        int compteur = 0;
        for (int i = 0; i < this.faceDeg1.getActiveDofCount(); i++) {
            // Test si l'on est sur une face exterieure active
            if (lienFaceElement[i][0] == -1 || lienFaceElement[i][1] == -1) {
                // On ajoute la partie dans le materiau
                matTopo[i][0] = i; // Numero de branche = numero de face
                matTopo[i][1] = lienFaceElement[i][0] == -1 ? lienFaceElement[i][1] : lienFaceElement[i][0]; // Numero du noeud interieur
                matTopo[i][2] = indNoeudLibre + compteur; // Numero du noeud sur la face
                compteur = compteur + 1;
                matTopo[i][3] = type; // Branche de type inductive
                if (verbose) {
                    System.out.println("Face exterieur d'indice " + i + " : N1 = " + matTopo[i][1] + " N2 = " + matTopo[i][2]);
                }
            } else {
                // On est sur un dof interne
                matTopo[i][0] = i; // Numero de branche = numero de face
                matTopo[i][1] = lienFaceElement[i][0] < lienFaceElement[i][1] ? lienFaceElement[i][0] : lienFaceElement[i][1]; // Numero du noeud interieur de plus petite valeure (normal allant du plus petit numero d'element au plus grand)
                matTopo[i][2] = lienFaceElement[i][0] > lienFaceElement[i][1] ? lienFaceElement[i][0] : lienFaceElement[i][1]; // Numero du noeud interieur de plus grande valeure (normal allant du plus petit numero d'element au plus grand)
                matTopo[i][3] = type; // Branche de type inductive
                if (verbose) {
                    System.out.println("Face interieure d'indice " + i + " : N1 = " + matTopo[i][1] + " N2 = " + matTopo[i][2]);
                }
            }
        }

        // Add the capacitive faces on the border
        if (this.formulation > 1) {
            //Ajoute la partie capa
            int offSet = faceDeg1.getActiveDofCount();
            int nbCapa = 0; // Count the number of capacitive faces
            int noeudCapaCommun = indNoeudLibre + compteur; // The last index of the inductives faces

            //noeudCapaCommun = borderCell.getActiveDofCount(); 
            compteur = 0;
            if (faceDeg1.getElementSet().isSurf()) {
                ElementSetHomogene elementSet = faceDeg1.getElementSet();
                for (int i = 0; i < elementSet.getNbElement(); i++) {
                    int numNode = elementSet.getElements(i).getNum();
                    matTopo[i + offSet][0] = offSet + nbCapa;
                    matTopo[i + offSet][1] = numNode;
                    matTopo[i + offSet][2] = noeudCapaCommun; // Numero du noeud commun      
                    matTopo[i + offSet][3] = -3; // Branche de type capacitive         
                    if (verbose) {
                        System.out.println("Face interieure CAPA d'indice " + (offSet + i) + " : N1 = " + matTopo[i + offSet][1] + " N2 = " + matTopo[i + offSet][2]);
                    }
                    nbCapa++;
                }
            } else {
                for (int i = 0; i < this.faceDeg1.getActiveDofCount(); i++) {
                    // Test si l'on est sur une face exterieure active
                    if (lienFaceElement[i][0] == -1 || lienFaceElement[i][1] == -1) {
                        // On ajoute la partie dans le materiau
                        matTopo[nbCapa + offSet][0] = offSet + nbCapa; // Numero de branche 
                        matTopo[nbCapa + offSet][1] = indNoeudLibre + compteur; // Numero du noeud interieur
                        matTopo[nbCapa + offSet][2] = noeudCapaCommun; // // Numero du noeud commun                    
                        matTopo[nbCapa + offSet][3] = -3; // Branche de type capacitive
                        if (verbose) {
                            System.out.println("Face exterieure CAPA d'indice " + (offSet + nbCapa) + " : N1 = " + matTopo[nbCapa + offSet][1] + " N2 = " + matTopo[nbCapa + offSet][2]);
                        }
//System.out.println(" facet Numero : "+j + " lien FaceElement 0 1 : " + lienFaceElement[j][0]+ " "+ lienFaceElement[j][1]);
                        nbCapa++;
                        compteur++;
                    }
                }
            }
        }
        nbBranch = matTopo.length;
        return matTopo;
    }

    public Matrix getP() {
        return P;
    }

    public int getSizeL() {
        return this.faceDeg1.getActiveDofCount();
    }

    public int getSizeC() {
        if (this.formulation > 1) {
            return this.borderCell.getActiveDofCount();
        } else {
            return 0;
        }
    }

    @Override
    public int getNbLignes() {
        if (this.faceDeg1 == null) {
            throw new UnsupportedOperationException("Aucun espace fonctionnel d'associe, employer la methode setMesh en amont");
        }
        return (getSizeL() + getSizeC());
    }

    @Override
    public double[][] getZbFull(double[][] zb, int indDepLigne, int indFinLigne, int indDebCol, int indFinCol, double f) {
        // On parcours toutes les lignes (faces donc branches) et on remplie la matrice pleine
        int nL = getSizeL();
        double omega = Math.PI * 2 * f;

        for (int i = indDepLigne; i <= indFinLigne; i++) {
            for (int j = indDebCol; j <= indFinCol; j++) {
                if (i < nL && j < nL) {
                    zb[i][2 * (j)] = this.R.getQuick(i, j);
                    zb[i][2 * (j) + 1] = L.getElement(i, j) * omega;
                } else if (i >= nL && j >= nL) {
                    zb[i][2 * (j) + 1] = -this.P.getElement(i - nL, j - nL) / omega;
                }
            }
        }

//        ecrireMat(L, omega, "d:/L.out");
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

    /**
     * Control if it print the informations of topology
     *
     * @param verb
     */
    public void setVerbose(boolean verb) {
        this.verbose = verb;
    }

    public void switch2FullAnalytical4Capa() {
        this.FA_capa = true;
    }

    public void setSourceI(double[] I) {
        SourceI = I;
    }

    public void setSourceU(double[] U) {
        SourceU = new double[2 * (getSizeC() + getSizeL())];

        System.arraycopy(U, 0, SourceU, 0, 2 * getSizeL());

        if (this.formulation > 1) {
            System.arraycopy(U, 2 * (getSizeL() - getSizeC()), SourceU, 2 * getSizeL(), 2 * getSizeC());
        }
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

        int sizeL = this.getSizeL();
        int sizeC = this.getSizeC();

        double LRe[] = new double[sizeL];
        double LIm[] = new double[sizeL];
        double CRe[] = new double[sizeC];
        double CIm[] = new double[sizeC];

        for (int i = 0; i < ib.length / 2; i++) {
            if (i < sizeL) {
                LRe[i] = ib[2 * i];
                LIm[i] = ib[2 * i + 1];
            } else {
                CRe[i - sizeL] = ib[2 * i];
                CIm[i - sizeL] = ib[2 * i + 1];
            }
        }

        double[] zbRLRe = new double[LRe.length];
        double[] zbRLIm = new double[LIm.length];
        // Resistives
        zbRLRe = R.product(LRe, zbRLRe);
        zbRLIm = R.product(LIm, zbRLIm);

        double[] zbLLRe = new double[LRe.length];
        double[] zbLLIm = new double[LIm.length];
        if (this.formulation > 0) { // Inductives
            GOTMatrix zbLprod = new GOTMatrix(L);
            zbLLRe = zbLprod.product(LRe, zbLLRe);
            zbLLIm = zbLprod.product(LIm, zbLLIm);
        }
        double[] zbCCRe = new double[CRe.length];
        double[] zbCCIm = new double[CIm.length];
        if (this.formulation > 1) { // Capacitives
            GOTMatrix zbCprod = new GOTMatrix(P);
            zbCCRe = zbCprod.product(CRe, zbCCRe);
            zbCCIm = zbCprod.product(CIm, zbCCIm);
        }
        double omega = 2 * Math.PI * f;
        // Assemblee the result
        for (int i = 0; i < vb.length / 2; i++) {
            if (i < sizeL) {
                vb[2 * i] = zbRLRe[i] - omega * zbLLIm[i];
                vb[2 * i + 1] = zbRLIm[i] + omega * zbLLRe[i];
            } else {
                vb[2 * i] = zbCCIm[i - sizeL] / omega;
                vb[2 * i + 1] = -zbCCRe[i - sizeL] / omega;
            }
        }
        return vb;
    }

//    @Override
//    public double[] product(double[] x, double[] mult) {
//        // MtI = M^t . x
//        double MtI[] = M.produitTransposeVecteurComplexe(x);
//        // Creer un stockage complexe
//        zb = new MatrixComplexPartReIm();
//        zb.setPartReIm(zbRe, new StorageCombLin(1, this.L, 1, this.P, offSet, offSet));
//        // ZbMtI = Z . MtI
//        double ZbMtI[] = zb.product(MtI, null);
//        // Return [M] [Z] [M^t] x
//        return M.produitVecteurComplexe(ZbMtI);
//    }
//
//    // PRECONDITIONNED WITH THE IDENTITY ! WHAT A JOKE
//    @Override
//    public double[] precond(double[] x, double[] mult) {
//        /*
//         mult = this.M.produitTransposeVecteurComplexe(x);
//         double tmp[] = precondDiag.precond(nbUnknown, mult, null);
//         mult = this.M.produitVecteurComplexe(tmp);
//         return mult;
//         /*/
//        if (mult == null) {
//            mult = new double[x.length];
//        }
//        System.arraycopy(x, 0, mult, 0, x.length);
//        return mult;
//        //*/
//    }
//
//    @Override
//    public void configuration(double[] config) {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//
//    @Override
//    public long getMemoryUsed() {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
//
//    @Override
//    public void setMatrix(Object mat, double[] configuration) {
//        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
//    }
    public MatriceIncidence getMI() {
        return M;
    }

    /**
     * Compute a cell that contains all the elements of the border of faceDeg1
     */
    protected final void computeBorderCell() {
        System.out.println(".......... Border cell creation ..........");
        if (faceDeg1.getElementSet().isVol()) {            
            borderCell = faceDeg1.generateCellSpaceFunction(0, faceDeg1.getImplicitConstraintAdresses().length-1);            
        } else if (faceDeg1.getElementSet().isSurf()) {
            borderCell = new Cell(faceDeg1.getElementSet());
            if (borderCell.getElementSet().getElements(0) instanceof TriangleDroit) {
                this.typeElement = 1;
            } else if (borderCell.getElementSet().getElements(0) instanceof QuadrangleDroit) {
                this.typeElement = 2;
            } else {
                throw new IllegalArgumentException("Only triangles and quadrangles are supported !");
            }
        } else {
            throw new IllegalArgumentException("Only volumic and surf meshes are supported !");
        }

        System.out.println(" ... End of border cell creation !!! ");
    }

    /**
     * Print on a gmsh format !
     *
     * @param fileName
     */
    public void printCellBorder(String fileName) {
        ExportGmsh export = new ExportGmsh(borderCell.getElementSet(), fileName);
    }
}
