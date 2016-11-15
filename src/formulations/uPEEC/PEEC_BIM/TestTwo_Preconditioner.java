/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.PEEC_BIM;

import g2elab.mipse.circuit.equationsCircuit.Circuit;
import g2elab.mipse.circuit.solverCircuitComplex.BlocCircuit;
import g2elab.mipse.circuit.solverCircuitComplex.SolveurMaillesIndependantes;
import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Surf.ExportTopologieGMSH;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshCell;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.contraints.ExternalDriveCellConstraint;
import g2elab.mipse.meshCore.contraints.ExternalDriveFaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceRealDirichlet;
import g2elab.mipse.meshCore.elements.ElementLinSetHomogene;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.elements.Node;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.functionSpace.NormalCell;
import g2elab.mipse.meshCore.geomElements.GeomElementSurf;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.meshCore.region.LineRegion;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.mipseCore.blockAssembly.real.RealMatrixBlocScale;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.*;
import g2elab.mipse.mipseCore.integralIntegration.integrationstrategies.FixedGaussSource;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.integralIntegration.kernel.NegMultDG;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.Storage;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import g2elab.mipse.mipseCore.storage.StorageSparse;
import g2elab.mipse.numericalTools.iterativeSolver.ProductComplex;
import g2elab.mipse.numericalTools.iterativeSolver.complex.FGMResComplex;
import g2elab.mipse.numericalTools.matrix.MatriceIncidence;
import g2elab.mipse.numericalTools.matrix.complex.ComplexOperator;
import g2elab.mipse.numericalTools.matrix.complex.MatrixBlocComplex;
import g2elab.mipse.numericalTools.matrix.complex.bloc.ComplexMatrixBloc;
import g2elab.mipse.numericalTools.matrix.complex.bloc.ComplexMatrixBlocAdd;
import g2elab.mipse.numericalTools.matrix.complex.bloc.ComplexMatrixBlocScale;
import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.SolverLUBasic2D;
import g2elab.mipse.numericalTools.matrix.real.dense.basic2D.Basic2D;
import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
import g2elab.mipse.numericalTools.preconditioner.ComplexPrecond;
import g2elab.mipse.numericalTools.vector.sparse.SparseVectorComplex;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.list.IntArrayList;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.Matrix;
import got.matrix.RowVector;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.logging.Level;
import java.util.logging.Logger;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.eps0;
import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.CONSTANTS.mu0;

/**
 *
 * @author jsiau
 */
public class TestTwo_Preconditioner implements BlocCircuit, ProductComplex, ComplexPrecond {

    /**
     * Use the Hmatrix compression or not for the integral matrices.
     */
    private boolean HmatCompression = false;
    private double epsHmat = 1e-4, eta = 2.0;
    private int nmin = 30, kmax = 50, order = 3;
    private boolean recomp = true;

    private final FaceDeg1 FDcond;
    private final Cell CellDielec;
    private final boolean borne;
    private final Cell CondExt;
    private final SolveurMaillesIndependantes solveurCircuit;
    private final double sigma;
    private final double ep;
    private final double epsilon;

    /**
     * Nombre de points de gauss source et cible pour l'integration de L.
     */
    private int nbPtGCiblesInduc, nbPtGSourcesInduc, nbPtGSourcesInducCorr = -1;
    /**
     * Nombre de points de gauss source et cible pour l'integration de P.
     */
    private int nbPtGCiblesCapa, nbPtGSourcesCapa;
    private SparseMatrixRowReal R;
    private Matrix L;
    private StorageHmatrix Lh;
    private double[] DielecInvSurf;
    private double[] CondInvSurf;
    private Storage Pc;
    private Storage B;
    private Storage A;
    private Storage Pd;
    private double[] Q;
    private double fcurr = -1;
    private MatrixBlocComplex Zm;
    private ProductComplex Z;
    private double[][] Zfull;
    private String exportTopo = null;
    private int[] equiPot = null;
    private boolean fullCapa = false;
    private boolean voisCorr = false;
    private final int indBornesCons;
    private final boolean isSurf;
    private boolean computePc = true;

    /**
     * Main constructor 101.
     *
     * @param conductor conductor region,
     * @param sigma its conductivity,
     * @param thickness its thickness,
     * @param posFlux positive terminal to plug a source (can be null)
     * @param negFlux negative terminal to plug a source (can be null)
     * @param dielectric the dielectric border,
     * @param epsilon dielectric's constant,
     * @param regCondDielec the region that would be at the intersection of the
     * dielectrric and the conductor, it must be on the conductor.<BR>
     * #WARNING: An inside (for the dielectric) node will be searched on the
     * conductor region. So it may fail if the conductor exceed from the
     * dielectric.
     */
    public TestTwo_Preconditioner(SurfaceRegion conductor, double sigma, double thickness, LineRegion posFlux, LineRegion negFlux,
            SurfaceRegion dielectric, double epsilon,
            SurfaceRegion regCondDielec) {
        this(conductor, sigma, thickness, posFlux, negFlux, dielectric, epsilon, regCondDielec, null, false);
    }

    /**
     * Main constructor 201.
     *
     * @param conductor conductor region,
     * @param sigma its conductivity,
     * @param thickness its thickness,
     * @param posFlux positive terminal to plug a source (can be null)
     * @param negFlux negative terminal to plug a source (can be null)
     * @param dielectric the dielectric border,
     * @param epsilon dielectric's constant,
     * @param regCondDielec the region that would be at the intersection of the
     * dielectrric and the conductor, it must be on the conductor.<BR>
     * #WARNING: An inside (for the dielectric) node will be searched on the
     * conductor region. So it may fail if the conductor exceed from the
     * dielectric.
     */
    public TestTwo_Preconditioner(VolumeRegion conductor, double sigma, double thickness, SurfaceRegion posFlux, SurfaceRegion negFlux,
            SurfaceRegion dielectric, double epsilon,
            SurfaceRegion regCondDielec) {
        this(conductor, sigma, thickness, posFlux, negFlux, dielectric, epsilon, regCondDielec, null, false);
    }

    /**
     * Main constructor 101.
     *
     * @param conductor conductor region,
     * @param sigma its conductivity,
     * @param thickness its thickness,
     * @param posFlux positive terminal to plug a source (can be null)
     * @param negFlux negative terminal to plug a source (can be null)
     * @param dielectric the dielectric border,
     * @param epsilon dielectric's constant,
     * @param regCondDielec the region that would be at the intersection of the
     * dielectrric and the conductor, it must be on the conductor;
     * @param refNode node used to set the normals,
     * @param inside true if the refNode is 'inside' the dielectric, false if
     * its exterior.
     */
    protected TestTwo_Preconditioner(Region conductor, double sigma, double thickness, Region posFlux, Region negFlux,
            SurfaceRegion dielectric, double epsilon,
            SurfaceRegion regCondDielec, Node refNode, boolean inside) {

        borne = posFlux != null && negFlux != null;

        if (conductor instanceof VolumeRegion) {
            isSurf = false;
            this.ep = -1.0;
            // Generate the border and take off the bornes if there are any
            SurfaceRegion border = (SurfaceRegion) conductor.generateBorder();
            if (borne) {
                border = new SurfaceRegion((ElementSurfSetHomogene) border.getComplement((SurfaceRegion) posFlux, 2));
                border = new SurfaceRegion((ElementSurfSetHomogene) border.getComplement((SurfaceRegion) negFlux, 2));
                border = new SurfaceRegion((ElementSurfSetHomogene) border.getComplement((SurfaceRegion) regCondDielec, 2));
                // Renumber the conductor with the [ interior, exterior, b1, b2 ] !
                FDcond = new FaceDeg1((ElementSetHomogene) conductor.getElementSet(),
                        new FaceConstraint[]{new ExternalDriveFaceConstraint(border),
                            new ExternalDriveFaceConstraint(regCondDielec),
                            new ExternalDriveFaceConstraint(posFlux), new ExternalDriveFaceConstraint(negFlux)});
                indBornesCons = 2;
                // The exterior domain, have to add again the bornes
//                border = (SurfaceRegion) border.generateUnionWith(posFlux);
//                border = (SurfaceRegion) border.generateUnionWith(negFlux);
                border = (SurfaceRegion) border.generateUnionWith(regCondDielec);
                CondExt = new Cell((ElementSetHomogene) border.getElementSet(), new ExternalDriveCellConstraint(regCondDielec));
            } else {
                indBornesCons = -1;
                FDcond = new FaceDeg1((ElementSetHomogene) conductor.getElementSet(), new ExternalDriveFaceConstraint(border));
                // The exterior domain
                CondExt = new Cell((ElementSetHomogene) border.getElementSet(), new ExternalDriveCellConstraint(regCondDielec));
            }
        } else if (conductor instanceof SurfaceRegion) {
            isSurf = true;
            this.ep = thickness;
            // Generate the border and take off the bornes if there are any
            LineRegion borderNull = (LineRegion) conductor.generateBorder();
            if (borne) {
                borderNull = new LineRegion((ElementLinSetHomogene) borderNull.getComplement((LineRegion) posFlux, 1));
                borderNull = new LineRegion((ElementLinSetHomogene) borderNull.getComplement((LineRegion) negFlux, 1));
                // When the conductor is a surface, we have to descativavte the border.
                FDcond = new FaceDeg1((ElementSetHomogene) conductor.getElementSet(),
                        new FaceConstraint[]{new FaceRealDirichlet(borderNull, 0.0),
                            new ExternalDriveFaceConstraint(posFlux), new ExternalDriveFaceConstraint(negFlux)});
                indBornesCons = 0;
            } else {
                indBornesCons = -1;
                FDcond = new FaceDeg1((ElementSetHomogene) conductor.getElementSet(), new FaceRealDirichlet(borderNull, 0.0));
            }
            // The exterior is the whole conductor
            CondExt = new Cell((ElementSetHomogene) conductor.getElementSet(), new ExternalDriveCellConstraint(regCondDielec));
        } else {
            throw new IllegalArgumentException("Conductor must be either a Surface or a Volume region !");
        }
        ElementSurfSetHomogene ESHdielec = (ElementSurfSetHomogene) dielectric.getElementSet();
        //
        // Reset the orientation of the dielectric mesh, taking as reference a
        // node on the conductor supposed not to be on the dielectrics and 'inside'.
        //
        if (refNode != null) {
            ESHdielec.surfOrientation(refNode.getCoord(), inside);
        } else {
            //*
            ESHdielec.surfOrientation(new double[]{1e6, 1e6, 1e6}, true);
            /*/
             ESHdielec.surfOrientation(new double[]{0.001, 0.001, 0.0}, false);
             //*/
        }

        CellDielec = new Cell(ESHdielec);

        ESHdielec.plotElementSet("d:/tmp/dielec.vtk");
        FDcond.getElementSet().plotElementSet("d:/tmp/cond.vtk");
        CondExt.getElementSet().plotElementSet("d:/tmp/condExt.vtk");

        System.out.println("===========================");
        System.out.println("== Function Spaces Info. ==");
        System.out.println("===========================");
        System.out.println("_Conductor:   " + FDcond.getNbElement() + " elements");
        System.out.println("              " + FDcond.getActiveDofCount() + " faces");
        System.out.println("              " + Arrays.toString(FDcond.getImplicitConstraintAdresses()));
        System.out.println("              " + Arrays.toString(FDcond.getImplicitConstraintCount()));
        System.out.println("");
        System.out.println("_Dielectric:  " + CellDielec.getNbElement() + " elements/faces");
        System.out.println("");
        System.out.println("_Ext. Border: " + CondExt.getNbElement() + " elements/faces");
        System.out.println("              " + Arrays.toString(CondExt.getImplicitConstraintAdresses()));
        System.out.println("              " + Arrays.toString(CondExt.getImplicitConstraintCount()));
        System.out.println("");

        this.sigma = sigma;
        this.epsilon = epsilon;

        this.solveurCircuit = new SolveurMaillesIndependantes(2);
        this.solveurCircuit.setBloc(this, 0, 0);
    }

    @Override
    public int[][] getTopologie() {
        ExportTopologieGMSH topo = null;
        if (this.exportTopo != null) {
            topo = new ExportTopologieGMSH(exportTopo, FDcond);
        }
        // Construction de la matrice topologique
        int matTopo[][] = new int[getNbLignes()][4];

        // Indice de noeud libre
        int indNextNode = this.FDcond.getElementSet().getNbElement();
        // Count the added nodes on the border
        int cpt = 0;
        // Recuperation des liens entre face et elements
        int lienFaceElement[][] = this.FDcond.getFacesSet().getPosElementFace();

        for (int i = 0; i < this.FDcond.getActiveDofCount(); i++) {
            // Check if we are on the border.
            if (lienFaceElement[i][0] == -1 || lienFaceElement[i][1] == -1) {

                if (borne && FDcond.getImplicitConstraintAdresses()[indBornesCons] <= i && i < FDcond.getImplicitConstraintAdresses()[indBornesCons + 1]) {
                    matTopo[i][0] = i; // Numero de branche = numero de face
                    matTopo[i][1] = lienFaceElement[i][0] == -1 ? lienFaceElement[i][1] : lienFaceElement[i][0]; // Numero du noeud interieur
                    matTopo[i][2] = indNextNode; // Numero du noeud sur la face
                    matTopo[i][3] = -2; // Branche de type inductive
                    System.out.println("Face de bord 'flux pos.' d'indice " + i + " : N1 = " + matTopo[i][1] + " N2 = " + matTopo[i][2]);
                } else if (borne && FDcond.getImplicitConstraintAdresses()[indBornesCons + 1] <= i) {
                    matTopo[i][0] = i; // Numero de branche = numero de face
                    matTopo[i][1] = lienFaceElement[i][0] == -1 ? lienFaceElement[i][1] : lienFaceElement[i][0]; // Numero du noeud interieur
                    matTopo[i][2] = indNextNode + 1; // Numero du noeud sur la face
                    matTopo[i][3] = -2; // Branche de type inductive
                    System.out.println("Face de bord 'flux neg.' d'indice " + i + " : N1 = " + matTopo[i][1] + " N2 = " + matTopo[i][2]);
                } else {
                    // Should happens on the volumic case !
                    matTopo[i][0] = i; // Numero de branche = numero de face
                    matTopo[i][1] = lienFaceElement[i][0] == -1 ? lienFaceElement[i][1] : lienFaceElement[i][0]; // Numero du noeud interieur
                    matTopo[i][2] = indNextNode + 2 + cpt; // Numero du noeud sur la face
                    matTopo[i][3] = -2; // Branche de type inductive   
                    System.out.println("Face de bord d'indice " + i + " : N1 = " + matTopo[i][1] + " N2 = " + matTopo[i][2]);
                    cpt++;
                }
                if (topo != null) {
                    topo.addBorderBranch(matTopo[i][0], matTopo[i][1], matTopo[i][2]);
                }
            } else {
                // On est sur un dof interne
                matTopo[i][0] = i; // Numero de branche = numero de face
                matTopo[i][1] = lienFaceElement[i][0] < lienFaceElement[i][1] ? lienFaceElement[i][0] : lienFaceElement[i][1]; // Numero du noeud interieur de plus petite valeure (normal allant du plus petit numero d'element au plus grand)
                matTopo[i][2] = lienFaceElement[i][0] > lienFaceElement[i][1] ? lienFaceElement[i][0] : lienFaceElement[i][1]; // Numero du noeud interieur de plus grande valeure (normal allant du plus petit numero d'element au plus grand)
                matTopo[i][3] = -2; // Branche de type inductive
                System.out.println("Face interieure d'indice " + i + " : N1 = " + matTopo[i][1] + " N2 = " + matTopo[i][2]);

                if (topo != null) {
                    topo.addInnerBranch(matTopo[i][0], matTopo[i][1], matTopo[i][2]);
                }
            }
        }
        // Voir une meilleur solution que le '+2' !
        indNextNode += 2 + cpt;
        //
        // CAPACITIVE PART.
        // Connect all the node to the infinite node
        //
        if (computePc) {
            if(topo!=null){
                topo.setNodeInfinite(indNextNode, new double[3]);
            }
            int offSet = FDcond.getActiveDofCount();
            if (isSurf) {
                // 
                // In the volumic case, only have to pull some capacitive branches
                // from all the elements
                //
                ElementSetHomogene elementSet = FDcond.getElementSet();
                for (int i = 0; i < elementSet.getNbElement(); i++) {
                    int numNode = elementSet.getElements(i).getNum();
                    matTopo[i + offSet][0] = offSet + i;
                    matTopo[i + offSet][1] = numNode;
                    matTopo[i + offSet][2] = indNextNode; // Numero du noeud commun      
                    matTopo[i + offSet][3] = -3; // Branche de type capacitive                
                    System.out.println("Face exterieur capacitive d'indice " + (offSet + i) + " : N1 = " + matTopo[i + offSet][1] + " N2 = " + matTopo[i + offSet][2]);

                    if (topo != null) {
                        topo.addCapacitiveBranch(matTopo[i + offSet][0], matTopo[i + offSet][1]);
                    }
                }
            } else {
                // The next node is the capacitive node
                int noeudCapaCommun = indNextNode;
                // Reste the counters for the next loop.
                indNextNode = FDcond.getNbElement();
                cpt = 0;
                for (int i = 0; i < this.FDcond.getActiveDofCount(); i++) {
                    // Check if we are on the border.
                    if (lienFaceElement[i][0] == -1 || lienFaceElement[i][1] == -1) {
                        if (borne && FDcond.getImplicitConstraintAdresses()[indBornesCons] <= i && i < FDcond.getImplicitConstraintAdresses()[indBornesCons + 1]) {
//                        matTopo[offSet][0] = offSet; // Numero de branche = numero de face
//                        matTopo[offSet][1] = indNextNode;
//                        matTopo[offSet][2] = noeudCapaCommun;
//                        matTopo[offSet][3] = -3; // Branche de type inductive
//                        System.out.println("Face exterieure de bord 'flux pos.' d'indice " + offSet + " : N1 = " + matTopo[offSet][1] + " N2 = " + matTopo[offSet][2]);
                        } else if (borne && FDcond.getImplicitConstraintAdresses()[indBornesCons + 1] <= i) {
//                        matTopo[offSet][0] = offSet; // Numero de branche = numero de face
//                        matTopo[offSet][1] = indNextNode + 1;
//                        matTopo[offSet][2] = noeudCapaCommun; // Numero du noeud sur la face
//                        matTopo[offSet][3] = -3; // Branche de type inductive
//                        System.out.println("Face exterieure de bord 'flux neg.' d'indice " + offSet + " : N1 = " + matTopo[offSet][1] + " N2 = " + matTopo[offSet][2]);
                        } else {
                            // Should happens on the volumic case !
                            matTopo[offSet][0] = offSet; // Numero de branche = numero de face
                            matTopo[offSet][1] = indNextNode + 2 + cpt;
                            matTopo[offSet][2] = noeudCapaCommun; // Numero du noeud sur la face
                            matTopo[offSet][3] = -3; // Branche de type inductive                    
                            cpt++;
                            System.out.println("Face exterieur capacitive d'indice " + offSet + " : N1 = " + matTopo[offSet][1] + " N2 = " + matTopo[offSet][2]);
                            offSet++;
                        }
                    }
                }
            }
        }
        if (topo != null) {
            topo.close();
        }
        if (equiPot != null) {
            int iBrRef = equiPot[0];
            // Maximum number of repetition
            int nbRep = FDcond.getElementSet().getElements(0).getGeom().getNbCote();
            // Loop on all the indices to put at the same potential
            for (int i = 1; i < equiPot.length; i++) {
                // Look on the topologie
                cpt = 0;// The counter will be uneffective only for the border nodes
                // While we haven't found all the shared branches and until the end
                for (int j = 0; j < matTopo.length && cpt < nbRep; j++) {
                    if (matTopo[j][1] == equiPot[i]) {
                        matTopo[j][1] = iBrRef;
                        cpt++;
                    } else if (matTopo[j][2] == equiPot[i]) {
                        matTopo[j][2] = iBrRef;
                        cpt++;
                    }
                }
            }
        }

        return matTopo;
    }

    public void setExportTopologie(String path) {
        this.exportTopo = path;
    }

    /**
     * Integration of the matrices [R], [L] and [Pc]
     */
    protected void integrationConductor() {
        //
        // R
        //
        System.out.println("===== Integration: [R] =====");
        long t = System.currentTimeMillis();
        RowVector v = new RowVector(FDcond.getNbElement());
        if (isSurf) {
            v.setAllElements(1 / (sigma * ep));
        } else {
            v.setAllElements(1 / sigma);
        }
        // Create the quanty 1/sigma 
        Cell support = new Cell(FDcond.getElementSet());
        RealScalarCellQuantity q = new RealScalarCellQuantity(support, v);
        // Assemblage de [R]
        FiniteElementFormulation EF = new FiniteElementFormulation(FDcond);
        EF.assembly(nbPtGCiblesInduc, q);
        R = ((StorageSparse) EF.getStore()).getSparseMat();
        System.out.println("Assembly Time: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");
        //
        // L 
        //
        System.out.println("===== Integration: [L] =====");
        IntegrationCorrectionStrategySource strategy;
        // Choose the integration strategy.
        System.out.println("Strategie d'integration: ");
        if (this.nbPtGSourcesInducCorr != -1) {
            System.out.println("Integrations avec une correction 'shift' de la diagonale: (" + nbPtGCiblesInduc + " , " + nbPtGSourcesInduc + ") -> (" + nbPtGCiblesInduc + " , " + nbPtGSourcesInducCorr + ") !");
            strategy = new SelfElementFixedGauss(nbPtGSourcesInduc, new InCreasedPGSourceNumber(nbPtGSourcesInducCorr));
        } else {
            System.out.println("Integrations sans correction ! Un 'shift' sur toute la matrice: (" + nbPtGCiblesInduc + " , " + nbPtGSourcesInduc + ") !");
            strategy = new NoCorrection(new FixedGaussSource(nbPtGSourcesInduc));
        }
        // The quantity is only mu_0
        v.setAllElements(mu0);
        // Create the quanty
        q = new RealScalarCellQuantity(support, v);
        // Assembly the matrix [L] as an H-matrix or Full-matrix
        GalerkinIntegralFormulation IV;
        t = System.currentTimeMillis();
        if (this.HmatCompression) {
            IV = new GalerkinIntegralFormulationHCA(this.FDcond, this.FDcond, new MultGvect(), strategy, this.nbPtGCiblesInduc, this.nbPtGCiblesInduc,
                    epsHmat, kmax, nmin, order, eta, recomp);
            IV.assembly(q);
            this.Lh = (StorageHmatrix) IV.getStore();
        } else {
            IV = new GalerkinIntegralFormulationFull(this.FDcond, this.FDcond, new MultGvect(), strategy, this.nbPtGCiblesInduc);
            IV.assembly(q);
            this.L = ((StorageFull) IV.getStore()).getMatrix();
        }
        System.out.println("Assembly Time: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");
        //
        // [Pc]
        //
        if (computePc) {
            System.out.println("===== Integration: [P_c] =====");
            t = System.currentTimeMillis();
            if (fullCapa) {
                System.out.println("Integrations fully analytical !");
                strategy = new FullCorrection(new AnalyticalCorrection());
            } else {
                System.out.println("Numerical integrations with an analytical correction on the diagonal !");
                strategy = new SelfElementFixedGauss(nbPtGSourcesCapa, new AnalyticalCorrection());
            }
            if (HmatCompression) {
                // For the Hmatrix assembly, we have to scale the rows and the columns by the surfaces.
                // We also have to scale the columns by 1 / eps, eps = eps0 on the interface with air.
                // 
                // 1 / (eps * Sj)
                double tab[] = new double[CondExt.getActiveDofCount()];
                for (int i = 0; i < CondExt.getImplicitConstraintAdresses()[0]; i++) {
                    tab[i] = CondInvSurf[i] / eps0;
                }
                // Interface conductor / dielectric.
                for (int i = CondExt.getImplicitConstraintAdresses()[0]; i < CondExt.getActiveDofCount(); i++) {
                    tab[i] = CondInvSurf[i] / epsilon;
                }
                v = new RowVector(tab);
                // Scale the columns
                q = new RealScalarCellQuantity(support, v);
                IV = new GalerkinIntegralFormulationHCA(CondExt, CondExt, new MultG(), strategy, nbPtGCiblesCapa, nbPtGSourcesCapa,
                        epsHmat, kmax, nmin, order, eta, recomp);
                IV.assembly(q);
                Pc = IV.getStore();
                // Finally sclae the rows by 1 / Si
                ((StorageHmatrix) Pc).scaleRows(CondInvSurf);
            } else {
                IV = new GalerkinIntegralFormulationFull(CondExt, CondExt, new MultG(), strategy, nbPtGCiblesCapa);
                IV.assembly();
                Matrix M = ((StorageFull) IV.getStore()).getMatrix();
                // Scale by 1/(eps0 * Si  * Sj) the part conductor / air
                for (int j = 0; j < CondExt.getImplicitConstraintAdresses()[0]; j++) {
                    for (int i = 0; i < CondExt.getActiveDofCount(); i++) {
                        M.setElement(i, j, CondInvSurf[j] * CondInvSurf[i] * M.getElement(i, j) / eps0);
                    }
                }
                // Scale by 1/(epsilon * Si  * Sj) the part conductor / dielectric
                double cst;
                if (isSurf) {
                    cst = (1 / epsilon);//(1 / eps0) - (1 / epsilon);
                } else {
                    cst = 1 / epsilon;
                }
                for (int j = CondExt.getImplicitConstraintAdresses()[0]; j < CondExt.getActiveDofCount(); j++) {
                    for (int i = 0; i < CondExt.getActiveDofCount(); i++) {
                        M.setElement(i, j, CondInvSurf[j] * CondInvSurf[i] * M.getElement(i, j) * cst);
                    }
                }
                // Re-assign as a Storage
                Pc = new StorageFull(M);
            }
            System.out.println("Assembly Time: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");
        }
    }

    public void setCorrectionVoisin(boolean bull) {
        voisCorr = bull;
    }

    protected void integrationDielectric() {
        // We have to create a new Cell, projected on the normals.
        NormalCell projNdielec = new NormalCell(CellDielec);
        ExportGmshCell exp = new ExportGmshCell(projNdielec, "D:/nProjDielec.msh");
        //
        // Choose the integration stategy.
        //
        IntegrationCorrectionStrategySource strategy;
        if (fullCapa) {
            System.out.println("Integrations fully analytical !");
            strategy = new FullCorrection(new AnalyticalCorrection());
        } else {
            if (voisCorr) {
                System.out.println("Numerical integrations with an analytical correction on the neighbors !");
                strategy = new NearElementSphereFixedGauss(nbPtGSourcesCapa, new AnalyticalCorrection(), 3.0);
            } else {
                System.out.println("Numerical integrations with an analytical correction on the diagonal !");
                strategy = new SelfElementFixedGauss(nbPtGSourcesCapa, new AnalyticalCorrection());
            }
        }
        //
        //
        // [B]
        //
        //
        double sign = -1;
        System.out.println("===== Integration: [B] =====");
        long t = System.currentTimeMillis();
        RowVector v = new RowVector(CondExt.getActiveDofCount());
        // Compute 1/ eps0 for the interface air / conductor
        for (int i = 0; i < CondExt.getImplicitConstraintAdresses()[0]; i++) {
            v.setElement(i, sign * CondInvSurf[i] / eps0);
        }
        // Compute 1/ epsilon for the interface dielectric / conductor
        double cstInter;
        if (isSurf) {
            cstInter = sign * (1 / epsilon);
        } else {
            cstInter = sign / epsilon;
        }
        for (int i = CondExt.getImplicitConstraintAdresses()[0]; i < CondExt.getActiveDofCount(); i++) {
            v.setElement(i, CondInvSurf[i] * cstInter);//* ((1 / epsilon) + ((1 / eps0))));
        }
        // Create the quanty 1/epsilon
        RealScalarCellQuantity q = new RealScalarCellQuantity(CondExt, v);

        GalerkinIntegralFormulation IF = new GalerkinIntegralFormulationFull(projNdielec, CondExt, new NegMultDG(), strategy, nbPtGCiblesCapa);
        IF.assembly(q);
        //
        // A voir s'il ne faut pas multiplier par les surfaces !
        /*
         Matrix Bm = ((StorageFull) IF.getStore()).getMatrix();
         Bm.scaleRow(DielecInvSurf);
         B = new StorageFull(Bm);
         /*/
        B = IF.getStore();
        //*/
        System.out.println("Assembly Time: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");
        //
        //
        // [A]
        //
        // a verifier lexpression !
        System.out.println("===== Integration: [A] =====");
        t = System.currentTimeMillis();
        v = new RowVector(CellDielec.getActiveDofCount());
        for (int i = 0; i < v.getColumnCount(); i++) {
            v.setElement(i, sign * DielecInvSurf[i] / eps0);
        }
        q = new RealScalarCellQuantity(CellDielec, v);
        IF = new GalerkinIntegralFormulationFull(projNdielec, CellDielec, new NegMultDG(),
                /*
                 new SelfElementAdaptiveGauss(4, 16, 1e-5, new InCreasedPGSourceNumber(25)),
                 /*/
                strategy,
                //*/
                nbPtGCiblesCapa);
        IF.assembly(q);
        //
        // We have to add something on the diagonale !
        //
        Matrix Am = ((StorageFull) IF.getStore()).getMatrix();
//        Am.scaleRow(DielecInvSurf);
        double cst = (0.5 / eps0) * (eps0 + epsilon) / (epsilon - eps0);
        for (int i = 0; i < CellDielec.getActiveDofCount(); i++) {
//            Am.setElement(i, i, Am.getElement(i, i) - cst);
            Am.setElement(i, i, -cst);
        }
        //
        // A voir s'il ne faut pas multiplier par les surfaces !
        //
        A = new StorageFull(Am);
        System.out.println("Assembly Time: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");
        //
        //
        // [P_d]
        //
        //
        System.out.println("===== Integration: [P_D] =====");
        t = System.currentTimeMillis();
        IF = new GalerkinIntegralFormulationFull(CondExt, CellDielec, new MultG(), strategy, nbPtGCiblesCapa);
        IF.assembly();
        Matrix PD = ((StorageFull) IF.getStore()).getMatrix();
        double tmp;
        for (int i = 0; i < PD.getRowCount(); i++) {
            tmp = CondInvSurf[i] / eps0;// Fausse bonne idée ?
            for (int j = 0; j < PD.getColumnCount(); j++) {
                PD.setElement(i, j, PD.getElement(i, j) * tmp * DielecInvSurf[j]);
            }
        }
        Pd = new StorageFull(PD);
        System.out.println("Assembly Time: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");
    }

    /**
     * Compute M_ij = s * M_ij / r[i] / c[j]
     *
     * @param M
     * @param s
     * @param r
     * @param c
     */
    protected void scaleRowCol(Matrix M, double s, double r[], double c[]) {
        for (int i = 0; i < M.getRowCount(); i++) {
            for (int j = 0; j < M.getColumnCount(); j++) {
                M.setElement(i, j, s * M.getElement(i, j) / r[i] / c[j]);
            }
        }
    }

    /**
     * Compute Jd using Q. Jd_j = j * w * [(eps1 * eps2) / (eps1 - eps2)] * Q_j
     *
     * @return Jd
     */
    public double[] getJd(double f) {
        double jd[] = new double[Q.length];
        double omega = 2 * Math.PI * f;
        // Pour le moment on considère un seul dielectrique.
        double cst = (epsilon * eps0) / (epsilon - eps0);
        for (int i = 0; i < CellDielec.getActiveDofCount(); i++) {
            jd[2 * i] = -omega * cst * Q[2 * i + 1];
            jd[2 * i + 1] = omega * cst * Q[2 * i];
        }
        return jd;
    }

    /**
     * Compute the integration of all the matrices.
     */
    public void integration() {
        computeSurfaces();
        System.out.println("_CONDUCTOR ASSEMBLIES");
        long t = System.currentTimeMillis();
        this.integrationConductor();
        System.out.println("\n_Assembly Time of the conductor: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");

        System.out.println("_DIELECTRIC ASSEMBLIES");
        t = System.currentTimeMillis();
        this.integrationDielectric();
        System.out.println("\n_Assembly Time of the dielectric: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");
    }

    /**
     * Check if the matrices contains any Infinity or NaN. The Finite Element
     * Matrix [R] is not checked.
     */
    public void checkMatricesSingularities() {
        if (HmatCompression) {
            System.out.println("a implementer !");
        } else {
            if (checkSingMat(L)) {
                System.err.println("Singularity found in [L] !");
            } else {
                System.out.println("No singularity in [L]");
            }
            if (checkSingMat(((StorageFull) Pc).getMatrix())) {
                System.err.println("Singularity found in [Pc] !");
            } else {
                System.out.println("No singularity in [Pc]");
            }
            if (checkSingMat(((StorageFull) B).getMatrix())) {
                System.err.println("Singularity found in [B] !");
            } else {
                System.out.println("No singularity in [B]");
            }
            if (checkSingMat(((StorageFull) A).getMatrix())) {
                System.err.println("Singularity found in [A] !");
            } else {
                System.out.println("No singularity in [A]");
            }
            if (checkSingMat(((StorageFull) Pd).getMatrix())) {
                System.err.println("Singularity found in [Pd] !");
            } else {
                System.out.println("No singularity in [Pd]");
            }
        }
    }

    protected boolean checkSingMat(Matrix M) {
        for (int i = 0; i < M.getRowCount(); i++) {
            for (int j = 0; j < M.getColumnCount(); j++) {
                if (Double.isInfinite(M.getElement(i, j))) {
                    System.out.println("Infinite found ! At (" + i + " , " + j + ") !");
                    return true;
                }
                if (Double.isNaN(M.getElement(i, j))) {
                    System.out.println("NaN found ! At (" + i + " , " + j + ") !");
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Compute the surfaces of the exterior faces.
     */
    protected void computeSurfaces() {
        // Table des surfaces exterieures conductrices
        CondInvSurf = new double[CondExt.getActiveDofCount()];
        for (int i = 0; i < CondExt.getActiveDofCount(); i++) {
            CondInvSurf[i] = 1 / ((GeomElementSurf) CondExt.getElementSet().getElements(i).getGeom()).getSurface();
        }
        // Table des surfaces dielectriques
        DielecInvSurf = new double[CellDielec.getActiveDofCount()];
        for (int i = 0; i < CellDielec.getActiveDofCount(); i++) {
            DielecInvSurf[i] = 1 / ((GeomElementSurf) CellDielec.getElementSet().getElements(i).getGeom()).getSurface();
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
     * @param nbPtsGaussSourcesCorrection Nombre de points de Gauss pour la
     * correction numérique de l'élément sur lui-même
     */
    public void setPtsDeGaussInductifs(int nbPtsGaussSourcesCorrection, int nbPtsGaussCibles, int nbPtsGaussSources) {
        this.nbPtGCiblesInduc = nbPtsGaussCibles;
        this.nbPtGSourcesInduc = nbPtsGaussSources;
        this.nbPtGSourcesInducCorr = nbPtsGaussSourcesCorrection;
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
     * Affectation du bloc circuit electrique pure (description des elements
     * autours de la partie PEEC : R, L , C , SU, SI
     *
     * @param circuitPur Description des elements de circuit exterieur a la
     * geometrie
     */
    public void setCircuitElectrique(BlocCircuit circuitPur) {
        System.out.println("===============================");
        System.out.println("     SET CIRCUIT ELECTRIC      ");
        System.out.println("===============================");
        System.out.println(circuitPur.toString());
        this.solveurCircuit.setBloc(circuitPur, 1, 1);
    }

    /**
     * Assembly the matrix blocs of Zm
     *
     * @param f current frequency
     */
    private void assembleZm(double f) {
        fcurr = f;
        System.out.println("f= "+f);
        System.out.println("ep= "+ep);
        double omega = 2 * Math.PI * f;
        // [  R + j w L   |        0      ]
        // [--------------+---------------]
        // [      0       |   1/(j w) Pc  ] 
        if (fcurr != 0 && computePc) {
            Zm = new MatrixBlocComplex(2, 2);
        } else {
            Zm = new MatrixBlocComplex(1, 1);
        }
        boolean isHF;
        double delta = 0;
        if (isSurf) {
            delta = Math.sqrt((1 / (sigma * Math.PI * f * mu0))); // epaisseur de peau
            isHF = delta < this.ep;
        } else {
            isHF = false;
        }
        if (isHF) {
            System.err.println("High Frequency ! Scale by 1/G ! ep= " + ep + " > delta= " + delta);
            double G[] = this.computeInverseG(f, ep, delta);
            // 1/G .* R + j omega L
            Zm.setBloc(0, 0, new ComplexMatrixBlocAdd(
                    new ComplexMatrixBlocScale(new ComplexMatrixBloc(R, null), G[0], G[1]),
                    new ComplexMatrixBloc(null, new RealMatrixBlocScale(omega, HmatCompression ? Lh : new StorageFull(L))))
            );
        } else {
            if (delta != 0) {
                System.err.println("'Low' Frequency ! ep= " + ep + " < delta= " + delta);
            }
            // R + j omega L
            Zm.setBloc(0, 0, new ComplexMatrixBloc(R, new RealMatrixBlocScale(omega, HmatCompression ? Lh : new StorageFull(L))));
        }
        if (computePc && fcurr != 0) {
            // 1 / (j  omega) Pc
            Zm.setBloc(1, 1, new ComplexMatrixBloc(null, new RealMatrixBlocScale(-1 / omega, Pc)));
        }
    }

    /**
     * Assemble the whole matrix.<BR>
     *
     * [ Zm , M*P_d ] <BR>
     * [ B*M^t , A ] <BR>
     *
     * @param f the current frequency.
     */
    private void assembleMatrixSolution(double f) {
        if (f != fcurr || Zfull == null) {
            Z = new g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D(getZ(f));
        } else {
            Z = new g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.Basic2D(Zfull);
        }
    }

    public double[][] resolutionDirectePureConductor(double f) {
        // Compute the matrices for the 1st time
        if (fcurr == -1) {
            computeSurfaces();
            System.out.println("_CONDUCTOR ASSEMBLIES");
            long t = System.currentTimeMillis();
            this.integrationConductor();
            System.out.println("\n_Assembly Time of the conductor: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");
        }
        // Only change to the current frequency
        if (fcurr != f) {
            assembleZm(f);
        }
        // Check the topologie before computing the matrices !
        if (solveurCircuit.getCircuit() == null) {
            solveurCircuit.makeCircuit();
        }
        return solveurCircuit.resolutionDirecte(new double[]{f}, true);

    }
    
    public double[][] resolutionIterativePureConductor(double f) {
        // Compute the matrices for the 1st time
        if (fcurr == -1) {
            computeSurfaces();
            System.out.println("_CONDUCTOR ASSEMBLIES");
            long t = System.currentTimeMillis();
            this.integrationConductor();
            System.out.println("\n_Assembly Time of the conductor: " + (System.currentTimeMillis() - t) * 1e-3 + " sec");
        }
        // Only change to the current frequency
        if (fcurr != f) {
            assembleZm(f);
        }
        // Check the topologie before computing the matrices !
        if (solveurCircuit.getCircuit() == null) {
            solveurCircuit.makeCircuit();
        }
        return solveurCircuit.resolutionIterative(new double[]{f}, true);

    }

    /**
     * Direct Solution using LU-factorization.
     *
     * @param f solution frequency
     * @return datas solution on the conductor. Use getQ() to get the results on
     * the dielectric.
     */
    public double[][] resolutionDirecte(double f) {
        // Check the topologie before computing the matrices !
        if (solveurCircuit.getCircuit() == null) {
            solveurCircuit.makeCircuit();
        }
        // Compute the matrices for the 1st time
        if (fcurr == -1) {
            integration();
        }
        // Only change to the current frequency
        if (fcurr != f) {
            assembleZm(f);
            Zfull = getZ(f);
        }
        //////////////////////////////
        //  Compute the RHS-vector  //
        //////////////////////////////
        // Conductor part
        double bCond[] = solveurCircuit.constructionSecondMembre(f);
        // Dielectric part
        // Decommenter si bDielec != 0
//        double bDielec[] = new double[2 * CellDielec.getActiveDofCount()];
        // Assembly the two parts
        double b[] = new double[bCond.length + 2 * CellDielec.getActiveDofCount()];
        System.arraycopy(bCond, 0, b, 0, bCond.length);
//        System.arraycopy(bDielec, 0, b, bCond.length, bDielec.length);

        ///////////////////////////////////
        //  Compute the Matrix to Solve  //
        ///////////////////////////////////
        // Decomposition LU de la matrice z
        long t = System.currentTimeMillis();
        SolverLUBasic2D solveur = new SolverLUBasic2D(Zfull, 1e-15, false, true, false);
        // Resolution des systemes triangulaires LU
        double sol[] = solveur.solveLU(b);
        System.out.println("  Duree decomposition LU et resolution  : " + (System.currentTimeMillis() - t) + " ms");
        //
        // Extract the DOF
        //
        // Courants de mailles        
        MatriceIncidence mi = solveurCircuit.getMI();
        int nbMI = mi.getNbLignes();
        double im[] = new double[2 * nbMI];
        System.arraycopy(sol, 0, im, 0, im.length);
        // 'Charges' sur le dielectrique
        Q = new double[2 * CellDielec.getActiveDofCount()];
        System.arraycopy(sol, im.length, Q, 0, Q.length);

        return postSolution(im, f);
    }

    public double[][] resolutionIterative(double f) {
        // Check the topologie before computing the matrices !
        if (solveurCircuit.getCircuit() == null) {
            solveurCircuit.makeCircuit();
        }
        // Compute the matrices for the 1st time
        if (fcurr == -1) {
            integration();
        }
        // Only change to the current frequency
        if (fcurr != f) {
            assembleZm(f);
            assembleMatrixSolution(f);
        }
        //////////////////////////////
        //  Compute the RHS-vector  //
        //////////////////////////////
        // Conductor part
        double bCond[] = solveurCircuit.constructionSecondMembre(f);
        // Dielectric part
        // Decommenter si bDielec != 0
//        double bDielec[] = new double[2 * CellDielec.getActiveDofCount()];
        // Assembly the two parts
        double b[] = new double[bCond.length + 2 * CellDielec.getActiveDofCount()];
        System.arraycopy(bCond, 0, b, 0, bCond.length);
//        System.arraycopy(bDielec, 0, b, bCond.length, bDielec.length);

        System.out.println("b= " + Arrays.toString(b));
        ///////////////////////////////////
        //  Compute the Matrix to Solve  //
        ///////////////////////////////////
        long t = System.currentTimeMillis();
        FGMResComplex solver = new FGMResComplex(this.Z, null);
        solver.setInfoResolution(new double[]{1, 1e-8, 1, -1});
        double sol[] = solver.solve(new double[b.length], b);
        System.out.println("  Duree resolution iterative  : " + (System.currentTimeMillis() - t) + " ms");
        //
        // Extract the DOF
        //
        // Courants de mailles        
        MatriceIncidence mi = solveurCircuit.getMI();
        int nbMI = mi.getNbLignes();
        double im[] = new double[2 * nbMI];
        System.arraycopy(sol, 0, im, 0, im.length);
        // 'Charges' sur le dielectrique
        Q = new double[2 * CellDielec.getActiveDofCount()];
        System.arraycopy(sol, im.length, Q, 0, Q.length);

        return postSolution(im, f);
    }

    public double[] getQ() {
        return Q;
    }

    /**
     * Compute the matrix to solve.
     *
     * @param f
     * @return Z
     */
    public double[][] getZ(double f) {
        int nbLigne = solveurCircuit.getNbBranches();
        int indAvtApr[] = solveurCircuit.getIndiceAvantApres();
        // The matrix is :
        // [      Zm     ;  M*Pd  ] { Im }
        // [  B*Mt/(j w) ;    A   ] {  Q }
        //
        // Get M
        MatriceIncidence mi = solveurCircuit.getMI();
        int nbMI = mi.getNbLignes();
        double z[][] = new double[nbMI + CellDielec.getActiveDofCount()][2 * (nbMI + CellDielec.getActiveDofCount())];

        ////////////////////////////////////////////////////////////////////////
        // Zm
        ////////////////////////////////////////////////////////////////////////
        long t = System.currentTimeMillis();
        /*
         // tmp is now Complex
         double tmp[][] = new double[nbLigne][2 * nbLigne];
         // Get Zb
         tmp = getZbFull(tmp, 0, nbLigne - 1, 0, nbLigne - 1, f);
         // Zm = M Zb Mt
         System.out.println("mi= " + mi.getNbColonnes());
         System.out.println("nbLigne= " + nbLigne);
         tmp = mi.changementBasePleinComplexe(tmp);
         /*/
        double tmp[][] = solveurCircuit.getMatriceAResoudreFull(f);//getMatriceAResoudreFull(f); ou getZMIFull
        //*/
        z = setSubMatrix(tmp, z, 0, 0);
        System.out.println("Time to compute [Zm]: " + (System.currentTimeMillis() - t) + " msec.");
        ////////////////////////////////////////////////////////////////////////
        // M * Pd
        ////////////////////////////////////////////////////////////////////////
        t = System.currentTimeMillis();
        // We have to convert the byte[][] to double[][] !!!
        double av[][] = new double[mi.getNbLignes()][];
        byte[][] vals = mi.getListeValeurs();
        for (int i = 0; i < vals.length; i++) {
            av[i] = new double[vals[i].length];
            for (int j = 0; j < vals[i].length; j++) {
                av[i][j] = (double) vals[i][j];
            }
        }
        // M = mi as a SparseMatrix !
        SparseMatrixRowReal Mi = new SparseMatrixRowReal(nbLigne, mi.getListeIndexes(), av);
        // Bigger matrix

        Matrix M = new Matrix(nbLigne, CellDielec.getActiveDofCount());
        // We have to place the Matrix Pd in a bigger matrix (full of zeros...)
        M.setSubmatrix(FDcond.getActiveDofCount(), 0, ((StorageFull) Pd).getMatrix());
        // tmp is now real
        M.permuteRows(indAvtApr, new int[nbLigne]);
        tmp = productSparseFull(Mi, M);
        z = setRealSubMatrix(tmp, z, 0, nbMI);
        System.out.println("Time to compute [M * Pd]: " + (System.currentTimeMillis() - t) + " msec.");
        ////////////////////////////////////////////////////////////////////////
        // 1/(j w) B . Mt
        ////////////////////////////////////////////////////////////////////////
        t = System.currentTimeMillis();
        SparseMatrixRowReal Mit = Mi.getTranspose();
        // We have to place the Matrix Pd in a bigger matrix (full of zeros...)
        M = new Matrix(CellDielec.getActiveDofCount(), nbLigne);
        M.setSubmatrix(0, FDcond.getActiveDofCount(), ((StorageFull) B).getMatrix());
        M.permuteColumns(indAvtApr, new int[nbLigne]);
        tmp = productFullSparse(M, Mit);
        // Scale by -j / omega. 
        z = setImagSubMatrix(tmp, z, nbMI, 0, -1 / (2 * Math.PI * f));
        System.out.println("Time to compute 1/(j w) .* [B * Mt]: " + (System.currentTimeMillis() - t) + " msec.");
        ////////////////////////////////////////////////////////////////////////
        // [A]
        //////////////////////////////////////////////////////////////////////// 
        // Nothing to do.
        t = System.currentTimeMillis();
        z = setRealSubMatrix(((StorageFull) A).getMatrix(), z, nbMI, nbMI);
        System.out.println("Time to compute [A]: " + (System.currentTimeMillis() - t) + " msec.");
        return z;
    }

    /**
     * Sparse x Full
     *
     * @param F
     * @param M
     * @return F * M
     */
    protected double[][] productSparseFull(SparseMatrixRowReal F, Matrix M) {
        if (F.getColumns() != M.getRowCount()) {
            throw new InternalError("Non compatible ! (" + F.getRows() + " , " + F.getColumns() + ") x (" + M.getRowCount() + " x " + M.getColumnCount() + ")");
        }
        double res[][] = new double[F.getRows()][M.getColumnCount()];
        double coli[] = new double[M.getRowCount()];
        double[] tmp;
        for (int j = 0; j < M.getColumnCount(); j++) {
            // Column j of M
            M.column(j).get(coli);
            // F * M(:,j)
            tmp = F.product(coli, new double[F.getRows()]);
            // store the result column
            for (int i = 0; i < res.length; i++) {
                res[i][j] = tmp[i];
            }
        }
        return res;
    }

    /**
     * M x F
     *
     * @param M
     * @param F
     * @return M * F
     */
    protected double[][] productFullSparse(Matrix M, SparseMatrixRowReal F) {
        if (F.getRows() != M.getColumnCount()) {
            throw new InternalError("Non compatible !");
        }
        double res[][] = new double[M.getRowCount()][F.getColumns()];
        SparseMatrixRowReal Ft = F.getTranspose();
        int[][] ai = Ft.getAi();
        double[][] av = Ft.getAv();
        // Loop on the columns of F
        for (int col = 0; col < F.getColumns(); col++) {
            //
            // Compute a product Matrix * ColumnSparse
            // 
            double tmp[] = new double[M.getRowCount()];
            // Loop on the Full rows 
            for (int i = 0; i < M.getRowCount(); i++) {
                for (int j = 0; j < ai[col].length; j++) {
                    tmp[i] += av[col][j] * M.getElement(i, ai[col][j]);
                }
            }
            // Recopie dans la colonne col
            for (int i = 0; i < tmp.length; i++) {
                res[i][col] = tmp[i];
            }
        }
        return res;
    }

    /**
     *
     * @param z
     * @param ind
     * @return
     */
    protected double[][] changementDeBaseBranches(double[][] z, int ind[]) {
        double t[][] = new double[z.length][];
        for (int i = 0; i < t.length; i++) {
            t[i] = new double[z[i].length];
            for (int j = 0; j < z[i].length; j++) {
                t[ind[i]][ind[j]] = z[i][j];
            }
        }
        return t;
    }

    /**
     * RE-organize and consider the sources.
     *
     * @param sol
     * @param f
     * @return
     */
    protected double[][] postSolution(double sol[], double f) {

        MatriceIncidence mi = solveurCircuit.getMI();
        Circuit circuit = solveurCircuit.getCircuit();
        int nbSI = circuit.getNbSourcesI();
        int[] indiceApresAvant = solveurCircuit.getIndiceApresAvant();
        int[] indiceAvantApres = solveurCircuit.getIndiceAvantApres();
        int[] startInd = solveurCircuit.getStartInd();
        BlocCircuit[][] matriceBloc = solveurCircuit.getMatriceBloc();

        double ibvb[][] = new double[2][2 * circuit.getNbBranches()];
        // Recuperation des informations sur les sources
        // Construction du vecteur contenant les sources de courants branches et les sources de tensions branches
        double sib[] = new double[mi.getNbColonnes() * 2];
        double sub[] = new double[mi.getNbColonnes() * 2];
        // On parcours tous les blocs diagonaux et l'on recupere les valeurs des si
        for (int k = 0; k < matriceBloc.length; k++) {
            // Recuperation du vecteur local
            double tempSI[] = matriceBloc[k][k].getValeursSourceI(f);
            double tempSU[] = matriceBloc[k][k].getValeursSourceU(f);
            // Recopie dans le vecteur global avec reorganisation des vecteurs
            for (int j = 0; j < tempSI.length / 2; j++) {
                sib[2 * indiceAvantApres[j + startInd[k]]] = tempSI[2 * j];
                sib[2 * indiceAvantApres[j + startInd[k]] + 1] = tempSI[2 * j + 1];
                sub[2 * indiceAvantApres[j + startInd[k]]] = tempSU[2 * j];
                sub[2 * indiceAvantApres[j + startInd[k]] + 1] = tempSU[2 * j + 1];
            }
        }
        // Construction du vecteur simi (source de courant dans les mi) pour rappel on a une seule si par mi par construction (sinon cela signifie qu'il y a un probleme topologique du a un circuit mal construit)
        double simi[] = new double[mi.getNbLignes() * 2];
        // Recuperation de la table indexation de la matrice mi' afin de connaitre toutes les branches se trouvant sur les mi facilement
        int indexesMI[][] = mi.getListeIndexes();
        byte valeursMI[][] = mi.getListeValeurs();
        for (int k = indexesMI.length - nbSI; k < indexesMI.length; k++) {
            // On effectue une somme de facon bete meme s'il y a plein de 0 dans les vecteurs
            for (int j = 0; j < indexesMI[k].length; j++) {
                simi[2 * k] += valeursMI[k][j] * sib[indexesMI[k][j] * 2];
                simi[2 * k + 1] += valeursMI[k][j] * sib[indexesMI[k][j] * 2 + 1];
            }
        }
        // Recopie et modification du vecteur solution afin de retrouver les courants de mi veritable (sans les inconnues en tension)
        double imi[] = new double[sol.length];
        System.arraycopy(sol, 0, imi, 0, sol.length - nbSI * 2);
        System.arraycopy(simi, sol.length - nbSI * 2, imi, sol.length - nbSI * 2, nbSI * 2);
        // Calcul des courants branches
        double ib[] = mi.produitTransposeVecteurComplexe(imi);
        // Recopie des termes dans le vecteur final avec re-organisation
        for (int k = 0; k < ib.length / 2; k++) {
            ibvb[0][2 * indiceApresAvant[k]] = -ib[2 * k];
            ibvb[0][2 * indiceApresAvant[k] + 1] = -ib[2 * k + 1];
        }
        // Calcul des tensions branches
        double vb[] = solveurCircuit.produitZbIb(ib, f);
        for (int k = 0; k < vb.length / 2; k++) {
            ibvb[1][2 * indiceApresAvant[k]] = -vb[2 * k];
            ibvb[1][2 * indiceApresAvant[k] + 1] = -vb[2 * k + 1];
        }
//            ibvb[f.length+i] = this.produitZbIb(ibvb[i], f[i]);
        // Ajout des termes issues de la resolution correspondant aux tensions aux bornes des SI
        for (int j = sol.length / 2 - nbSI; j < sol.length / 2; j++) {
            // Recherche de l'indice de la source de courant se trouvant sur la mi en cours de traitement
            int indSIb = circuit.getIndiceSourceIMI(mi, j);
            int indIndexe = IntArrayList.binarySearch(indexesMI[j], 0, indexesMI[j].length, indSIb);
            ibvb[1][2 * indiceApresAvant[indexesMI[j][indIndexe]]] = sol[2 * j] * valeursMI[j][indIndexe];
            ibvb[1][2 * indiceApresAvant[indexesMI[j][indIndexe]] + 1] = sol[2 * j + 1] * valeursMI[j][indIndexe];
        }
        // Ajout des SU
        for (int j = 0; j < sub.length / 2; j++) {
            if (sub[2 * j] != 0 || sub[2 * j + 1] != 0) {
                ibvb[1][2 * indiceApresAvant[j]] = sub[2 * j];
                ibvb[1][2 * indiceApresAvant[j] + 1] = sub[2 * j + 1];
            }
        }
        return ibvb;
    }

    @Override
    public double[] produit(double[] ib, double[] vb, double f) {
        return this.Zm.product(ib, vb);
    }

    @Override
    public int getNbLignes() {
        if (!computePc) {
            return FDcond.getActiveDofCount();
        }

        return FDcond.getActiveDofCount() + CondExt.getActiveDofCount();
    }

    /**
     * Compute the constant 1/G to take in consideration the 'Courant de
     * Foucault'.
     *
     * @param f the frequency
     * @param ep the thickness
     * @param delta the skin depth
     * @return 1/G (complex) as an array
     */
    private double[] computeInverseG(double f, double ep, double delta) {
        double e2delta = ep / 2 / delta;
        double thRe = Math.sinh(2 * e2delta) / (Math.cosh(2 * e2delta) + Math.cos(2 * e2delta));
        double thIm = Math.sin(2 * e2delta) / (Math.cosh(2 * e2delta) + Math.cos(2 * e2delta));
        double[] one = {1, 1};
        double[] invG = ComplexOperator.div(one, thRe, thIm);
        invG[0] *= e2delta;
        invG[1] *= e2delta;
        return invG;
    }

    @Override
    public double[][] getZbFull(double[][] zb, int indDepLigne, int indFinLigne, int indDebCol, int indFinCol, double f) {
        double omega = 2 * Math.PI * f;
        System.out.println("f= "+f);
        boolean isHF;
        double G[] = null;
        if (isSurf) {
            double delta = Math.sqrt((1 / (sigma * Math.PI * f * mu0))); // epaisseur de peau
            isHF = delta < this.ep;
            if (isHF) {
                G = this.computeInverseG(f, ep, delta);
            }
        } else {
            isHF = false;
        }
        Matrix Pl = null;
        if (computePc) {
            Pl = ((StorageFull) Pc).getMatrix();
        }
        for (int i = indDepLigne; i <= indFinLigne; i++) {
            for (int j = indDebCol; j <= indFinCol; j++) {
                if (i < FDcond.getActiveDofCount() && j < FDcond.getActiveDofCount()) {
                    if (isHF) {
                        // 1/G .* R + j * omega .* L
                        double Rij = R.getQuick(i, j);
                        zb[i][2 * j] = G[0] * Rij;
                        zb[i][2 * j + 1] = G[1] * Rij + omega * L.getElement(i, j);
                    } else {
                        // R + j * omega * L
                        zb[i][2 * j] = R.getQuick(i, j);
                        zb[i][2 * j + 1] = omega * L.getElement(i, j);
                    }
                } else if (computePc && f != 0 && i >= FDcond.getActiveDofCount() && j >= FDcond.getActiveDofCount()) {
                    // 1/(j omega) Pc
                    zb[i][2 * j + 1] = -Pl.getElement(i - FDcond.getActiveDofCount(), j - FDcond.getActiveDofCount()) / omega;
                }
            }
        }
        return zb;
    }

    @Override
    public SparseVectorComplex getTermesExplicites(int indL, double f) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double[] product(double[] x, double[] res) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    /**
     * Put the matrix M in the Matrix Z, with the 1st index at (i, j).
     *
     * @param M matrix (complex) to copy in Z
     * @param Z matrix (complex) contains your mom
     * @param i row index
     * @param j column index
     * @return Z
     */
    private double[][] setSubMatrix(double[][] M, double[][] Z, int i, int j) {
        for (int k = 0; k < M.length; k++) {
            for (int l = 0; l < M[k].length / 2; l++) {
                Z[k + i][2 * (l + j)] = M[k][2 * l]; // Real Part
                Z[k + i][2 * (l + j) + 1] = M[k][2 * l + 1]; // Imaginar Part
            }
        }
        return Z;
    }

    /**
     * Put the matrix M in the real part of the Matrix Z, with the 1st index at
     * (i, j).
     *
     * @param M matrix (real) to copy in Z
     * @param Z matrix (complex) contains your mom
     * @param i row index
     * @param j column index
     * @return Z
     */
    private double[][] setRealSubMatrix(double[][] M, double[][] Z, int i, int j) {
        for (int k = 0; k < M.length; k++) {
            for (int l = 0; l < M[k].length; l++) {
                Z[k + i][2 * (l + j)] = M[k][l]; // Real Part
            }
        }
        return Z;
    }

    /**
     * Put the matrix M in the real part of the Matrix Z, with the 1st index at
     * (i, j).
     *
     * @param M matrix (real) to copy in Z
     * @param Z matrix (complex) contains your mom
     * @param i row index
     * @param j column index
     * @return Z
     */
    private double[][] setRealSubMatrix(Matrix M, double[][] Z, int i, int j) {
        for (int k = 0; k < M.getRowCount(); k++) {
            for (int l = 0; l < M.getColumnCount(); l++) {
                Z[k + i][2 * (l + j)] = M.getElement(k, l); // Real Part
            }
        }
        return Z;
    }

    /**
     * Put the matrix M in the imaginary part of the Matrix Z, with the 1st
     * index at (i, j).
     *
     * @param M matrix (real) to copy in Z
     * @param Z matrix (complex) contains your mom
     * @param i row index
     * @param j column index
     * @return Z
     */
    private double[][] setImagSubMatrix(double[][] M, double[][] Z, int i, int j, double scale) {
        for (int k = 0; k < M.length; k++) {
            for (int l = 0; l < M[k].length; l++) {
                Z[k + i][2 * (l + j) + 1] = scale * M[k][l]; // Complex Part
            }
        }
        return Z;
    }

    @Override
    public double[] getValeursSourceI(double f) {
        double t[] = new double[this.getNbLignes() * 2];
//        t[2 * 125] = 1.0;
//        t[2 * 130] = -1.0;
        return t;

    }

    @Override
    public double[] getValeursSourceU(double f) {
        return new double[2 * this.getNbLignes()];

    }

    /**
     *
     * @return MI solver.
     */
    public SolveurMaillesIndependantes getSolveurCircuit() {
        return solveurCircuit;
    }

    /**
     * Get the facet's function space of the conductor.
     *
     * @return FaceDeg1 of the conductor
     */
    public FaceDeg1 getFD1() {
        return FDcond;
    }

    /**
     * Get the Cell used to compute the dielectric.
     *
     * @return dielectric function space.
     */
    public Cell getDielecCell() {
        return CellDielec;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        GestionnaireTaches.getGestionnaireTaches().setNbCPU(1);
        String meshDir = null;
        try {
            meshDir = new java.io.File(".").getCanonicalPath();
        } catch (IOException ex) {
            Logger.getLogger(TestTwo_Preconditioner.class.getName()).log(Level.SEVERE, null, ex);
        }
        meshDir += "/src/g2elab/mipse/formulationInProgress";
        ImportFlux mesh = new ImportFlux(meshDir + "/magnetodynamic/U_PEEC_DIELECTRIC/Cond_n_Dielec_Surf/FS_T1.DEC");
        /*
         *******************************
         ***  Import du fichier Flux ***
         *******************************
         Nombre de regions importees : 5
         Region 0, Nom : CU, type : 2, 125 elements
         Region 1, Nom : DIELEC, type : 2, 245 elements
         Region 2, Nom : B1, type : 1, 5 elements
         Region 3, Nom : B2, type : 1, 5 elements
         Region 4, Nom : NULL, type : 1, 50 elements
         *******************************
         ***       Fin Import        ***
         *******************************
         */
        TestTwo_Preconditioner solP = new TestTwo_Preconditioner((SurfaceRegion) mesh.getRegion(0), 1 / 1.68e-8, 35e-6, (LineRegion) mesh.getRegion(2), (LineRegion) mesh.getRegion(3),
                (SurfaceRegion) mesh.getRegion(1), 4.7 * eps0,
                (SurfaceRegion) mesh.getRegion(0));

        solP.getTopologie();
        solP.setPtsDeGaussInductifs(9, 4, 4);
        solP.setPtsDeGaussCapacitifs(4, 4);

        FaceDeg1 fd = solP.getFD1();

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();
//        circuitPur.addSourceISimple(400, 125, 130, "t", 1.0, 0.0);
        int indNoeudBord = fd.getNbElement();

        circuitPur.addSourceUSimple(solP.getNbLignes(), indNoeudBord, indNoeudBord + 1, "Source", 240.0, 0.0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

//        solP.integration();
//        solP.checkMatricesSingularities();
//        solP.getZ(10);
        double[][] ib = solP.resolutionIterative(1e6);

        Matrix res = new Matrix(2, solP.getNbLignes());
        for (int i = 0; i < res.getColumnCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i] / 35e-6);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1] / 35e-6);
        }

        exportRes("D:/jsiau/_Backup_Sources/Resultats/Dielectric/Fully_Surf/Test", fd, res, solP.getDielecCell(), solP.getQ(), true);

        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    protected void checkProductFunctions() {
        // Get the full matrix of R
        Basic2D rm = R.getFullMatrix();
        // Set it at the Matrix (GOT) format
        Matrix rM = new Matrix(rm.getRows(), rm.getColumns());
        for (int i = 0; i < rm.getRows(); i++) {
            rM.setRow(i, rm.getArray()[i]);
        }
        // Compute the sparse product
        double[][] m = productSparseFull(R, L);
        // Set it as a MAtrix
        Matrix M2 = new Matrix(rm.getRows(), rm.getColumns());
        for (int i = 0; i < rm.getRows(); i++) {
            M2.setRow(i, m[i]);
        }
        // Compute the full product
        Matrix M = new Matrix(rm.getRows(), rm.getColumns());
        M.mul(rM, L);
        // Check the error
        M2.sub(M);
        System.out.println("Error productSparseFull= " + (M2.norm() / M.norm()));
        //
        //
        // Full product
        M.mul(L, rM);
        // Sparse product
        m = productFullSparse(L, R);
        M2 = new Matrix(rm.getRows(), rm.getColumns());
        for (int i = 0; i < rm.getRows(); i++) {
            M2.setRow(i, m[i]);
        }
        // Check error
        M2.sub(M);
        System.out.println("Error productFullSparse= " + (M2.norm() / M.norm()));
    }

    /**
     * Export the gmsh files to path+"_Q.msh"
     *
     * @param path the path + the name of the export (without extansion)
     * @param fd
     * @param res
     * @param dielec
     * @param Q
     * @param openFile
     */
    public static void exportRes(String path, FaceDeg1 fd, Matrix res, Cell dielec, double Q[], boolean openFile) {

        RealFaceQuantity Jreal = new RealFaceQuantity(fd, res.row(0).subvector(0, fd.getActiveDofCount()));
        RealFaceQuantity Jimag = new RealFaceQuantity(fd, res.row(1).subvector(0, fd.getActiveDofCount()));

        String path_Jreal = path + "_RE.msh";
        String path_Jim = path + "_IM.msh";
        String path_Jmod = path + "_MOD.msh";

        ExportGmshHdiv exportJreal = new ExportGmshHdiv(fd, path_Jreal);
        exportJreal.addQuantity(Jreal, "Jreal");
        ExportGmshHdiv exportJimag = new ExportGmshHdiv(fd, path_Jim);
        exportJimag.addQuantity(Jimag, "Jimag");
        ComplexFaceQuantity J = new ComplexFaceQuantity(fd, res.submatrix(0, 0, 2, fd.getActiveDofCount()));
        exportJreal = new ExportGmshHdiv(fd, path_Jmod);
        exportJreal.addQuantityExportMod(J, "Jmod");

        String path_Q = null;
        if (Q != null) {
            path_Q = path + "_Q.msh";
            ExportGmshCell expQ = new ExportGmshCell(dielec, path_Q);

            RowVector v = new RowVector(dielec.getActiveDofCount());
            // Compute module
            for (int i = 0; i < dielec.getActiveDofCount(); i++) {
                v.setElement(i, Math.hypot(Q[2 * i], Q[2 * i + 1]));
            }
            RealScalarCellQuantity Qmod = new RealScalarCellQuantity(dielec, v);
            expQ.addQuantity(Qmod, "Q");
        }
        String path_Merge = path + "_MERGE.msh";
        try {
            Ecriture merge = new Ecriture(path_Merge);
            merge.ecrire("Merge '" + path_Jreal + "';\n");
            merge.ecrire("Merge '" + path_Jim + "';\n");
            merge.ecrire("Merge '" + path_Jmod + "';\n");
            if (Q != null) {
                merge.ecrire("Merge '" + path_Q + "';\n");
            }
            merge.close();

            if (openFile) {
                File out_file = new File(path_Merge);
                Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file);
            }
        } catch (IOException ex) {
            Logger.getLogger(TestTwo_Preconditioner.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public void setEquiPot(int[] ind) {
        this.equiPot = ind;
    }

    void setAnalyticalIntegrationCapa(boolean fullAnalytical) {
        this.fullCapa = fullAnalytical;
    }

    @Override
    public double[] precond(double[] x, double[] res) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void configuration(double[] config) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public long getMemoryUsed() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setMatrix(Object mat, double[] configuration) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
