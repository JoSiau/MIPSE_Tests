package hMatrix.real.Magnetostatic;

///*
// * To change this license header, choose License Headers in Project Properties.
// * To change this template file, choose Tools | Templates
// * and open the template in the editor.
// */
//package g2elab.mipse.compressionMatricielle.Hmatrix.Tests.Magnetostatic;
//
//import g2elab.mipse.compressionMatricielle.octree.repartitionElements.RepartitionElemNbMaxElem;
//import g2elab.mipse.compressionMatricielle.octree.repartitionElements.RepartitionElemNbMinElem;
//import g2elab.mipse.compressionMatricielle.octree.repartitionElements.RepartitionElemNiveauConstant;
//import g2elab.mipse.meshCore.IO.flux.ImportFlux;
//import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHgrad;
//import g2elab.mipse.meshCore.elements.ElementSetHomogene;
//import g2elab.mipse.meshCore.functionSpace.FunctionSpace;
//import g2elab.mipse.meshCore.functionSpace.GradNodalDeg1;
//import g2elab.mipse.meshCore.functionSpace.Hgrad;
//import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
//import g2elab.mipse.meshCore.gaussValues.GaussPointsPositions;
//import g2elab.mipse.meshCore.gaussValues.GaussPointsValues;
//import g2elab.mipse.meshCore.gaussValues.GaussPointsWeights;
//import g2elab.mipse.meshCore.gaussValues.VectorGaussPointsValues;
//import g2elab.mipse.meshCore.quantity.RealNodalQuantity;
//import g2elab.mipse.method.correctionstrategies.typecorrection.AnalyticalCorrection;
//import g2elab.mipse.method.integrationcorrectionstrategies.IntegrationCorrectionStrategySource;
//import g2elab.mipse.method.integrationcorrectionstrategies.SelfElementFixedGauss;
//import g2elab.mipse.method.numerical.assembly.FiniteElementFormulation;
//import g2elab.mipse.method.numerical.assembly.galerkine.GalerkinIntegralFormulation;
//import g2elab.mipse.method.numerical.assembly.galerkine.GalerkinIntegralFormulationACA;
//import g2elab.mipse.method.numerical.assembly.galerkine.GalerkinIntegralFormulationFMM;
//import g2elab.mipse.method.numerical.assembly.galerkine.GalerkinIntegralFormulationFull;
//import g2elab.mipse.method.numerical.assembly.galerkine.GalerkinIntegralFormulationHCA;
//import g2elab.mipse.method.numerical.tools.kernel.Dipole;
//import g2elab.mipse.method.numerical.tools.kernel.KernelInterface;
//import g2elab.mipse.method.numerical.tools.kernel.MultG;
//import g2elab.mipse.method.numerical.tools.kernel.NegDotDG;
//import g2elab.mipse.method.numerical.tools.storage.StorageBasic2D;
//import g2elab.mipse.method.numerical.tools.storage.StorageFMMGalerkine;
//import g2elab.mipse.method.numerical.tools.storage.StorageFull;
//import g2elab.mipse.method.numerical.tools.storage.StorageHmatrix;
//import g2elab.mipse.method.numerical.tools.storage.StorageSparse;
//import g2elab.mipse.numericalTools.iterativeSolver.ProductReal;
//import g2elab.mipse.numericalTools.iterativeSolver.real.FGMResReal;
//import g2elab.mipse.numericalTools.matrix.real.dense.basic2D.Basic2D;
//import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
//import g2elab.mipse.outils.fichier.Ecriture;
//import g2elab.mipse.outils.fichier.Lecture;
//import g2elab.mipse.outils.multiTaches.GestionnaireTaches;
//import got.matrix.ColumnVector;
//import got.matrix.Matrix;
//import got.matrix.RowVector;
//import java.io.IOException;
//import java.util.Scanner;
//import java.util.logging.Level;
//import java.util.logging.Logger;
//
///**
// *
// * @author jsiau
// */
//public class CEFC_FullPaper1TestBB implements ProductReal {
//
//    final int nDof;
//
//    //
//    // Hmatrix arguments
//    //
//    final int nmin = 60;
//    final int kmax = 50;
//    double eps = 1e-4;
//    int ordre = 3;
//    final double eta = 2.0;
//    final boolean recomp = false;
//    final boolean agglo = false;
//
//    double epsSolver = 1e-11;
//    //
//    // STORAGES
//    //
//    boolean isHmatrix = true;
//    private SparseMatrixRowReal SFE;
//    private StorageHmatrix H;
//    private StorageFMMGalerkine FMM;
//
//    KernelInterface Kernel = new NegDotDG();
//
//    double Qi = 999;
//
//    double muR = Qi + 1;
//
//    final boolean electroStatic;
//
//    int nbGauss = 4;
//
//    private Hgrad Source;
//
//    private Hgrad Target;
//
//    IntegrationCorrectionStrategySource strategy = new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection());
//
//    int nbGaussImproved = 15;
//
//    double mu0 = 4 * Math.PI * 1e-7;
//
//    final String pathMesh;
//
//    public static void main(String[] args) {
//        CEFC_FullPaper1TestBB cefc_FullPaper = new CEFC_FullPaper1TestBB();
//    }
//    private double a;
//
//    /**
//     * Number of gauss point for the Right-hand-side vector integration
//     */
//    int nbGaussSM = 100;// 144
//
//    final boolean ComputeResolution;
//
//    private final String pathOut = "Y:/Resultats/CEFC14/";
//
//    public CEFC_FullPaper1TestBB() {
//        GestionnaireTaches.getGestionnaireTaches().setNbCPU(1);
//        // Sphere Creuse
//        System.out.println("=== CEFC_fullpaper.java ===");
//        System.out.println("Quelle formulation ? \n 1- Electrostatique \n 2- Magnetostatique");
//        Scanner sc = new Scanner(System.in);
//        int f = sc.nextInt();
//
//        boolean chooseSCsize;
//        if (f == 1) {
//            chooseSCsize = true;
//            electroStatic = true;
//        } else {
//            electroStatic = false;
//
//            System.out.println("Quel maillage ? \n 1- Sphere epaisse \n 2- Problem 13");
//            int cc = sc.nextInt();
//            if (cc == 1) {
//                chooseSCsize = true;
//            } else {
//                chooseSCsize = false;
//            }
//        }
//
//        if (chooseSCsize) {
//            if (electroStatic) {
//                System.out.println("Selectionner la taille du maillage:\n0= 4098\n1= 11426\n2= 45156");
//                int scSize = sc.nextInt();
////            4098 12938 26420 48154 79380
//                switch (scSize) {
//                    case 0:
//                        pathMesh = "D:/Meshs/Sphere/SPHERE_" + "4098.DEC";
//                        nDof = 4098;
//                        break;
//                    case 1:
//                        pathMesh = "D:/Meshs/Sphere/SPHERE_" + "11426.DEC";
//                        nDof = 11426;
//                        break;
//                    case 2:
//                        pathMesh = "D:/Meshs/Sphere/SPHERE_" + "45156.DEC";
//                        nDof = 45156;
//                        break;
//                    default:
//                        throw new InternalError("Miss tip ?");
//                }
//            } else {
//                System.out.println("Selectionner la taille du maillage:\n0= 3846\n1= 12938\n2= 48154\n3= 79380");
//                int scSize = sc.nextInt();
////            pathSC = "D:/Meshs/SphereCreuse/SC_" + "3846.DEC";// 3846 12938 26420 48154 79380
//                switch (scSize) {
//                    case 0:
//                        pathMesh = "D:/Meshs/SphereCreuse/SC_" + "3846.DEC";
//                        nDof = 3846;
//                        break;
//                    case 1:
//                        pathMesh = "D:/Meshs/SphereCreuse/SC_" + "12938.DEC";
//                        nDof = 12938;
//                        break;
//                    case 2:
//                        pathMesh = "D:/Meshs/SphereCreuse/SC_" + "48154.DEC";
//                        nDof = 48154;
//                        break;
//                    case 3:
//                        pathMesh = "D:/Meshs/SphereCreuse/SC_" + "79380.DEC";
//                        nDof = 79380;
//                        break;
//                    default:
//                        throw new InternalError("Miss tip ?");
//                }
//            }
//            ComputeResolution = true;
//        } else {
//            nDof = 29362;
//            pathMesh = "D:/Meshs/CEFC14/Problem13/P13_30KP.DEC";
//            ComputeResolution = false; // On assemble uniquement pour le problem 13 !
//
//            System.out.println("Choisissez la précision: \n 0- 1e-2\n 1- 1e-3\n 2- 1e-4");
//
//            int Neps = sc.nextInt();
//            switch (Neps) {
//                case 0:
//                    pathMesh = "Y:/Maillages/SC_" + "200P.DEC";
//                    nDof = 252;
//                    break;
//                case 1:
//                    this.eps = 1e-3;
//                    this.ordre = 2;
//                    break;
//                case 2:
//                    this.eps = 1e-4;
//                    this.ordre = 3;
//                    break;
//                case 3:
//                    this.eps = 1e-5;
//                    this.ordre = 4;
//                    break;
//                case 4:
//                    this.eps = 1e-6;
//                    this.ordre = 5;
//                    break;
//                default:
//                    throw new InternalError("Miss tip ?");
//            }
//            System.out.println("eps= "+eps+"\t ordre= "+ordre);
//        }
//        int cpt = 0;
//        System.out.println("Choisissez quel compression: \n " + (cpt++) + "- Full\n " + (cpt++) + "- ACA\n " + (cpt++) + "- HCA" + "\n " + (cpt++) + "- FMM"  + "\n " + (cpt++) + "- Compare all");
//        int cx = sc.nextInt();
//        switch (cx) {
//            case 0:
//                ResFullSC();
//                break;
//            case 1:
//                ResACASC();
//                break;
//            case 2:
//                ResHCASC();
//                break;
//            case 3:
////                if (electroStatic) {
//                    ResFMMSC();
////                } else {
////                    autoCompareSC();
////                }
//                break;
//            case 4:
////                if (electroStatic) {
//                    autoCompareSC();
//                    break;
////                }
//            default:
//                throw new InternalError("Miss tip ?");
//        }
//
//        GestionnaireTaches.getGestionnaireTaches().stop();
//    }
//
//    protected void autoCompareSC() {
//        //
//        // INITIALISATION DES DONNEES
//        //
//        ColumnVector resFull = new ColumnVector(nDof);
//        ColumnVector resACA;
//        ColumnVector resFMM;
//        ColumnVector resHCA;
//
//        //
//        // FULL
//        //
//        System.out.println("========== Résolution Pleine ==========");
////        try {
////            String doss = electroStatic ? "electrostatic" : "magnetostatic_sc";
////            Lecture fic = new Lecture(pathOut + doss + "/vecResRes_" + nDof + ".out");
////            System.out.println("Reading the file");
////            for (int i = 0; i < nDof; i++) {
////                resFull.setElement(i, fic.lireLigneDouble());
////            }
////            fic.close();
////        } catch (IOException ex) {
//            System.out.println("Compute the full matrix");
//            resFull = ResFullSC();
////        }
//        double normFull = resFull.norm();
//        //
//        // ACA
//        //
//        System.out.println("========== Resolution ACA ==========");
////        try {
////            String doss = electroStatic ? "electrostatic" : "magnetostatic_sc";
////            Lecture fic = new Lecture(pathOut + doss + "/ACA_" + nDof + ".out");
////            System.out.println("Reading the file");
////            resACA = new ColumnVector(nDof);
////            for (int i = 0; i < nDof; i++) {
////                resACA.setElement(i, fic.lireLigneDoubleVecteur(' ')[1]);
////            }
////            fic.close();
////        } catch (IOException ex) {
//        resACA = ResACASC();
////        }
//        resACA.sub(resFull);
//        System.out.println("Erreur relative ACA = " + resACA.norm() / normFull);
//
//        //
//        // FMM
//        //
//        //if (electroStatic) {
//            System.out.println("========== Resolution FMM ==========");
//            resFMM = ResFMMSC();
//
//            resFMM.sub(resFull);
//            System.out.println("Erreur relative FMM = " + resFMM.norm() / normFull);
//        //}
//        //
//        // HCA
//        //
//        System.out.println("========== Resolution HCA ==========");
////        try {
////            String doss = electroStatic ? "electrostatic" : "magnetostatic_sc";
////            Lecture fic = new Lecture(pathOut + doss + "/HCA_" + nDof + ".out");
////            System.out.println("Reading the file");
////            resHCA = new ColumnVector(nDof);
////            for (int i = 0; i < nDof; i++) {
////                resHCA.setElement(i, fic.lireLigneDoubleVecteur(' ')[1]);
////            }
////            fic.close();
////        } catch (IOException ex) {
//        resHCA = ResHCASC();
////        }
//
//        resHCA.sub(resFull);
//        System.out.println("Erreur relative HCA = " + resHCA.norm() / normFull);
//
//    }
//
//    protected ColumnVector ResFullSC() {
//        double deb, fin;
//        ImportFlux mesh = new ImportFlux(pathMesh);
//        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegions(0).getElementSet();
//
//        System.out.println("nbNoeuds= " + ES.getNbNoeud() + "\t nbElmts= " + ES.getNbElement());
//
//        String doss;
//        if (electroStatic) {
//            doss = "electrostatic";
//            System.out.println("Electrostatic");
//            Target = new NodalDeg1(ES);
//            Source = new NodalDeg1(ES);
//            strategy = new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection());
//            Kernel = new MultG();
//        } else {
//            System.out.println("Magnetostatic");
//            doss = "magnetostatic_sc";
//            Target = new NodalDeg1(ES);
//            Source = new GradNodalDeg1(ES);
//        }
////        getSecondMembre(Target, nDof);
//
//        ColumnVector res = new ColumnVector(nDof);
//
//        GalerkinIntegralFormulationFull f = new GalerkinIntegralFormulationFull(Target, Source, Kernel, strategy, nbGauss, new Basic2D(new double[1][1]));
//        deb = System.nanoTime();
//        f.assembly();
//        fin = System.nanoTime();
//        System.out.println("Time to assembly the Full-Matrix: " + (fin - deb) / 1e9);
//
//        Basic2D Mgot = ((StorageBasic2D) f.getStore()).getBasic2D();
//
//        if (ComputeResolution) {
//            double M[][] = Mgot.getArray();
//            // Mgot.scale(Qi);
//            // Un scale tres mauvais.
//            if (!electroStatic) {
//                this.computeEF(Target, nbGauss);
//
//                for (int i = 0; i < M.length; i++) {
//                    for (int j = 0; j < M[i].length; j++) {
//                        M[i][j] *= Qi;
//                    }
//                }
//
//                int ai[][] = SFE.getAi();
//                double vi[][] = SFE.getAv();
//                for (int i = 0; i < ai.length; i++) {
//                    for (int j = 0; j < ai[i].length; j++) {
//                        M[i][ai[i][j]] += vi[i][j];
//                    }
//                }
//            }
//
//            deb = System.currentTimeMillis();
//            double[] secondMembre = getSecondMembre(Target, nDof);
//
//            FGMResReal FGMRes = new FGMResReal(Mgot, null);
//            FGMRes.setInfoResolution(new double[]{1000, epsSolver, 1, -50});
//            
//            //SolverLUBasic2D lu = new SolverLUBasic2D(M, 1e-40, false, true, true);
//            //res.set(lu.solveLU(secondMembre));
//            System.out.println("JERENTREDEDANS");
//            res.set(FGMRes.solve(new double[secondMembre.length], secondMembre));
//            fin = System.currentTimeMillis();
//            System.out.println("Solving time= " + (fin - deb) * 1e-3 + " sec.");
////             try {
////             Lecture ficLref = new Lecture(pathOut+doss+"/vecResRes_" + nDof + ".out");// Vecteur resultat de la resolution
////
////             for (int i = 0; i < nDof; i++) {
////             res.setElement(i, ficLref.lireLigneDouble());
////             }
////             ficLref.close();
////             } catch (IOException ex) {
////             System.out.println("The full resolution is not computed yet.");
////             throw new InternalError();
////             }
//            RealNodalQuantity phi = new RealNodalQuantity((Hgrad) Target, res.transpose());
//            /*
//             * Export
//             */
//            ExportGmshHgrad exportPhiH = new ExportGmshHgrad((Hgrad) Source, pathOut + doss + "/phi_" + nDof + ".msh");
//            exportPhiH.addQuantity(phi, "phi");
//
//            try {
//                Ecriture ficE = new Ecriture(pathOut + doss + "/vecResRes_" + nDof + ".out");
//                System.out.println("Ecriture de la matrice a : " + (pathOut + doss + "/vecResRes_" + nDof + ".out"));
//                for (int i = 0; i < nDof; i++) {
//                    ficE.ecrire(res.getElement(i) + "\n");
//                }
//                ficE.close();
//            } catch (IOException ex) {
//                Logger.getLogger(CEFC_FullPaper1TestBB.class.getName()).log(Level.SEVERE, null, ex);
//            }
//            return res;
//        } else {
//            double[] b = getRandomVector();
//            double[] r = Mgot.product(b, new double[nDof]);
//
//            try {
//                Ecriture ficE = new Ecriture(pathOut + "magnetostatic_p13/vecRes_" + nDof + ".out");
//                for (int i = 0; i < nDof; i++) {
//                    ficE.ecrire(r[i] + "\n");
//                }
//                ficE.close();
//            } catch (IOException ex) {
//                Logger.getLogger(CEFC_FullPaper1TestBB.class.getName()).log(Level.SEVERE, null, ex);
//            }
//
//            return new ColumnVector(r);
//        }
//    }
//
//    /**
//     * (COMPUTE AND SAVE) OR READ A RANDOM VECTOR. Only for the Problem13
//     *
//     * @return
//     */
//    protected double[] getRandomVector() {
//        double[] resFull = new double[nDof];
//        try {
//            Lecture ficLref = new Lecture(pathOut + "magnetostatic_p13/vecRdm_" + nDof + ".out");// Vecteur resultat de la resolution
//            System.out.println("Reading a rdm vector");
//            for (int i = 0; i < nDof; i++) {
//                resFull[i] = ficLref.lireLigneDouble();
//            }
//            ficLref.close();
//        } catch (IOException ex) {
//            System.out.println("Computing a new random vector !");
//
//            for (int i = 0; i < nDof; i++) {
//                resFull[i] = Math.random() * 100;
//            }
//            try {
//                Ecriture fic = new Ecriture(pathOut + "magnetostatic_p13/vecRdm_" + nDof + ".out");
//                for (int i = 0; i < nDof; i++) {
//                    fic.ecrire(resFull[i] + "\n");
//                }
//                fic.close();
//            } catch (IOException ex1) {
//                System.out.println("Couldn't write !");
//                Logger.getLogger(CEFC_FullPaper1TestBB.class.getName()).log(Level.SEVERE, null, ex1);
//            }
//        }
//        return resFull;
//    }
//
//    protected double[] getSecondMembre(Hgrad Target, int n) {
//        double[] secondMembre;
//        if (electroStatic) {
//            // V0 = 1
//            RowVector V0 = new RowVector(n);
//            //*
//            V0.setAllElements(1);
//            /*/
//             for (int i = 0; i < n; i++) {
//             V0.setElement(i, Math.random()/1000);
//             }
//             //*/
//            // int_\Omega V0 
//            RealNodalQuantity V0p1 = new RealNodalQuantity(Target, V0);
//            RowVector bg1 = (RowVector) V0p1.projWithExplicitDof(Target, nbGaussSM);
//
//            secondMembre = new double[n];
//            bg1.get(secondMembre);
//
//        } else {
//            a = 2.5;
//            // Compute phi0(X) = -a*x with X = (x, y, z);
//            secondMembre = new double[n];
//            for (int i = 0; i < n; i++) {
//                secondMembre[i] = -a * Target.getElementSet().getLocalNodes()[i].getCoord(0);
//            }
//
//            // Compute \int_\Omega a_i(r) phi0(r) dr
//            RealNodalQuantity V0p1 = new RealNodalQuantity(Target, new RowVector(secondMembre));
//            RowVector bg1 = (RowVector) V0p1.projWithExplicitDof(Target, nbGauss);
//            System.out.println("bg1= " + bg1);
//            // Recupere le resultat
//            secondMembre = new double[n];
//            bg1.get(secondMembre);
//        }
//        return secondMembre;
//    }
//
//    protected void checkError(RowVector resH) {
//        int n = resH.getColumnCount();
//        RowVector resFull;
//
//        String doss = ComputeResolution ? (electroStatic ? "electrostatic" : "magnetostatic_sc") : "magnetostatic_p13";
//
//        try {
//            Lecture ficLref = new Lecture(pathOut + doss + "/vecResRes_" + n + ".out");// Vecteur resultat de la resolution
//
//            System.out.println("Lecture de la matrice : " + pathOut + doss + "/vecResRes_" + n + ".out");
//
//            resFull = new RowVector(n);
//            for (int i = 0; i < n; i++) {
//                resFull.setElement(i, ficLref.lireLigneDouble());
//            }
//            ficLref.close();
//        } catch (IOException ex) {
//            System.out.println("The full resolution is not computed yet.");
//            return;
//        }
//        RowVector tmp = new RowVector(n);
//        tmp.sub(resH, resFull);
//        System.out.println("Erreur relative= " + tmp.norm() / resFull.norm());
//    }
//
//    protected ColumnVector ResHCASC() {
//        double deb, fin;
//
//        ImportFlux mesh = new ImportFlux(pathMesh);
//        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegions(0).getElementSet();
//
//        System.out.println("nbNoeuds= " + ES.getNbNoeud() + "\t nbElmts= " + ES.getNbElement());
//
//        String doss;
//        if (electroStatic) {
//            doss = "electrostatic";
//            System.out.println("Electrostatic");
//            Target = new NodalDeg1(ES);
//            Source = new NodalDeg1(ES);
//            strategy = new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection());
//            Kernel = new MultG();
//        } else {
//            System.out.println("Magnetostatic");
//            doss = "magnetostatic_sc";
//            Target = new NodalDeg1(ES);
//            Source = new GradNodalDeg1(ES);
//        }
//
//        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(Target, Source, Kernel, strategy, nbGauss, nbGauss,
//                eps, kmax, nmin, ordre, eta, recomp);
//        deb = System.nanoTime();
//        f.assembly();
//        H = (StorageHmatrix) f.getStore();
//        if(agglo && ComputeResolution){
//            H.Agglomerate(new TruncationControl("rel", eps));
//        }
//        fin = System.nanoTime();
//        System.out.println("Time to assembly the Hmatrix HCA: " + (fin - deb) / 1e9);
//        System.out.println("HCA memory use= " + H.getMemoryUsed());
//
//        if (ComputeResolution) {
//            isHmatrix = true;
//            if (!electroStatic) {
//                computeEF(Target, nbGauss);
//                H.scale(Qi);
//            }
//            FGMResReal solver = new FGMResReal(this, null);
//            solver.setInfoResolution(new double[]{1000, epsSolver, 1, -50});
//            double[] x = new double[nDof];
//            double[] secondMembre = getSecondMembre(Target, nDof);
//            deb = System.nanoTime();
//            solver.solve(x, secondMembre);
//            fin = System.nanoTime();
//            System.out.println("Time to solve the Hmatrix HCA: " + (fin - deb) / 1e9);
//            RowVector resH = new RowVector(x);
//            RealNodalQuantity phiH = new RealNodalQuantity((Hgrad) Target, resH);
//
//            /*
//             * Export
//             */
//            ExportGmshHgrad exportPhiH = new ExportGmshHgrad((Hgrad) Source, pathOut + doss + "/phiHCA_" + nDof + ".msh");
//            exportPhiH.addQuantity(phiH, "phiH");
//
//            checkError(resH);
//
//            return resH.transpose();
//        } else {
//
//            double[] b = getRandomVector();
//            double[] r = H.product(b, new double[nDof]);
//
//            this.checkError(new RowVector(r));
//
//            deb = System.currentTimeMillis();
//            H.CheckError();
//            fin = System.currentTimeMillis();
//            System.out.println("Time to compute the approximated error= " + (fin - deb) * 1e-3);
//
//            return new ColumnVector(r);
//        }
//
//    }
//
//    protected void computeEF(FunctionSpace Target, int nbGauss) {
//        System.out.print("Compute the mass matrix (EF): ");
//        long deb = System.currentTimeMillis();
//        FiniteElementFormulation FE = new FiniteElementFormulation(Target);
//        // Nombre de point de gauss pour integration element finis 4  minimum  (on prends les pt de gauss des sources)
//        FE.assembly(nbGauss);
//        System.out.println((System.currentTimeMillis() - deb) + " ms");
//
//        SFE = ((StorageSparse) FE.getStore()).getMatrixPrecond(null);
//    }
//
//    /**
//     * Calcul le champ H au point dans l'air sans correction annalytique Cette
//     * methode est pour les points qui situent loins de region active
//     *
//     * @param MPG
//     * @param cibleCood
//     * @return
//     */
//    public ColumnVector FieldPointAir(VectorGaussPointsValues MPG, double[] cibleCood) {
//        int nbGaussSource = MPG.getNbGauss();
//        ColumnVector fieldPoint = new ColumnVector(3);
//        Dipole dipole = new Dipole();
//        // Recuperation maillage
//        ElementSetHomogene mesh = MPG.getElementSet();
//
//        // Gauss point postion computation
//        GaussPointsPositions GPP = new GaussPointsPositions(mesh, nbGaussSource);
//
//        // Gauss points weight computation
//        GaussPointsWeights GPWSource = new GaussPointsWeights(mesh, nbGaussSource);
//
//        // Values de kernel aux points de Gauss
//        GaussPointsValues kerGaussPoints = dipole.getKernelGPV(mesh, MPG, new ColumnVector(cibleCood), GPP);
//        //Passe au repere global
//        kerGaussPoints.scale(GPWSource);
//
//        fieldPoint = kerGaussPoints.getInfluenceM();
//
//        return fieldPoint;
//    }
//
//    private ColumnVector ResACASC() {
//        double deb, fin;
//
//        ImportFlux mesh = new ImportFlux(pathMesh);
//        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegions(0).getElementSet();
//
//        System.out.println("nbNoeuds= " + ES.getNbNoeud() + "\t nbElmts= " + ES.getNbElement());
//
//        String doss;
//        if (electroStatic) {
//            System.out.println("Electrostatic");
//            doss = "electrostatic";
//            Target = new NodalDeg1(ES);
//            Source = new NodalDeg1(ES);
//            strategy = new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection());
//            Kernel = new MultG();
//        } else {
//            System.out.println("Magnetostatic");
//            doss = "magnetostatic_sc";
//            Target = new NodalDeg1(ES);
//            Source = new GradNodalDeg1(ES);
//        }
//
//        deb = System.nanoTime();
//        GalerkinIntegralFormulationACA f = new GalerkinIntegralFormulationACA(Target, Source, Kernel, strategy, nbGauss,
//                eps, kmax, nmin, eta, recomp);
//        f.assembly();
//        H = (StorageHmatrix) f.getStore();
//        if(agglo && ComputeResolution){
//            H.Agglomerate(new TruncationControl("rel", eps));
//        }
//        fin = System.nanoTime();
//        System.err.println("Time to assembly the Hmatrix ACA: " + (fin - deb) / 1e9);
//        System.out.println("ACA memory use= " + H.getMemoryUsed());
//
//        if (ComputeResolution) {
//            isHmatrix = true;
//            if (!electroStatic) {
//                computeEF(Target, nbGauss);
//                H.scale(Qi);
//            }
//            FGMResReal solver = new FGMResReal(this, null);
//            solver.setInfoResolution(new double[]{1000, epsSolver, 1, -50});
//            double[] x = new double[nDof];
//            double[] secondMembre = getSecondMembre(Target, nDof);
//            deb = System.nanoTime();
//            solver.solve(x, secondMembre);
//            fin = System.nanoTime();
//            System.err.println("Time to solve the Hmatrix ACA: " + (fin - deb) / 1e9);
//            RowVector resH = new RowVector(x);
//            RealNodalQuantity phiH = new RealNodalQuantity((Hgrad) Target, resH);
//
//            /*
//             * Export
//             */
//            ExportGmshHgrad exportPhiH = new ExportGmshHgrad((Hgrad) Source, pathOut + doss + "/phiACA_" + nDof + ".msh");
//            exportPhiH.addQuantity(phiH, "phiH");
//
//            checkError(resH);
//
//            return resH.transpose();
//        } else {
//            double[] b = getRandomVector();
//            double[] r = H.product(b, new double[nDof]);
//
//            this.checkError(new RowVector(r));
//
//            deb = System.currentTimeMillis();
//            H.CheckError();
//            fin = System.currentTimeMillis();
//            System.out.println("Time to compute the approximated error= " + (fin - deb) * 1e-3);
//
//            return new ColumnVector(r);
//        }
//    }
//
//    private ColumnVector ResFMMSC() {
//        double deb, fin;
//
//        ImportFlux mesh = new ImportFlux(pathMesh);
//        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegions(0).getElementSet();
//
//        System.out.println("nbNoeuds= " + ES.getNbNoeud() + "\t nbElmts= " + ES.getNbElement());
//
//        String doss;
//        if (electroStatic) {
//            System.out.println("Electrostatic");
//            doss = "electrostatic";
//            Target = new NodalDeg1(ES);
//            Source = new NodalDeg1(ES);
//            strategy = new SelfElementFixedGauss(nbGauss, new AnalyticalCorrection());
//            Kernel = new MultG();
//        } else {
//            System.out.println("Magnetostatic");
//            doss = "magnetostatic_sc";
//            Target = new NodalDeg1(ES);
//            Source = new GradNodalDeg1(ES);
//            computeEF(Target, nbGauss);
//        }
//
//
//        GalerkinIntegralFormulationFMM f = new GalerkinIntegralFormulationFMM(Target, Source, Kernel, strategy, nbGauss, nbGauss,
//                new RepartitionElemNbMaxElem(nmin), 1.0, 5, 1, 1);
//        deb = System.nanoTime();
//        f.assembly();
//        FMM = (StorageFMMGalerkine) f.getStore();
//        fin = System.nanoTime();
//        System.out.println("Time to assembly the FMM: " + (fin - deb) / 1e9);
//        System.out.println("FMM memory= " + FMM.getMemoryUsed());
//
//
//        if (ComputeResolution) {
//            isHmatrix = false;
//            FGMResReal solver = new FGMResReal(this, null);
//            solver.setInfoResolution(new double[]{1000, epsSolver, 1, -50});
//            double[] x = new double[nDof];
//            double[] secondMembre = getSecondMembre(Target, nDof);
//            deb = System.nanoTime();
//            solver.solve(x, secondMembre);
//            fin = System.nanoTime();
//            System.err.println("Time to solve the FMM: " + (fin - deb) / 1e9);
//            RowVector resH = new RowVector(x);
//            RealNodalQuantity phiH = new RealNodalQuantity((Hgrad) Target, resH);
//
//            /*
//             * Export
//             */
//            ExportGmshHgrad exportPhiH = new ExportGmshHgrad((Hgrad) Source, pathOut + doss + "/phiFMM_" + nDof + ".msh");
//            exportPhiH.addQuantity(phiH, "phiFMM");
//
//            checkError(resH);
//
//            return resH.transpose();
//        } else {
//            double[] b = getRandomVector();
//            double[] r = FMM.product(b, new double[nDof]);
//            return new ColumnVector(r);
//        }
//    }
//
//    @Override
//    public double[] product(double[] x, double[] res) {
//        if (res == null) {
//            res = new double[nDof];
//        }
//        if (!electroStatic) {
//            res = this.SFE.product(x, res);
//        }
//        if (isHmatrix) {
//            res = this.H.product(x, res);
//        } else {
//            double tmp[] = this.FMM.product(x, new double[nDof]);
//            for (int i = 0; i < tmp.length; i++) {
//                res[i] += !electroStatic? Qi*tmp[i] : tmp[i];
//            }
//        }
//        return res;
//    }
//
//}
