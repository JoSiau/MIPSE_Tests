/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHdiv;
import g2elab.mipse.meshCore.IO.paraview.ExportVtkHdiv;
import g2elab.mipse.meshCore.contraints.ExternalDriveFaceConstraint;
import g2elab.mipse.meshCore.contraints.FaceConstraint;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.elements.ElementVolSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.quantity.ComplexFaceQuantity;
import g2elab.mipse.meshCore.quantity.RealFaceQuantity;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.numericalTools.matrix.real.sparse.rowReal.SparseMatrixRowReal;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class LoopAntenna {

    PEECwithCapa_Full solPEEC;

    PEEC_RLMPC_ConducEquiv solP;

    public LoopAntenna(double f) {
        int nbCPU = 0;
        System.out.println("Entrez:\n 0: RLMC \n 1: RLMPC \n 2: Comparaison ");
        Scanner sc = new Scanner(System.in);
        nbCPU = sc.nextInt();

        File fil = new File("");
        String path = fil.getAbsolutePath();
        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/LOOPANTENNA_DIELEC_VOL.DEC";

        if (nbCPU == 0) {
            this.RLMC(f, file);
        } else if (nbCPU == 1) {
            this.RLMPC(f, file);
        } else {
            Matrix r1 = this.RLMC(f, file);
            Matrix r2 = this.RLMPC(f, file);

            Matrix L = solPEEC.getL();

            Matrix Lp[] = solP.getL();

            System.out.println("L.size= " + L.getRowCount() + " , " + L.getColumnCount());
            System.out.println("L2.size= " + Lp[1].getRowCount() + " , " + Lp[1].getColumnCount());
            System.out.println("norm(L[0])= " + Lp[0].norm() + "norm(L[1])= " + Lp[1].norm());
            L.scale(2 * Math.PI * f);
            Lp[1].sub(L);
            System.out.println("Erreur relative= " + Lp[1].norm() / L.norm());

            SparseMatrixRowReal R = solPEEC.getR();
            SparseMatrixRowReal Rp[] = solP.getR();

            double[] b = new double[R.getColumns()];
            for (int i = 0; i < b.length; i++) {
                b[i] = Math.random();
            }

            double[] v1 = new double[R.getRows()];
            v1 = R.product(b, v1);

            double[] v2 = new double[R.getRows()];
            v2 = Rp[0].product(b, v2);

            ColumnVector rv1 = new ColumnVector(v1);
            ColumnVector rv2 = new ColumnVector(v2);
            rv2.sub(rv1);
            System.out.println("Erreur relative PMV de R =" + rv2.norm() / rv1.norm());

            Matrix P = solPEEC.getP();
            P.scale(-1 / (2 * Math.PI * f));

            ecrireMat(P, "d:/P.out");
            Matrix P2[] = solP.getPe();
            ecrireMat(P2[0], "d:/PeRe.out");
            ecrireMat(P2[1], "d:/PeIm.out");

            System.out.println("norm(Pe_Re)= " + P2[0].norm() + "\t norm(Pe_Im)= " + P2[1].norm() + "\t norm(P)= " + P.norm());

            Matrix Q = P2[1].copy();
            Q.sub(P);
            System.out.println("Erreur relative de P= " + Q.norm() / P.norm());

            Cell G = solPEEC.getBorderCell();
            Cell Ge = solP.getGammaExt();

            G.getElementSet().plotElementSet("d:/G.vtk");
            Ge.getElementSet().plotElementSet("d:/Ge.vtk");

            if (G.getElementSet() == Ge.getElementSet()) {
                System.out.println("Les ElementSet sont identiques !!!");
            } else {
                boolean identical = true;
                for (int i = 0; i < G.getNbElement(); i++) {
                    if (G.getElementSet().getElements(i).compareElement(Ge.getElementSet().getElements(i))) {
                    } else {
                        System.out.println("Les elements " + i + " ne sont pas identiques !\nG= " + G.getElementSet().getElements(i) + "\nGe= " + Ge.getElementSet().getElements(i));
                        identical = false;
                    }
                }
                if (!identical) {
                    System.out.println("Les elements de set ne sont pas identiques !!!");
                } else {
                    System.out.println("Les elements de set sont identiques !!!");
                }
            }
            P2[1].swapColumns(1918, 1921);
            P2[1].swapRows(1918, 1921);

            P2[1].swapColumns(1923, 1924);
            P2[1].swapRows(1923, 1924);
            P2[1].sub(P);
            System.out.println("Erreur relative 2 de P = " + P2[1].norm() / P.norm());

            System.out.println("\n*************************************************\n");

            FaceDeg1 F = solPEEC.getFunctionSpace();
            FaceDeg1 Fe = solP.getFunctionSpace();
            if (F.getElementSet() == Fe.getElementSet()) {
                System.out.println("Les ElementSet sont identiques !!!");
            } else {
                boolean identical = true;
                for (int i = 0; i < F.getNbElement(); i++) {
                    if (F.getElementSet().getElements(i).compareElement(Fe.getElementSet().getElements(i))) {
                    } else {
                        System.out.println("Les elements " + i + " ne sont pas identiques !\nF= " + F.getElementSet().getElements(i) + "\nFe= " + Fe.getElementSet().getElements(i));
                        identical = false;
                    }
                }
                if (!identical) {
                    System.out.println("Les elements de set ne sont pas identiques !!!");
                } else {
                    System.out.println("Les elements de set sont identiques !!!");
                }
            }
            // 
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.out.println("Entrez la fréquence: ");
        Scanner sc = new Scanner(System.in);
        double f = sc.nextDouble();
        LoopAntenna lp = new LoopAntenna(f);
        GestionnaireTaches.getGestionnaireTaches().stop();
    }

    protected Matrix RLMC(double f, String file) {

        ImportFlux IF = new ImportFlux(file);
        int[] indices = {2, 3, 1};
        ExternalDriveFaceConstraint border = new ExternalDriveFaceConstraint(IF.getRegion(1));
        ExternalDriveFaceConstraint FluxNeg = new ExternalDriveFaceConstraint(IF.getRegion(2));
        ExternalDriveFaceConstraint FluxPos = new ExternalDriveFaceConstraint(IF.getRegion(3));
        FaceConstraint[] constraints = new FaceConstraint[]{border, FluxPos, FluxNeg};

        ElementSetHomogene mesh = new ElementVolSetHomogene(IF.getRegion(0).getElementSet().getElements());
        FaceDeg1 FS = new FaceDeg1(mesh, constraints);
        int nbActiveDof = FS.getActiveDofCount();
        System.out.println("nbActiveDof= " + nbActiveDof);

        solPEEC = new PEECwithCapa_Full(FS, 2);
        // Affectation despoints de Gauss
        solPEEC.setPtsDeGaussInductifs(15, 4);
        solPEEC.setPtsDeGaussCapacitifs(4, 4);

        // Construction du circuit electrique exterieur
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        solPEEC.setRho(2.836e-8);
        int[] num = FS.getImplicitConstraintCount();
        int nbE = mesh.getNbElement();
        int nbCapa = num[0];
        int num0 = num[1];
        int num1 = num[2];

        int nbBranches = solPEEC.getNbLignes();
////     noeud commun des capas (noeud GND)   
        int indCapa = 1;
        int C1 = nbE + nbCapa + num0 + num1 + 1 + indCapa;
        int C2 = C1 + 1;
////        
        for (int in = 0; in < num0; in++) {
            circuitPur.addSourceUSimple(nbBranches + in, nbE + nbCapa + 1 + in, C1, "CommunS1", 0, 0.0);
        }
//
        for (int out = 0; out < num1; out++) {
            circuitPur.addSourceUSimple(nbBranches + num0 + out, nbE + nbCapa + 1 + num0 + out, C2, "CommunS2", 0, 0.0);
        }

        circuitPur.addSourceISimple(nbBranches + num0 + num1, C2, C1, "SourceI", 1, 0.0);

        circuitPur.finSaisie(); // Methode a appeler necessairement a la fin de la saisie

        circuitPur.toString();
        // Affectation du circuit exterieur
        solPEEC.setCircuitElectrique(circuitPur);

        // Resolution
        long tres0 = System.currentTimeMillis();
        double ib[][] = solPEEC.resolutionDirecte(f);
        long tres1 = System.currentTimeMillis();
        System.out.println("Temps de resolution:" + (tres1 - tres0));
        // On reorganise les dof   
        Matrix res = new Matrix(2, FS.getActiveDofCount());
        for (int i = 0; i < FS.getActiveDofCount(); i++) {
            res.setElement(0, solPEEC.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solPEEC.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }

        double[] z = new double[]{ib[1][ib.length - 2], ib[1][ib.length - 1]};
        System.out.println("Z(" + f + ")= " + z[0] + " + i* " + z[1]);
        System.out.println("|Z|= " + Math.hypot(z[0], z[1]));
        System.out.println("I= " + ib[0][ib.length - 2] + " +i* " + ib[0][ib.length - 1]);

        // On reorganise les dof   
        RealFaceQuantity Jreal = new RealFaceQuantity(FS, res.row(0));
        RealFaceQuantity Jimag = new RealFaceQuantity(FS, res.row(1));

        ExportGmshHdiv exportJreal = new ExportGmshHdiv(FS, "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMC_f" + f + "_nDof" + nbActiveDof + "_Re.msh");
        exportJreal.addQuantity(Jreal, "Jreal");

        ExportGmshHdiv exportJimag = new ExportGmshHdiv(FS, "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMC_f" + f + "_nDof" + nbActiveDof + "_Im.msh");
        exportJimag.addQuantity(Jimag, "Jimag");

        ComplexFaceQuantity J = new ComplexFaceQuantity(FS, res);
        exportJreal = new ExportGmshHdiv(FS, "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMC_f" + f + "_nDof" + nbActiveDof + "_Jmod.msh");
        exportJreal.addQuantityExportMod(J, "Jmod");

        ExportVtkHdiv ez = new ExportVtkHdiv(FS, "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMC_f" + f + "_nDof" + nbActiveDof + ".vtk");
        ez.addQuantity(J, "J");

        String output_file = "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMC_f" + f + "_nDof" + nbActiveDof + "_Jmod.msh";
        File out_file = new File(output_file);
        File out_file1 = new File("D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMC_f" + f + "_nDof" + nbActiveDof + "_Re.msh");
        File out_file2 = new File("D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMC_f" + f + "_nDof" + nbActiveDof + "_Im.msh");
        try {
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);
        } catch (IOException ex) {
            Logger.getLogger(LoopAntenna.class.getName()).log(Level.SEVERE, null, ex);
        }

        return res;
    }

    protected Matrix RLMPC(double f, String file) {

        /*
         Nombre de regions importees : 4
         Region 0, Nom : LOOPANTENA, type : 3, 3907 elements
         Region 1, Nom : FLUXNULL, type : 2, 3160 elements
         Region 2, Nom : FLUXNEG, type : 2, 8 elements
         Region 3, Nom : FLUXPOS, type : 2, 8 elements
         */
        ImportFlux IF = new ImportFlux(file);
        VolumeRegion rV[] = new VolumeRegion[]{(VolumeRegion) IF.getRegion(0)};
        int nbEl = rV[0].getElementSet().getNbElement();
        double Sigma[][] = new double[2][nbEl];// Size = [2][nbElmttotal]
        // Conducteur
        double Re = 1 / 2.836e-8;
        double Im = PEEC_RLMPC_ConducEquiv.eps0;
        for (int i = 0; i < nbEl; i++) {
            Sigma[0][i] = Re;// Conductivite du cuivre
            Sigma[1][i] = Im;// Epsilon0
        }

        solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma, (SurfaceRegion) IF.getRegion(3), (SurfaceRegion) IF.getRegion(2));
        solP.setPtsDeGaussInductifs(15, 4);
        solP.setPtsDeGaussCapacitifs(3, 3);
        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        FaceDeg1 Omega = solP.getFunctionSpace();
        System.out.println("Activedofcount= " + Omega.getActiveDofCount());

        ImportFlux ifl = new ImportFlux(file);
        SurfaceRegion B1 = (SurfaceRegion) ifl.getRegion(3);
        SurfaceRegion B2 = (SurfaceRegion) ifl.getRegion(2);
        Cell G = solP.getGammaExt();
//        System.out.println("\n Bornes check");
//        System.out.println("B1");
//        for (int i = 0; i < 2; i++) {
//            if(G.getElementSet().getElements(i+1918).compareElement(B1.getElementSet().getElements(i))){
//                
//            }else{
//                System.out.println("Borne1: Les elements "+i+" ne sont pas identiques !\nB1= "+B1.getElementSet().getElements(i)+"\nG= "+G.getElementSet().getElements(i+1918));
//            }   
//        }
//        System.out.println("B2");
//        for (int i = 0; i <2; i++) {
//            if(G.getElementSet().getElements(i+1922).compareElement(B2.getElementSet().getElements(i))){
//                
//            }else{
//                System.out.println("Borne2: Les elements "+i+" ne sont pas identiques !\nB1= "+B2.getElementSet().getElements(i)+"\nG= "+G.getElementSet().getElements(i+1922));
//            }   
//        }
//        
        /*
         CONNEXION DES BORNES
         */
//        int indConsB = solP.getIndexConstraintBorne();
//        System.out.println("indConsB= "+indConsB);
//        int nbB1 = Omega.getImplicitConstraintCount()[indConsB];
//        int indB1 = Omega.getNbElement() + 1;
//        System.out.println("indB1= "+indB1);
//        int nbB2 = Omega.getImplicitConstraintCount()[indConsB+1];
//        int indB2 = indB1 + nbB1;
//        System.out.println("indB2= "+indB2);
//        int indCommun1 = indB2 + nbB2;
//        int indCommun2 = indCommun1 + 1;
//        
//        for (int i = 0; i < nbB1; i++) {
//            circuitPur.addSourceUSimple(Omega.getActiveDofCount() + i, indB1 + i, indCommun1, "Commun1", 0, 0);
//        }
//        for (int i = 0; i < nbB2; i++) {
//            circuitPur.addSourceUSimple(Omega.getActiveDofCount() + nbB1 + i, indB2 + i, indCommun2, "Commun2", 0, 0);
//        }
//        circuitPur.addSourceUSimple(Omega.getActiveDofCount() + nbB1 + nbB2, indCommun2, indCommun1, "SourceI", 1, 0);

        int N1 = Omega.getElementSet().getNbElement() + 1;
        int N2 = N1 + 1;
//        int iB1[] = PEEC_RLMPC_ConducEquiv.getIndicesBornes(Omega, (SurfaceRegion) ifl.getRegion(3));
//        int iB2[] = PEEC_RLMPC_ConducEquiv.getIndicesBornes(Omega, (SurfaceRegion) ifl.getRegion(2));
//
//        for (int i = 0; i < iB1.length; i++) {
//            circuitPur.addSourceUSimple(Omega.getActiveDofCount() + i, iB1[i], N1, "Commun1", 0, 0);
//        }
//        for (int i = 0; i < iB2.length; i++) {
//            circuitPur.addSourceUSimple(Omega.getActiveDofCount() + iB1.length + i, iB2[i], N2, "Commun2", 0, 0);
//        }
        circuitPur.addSourceISimple(solP.getNbLignes(), N2, N1, "Source I", 1, 0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

        // Resolution
        double ib[][] = solP.resolutionDirecte(f);
        Matrix res = new Matrix(2, Omega.getActiveDofCount());
        for (int i = 0; i < Omega.getActiveDofCount(); i++) {
            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
        }
        
        System.out.println("I= "+ib[0][2*solP.getNbLignes()]+" + j* "+ib[0][2*solP.getNbLignes()+1]);
        System.out.println("U= "+ib[1][2*solP.getNbLignes()]+" + j* "+ib[1][2*solP.getNbLignes()+1]);

//        double[] z = new double[]{ib[1][2 * nbRef], ib[1][2 * nbRef - 1]};
//        System.out.println("Z(" + f + ")= " + z[0] + " + i* " + z[1]);
//        System.out.println("|Z|= " + Math.hypot(z[0], z[1]));
//        System.out.println("I= " + ib[0][2 * nbRef] + " +i* " + ib[0][2 * nbRef - 1]);
        // On reorganise les dof   
        RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
        RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));

        ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_re.msh");
        exportJreal.addQuantity(Jreal, "Jreal");

        ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_im.msh");
        exportJimag.addQuantity(Jimag, "Jimag");

        ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
        exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_Jmod.msh");
        exportJreal.addQuantityExportMod(J, "Jmod");

        File out_file = new File("D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_Jmod.msh");
        File out_file1 = new File("D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_im.msh");
        File out_file2 = new File("D:/jsiau/_Backup_Sources/Resultats/LoopAntenna/RLMPC_f" + f + "_nDof" + Omega.getActiveDofCount() + "_re.msh");
        try {
            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);
        } catch (IOException ex) {
            Logger.getLogger(LoopAntenna.class.getName()).log(Level.SEVERE, null, ex);
        }
        return res;
    }

    private void ecrireMat(Matrix r1, String name) {
        Ecriture e1 = null;
        System.out.println("Ecriture de la matrix (" + r1.getRowCount() + " , " + r1.getColumnCount() + ")");
        try {
            e1 = new Ecriture(name);
            double tab[][] = new double[r1.getRowCount()][r1.getColumnCount()];
            for (int i = 0; i < r1.getRowCount(); i++) {
                for (int j = 0; j < r1.getColumnCount(); j++) {
                    tab[i][j] = r1.getElement(i, j);
                }
            }
            e1.ecrire(tab, ' ');
            e1.close();
        } catch (IOException ex) {
            Logger.getLogger(LoopAntenna.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
