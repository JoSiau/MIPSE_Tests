/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.impedanceCurves;

import g2elab.mipse.circuit.solverCircuitComplex.electricalBlock.BlocElectriqueBasique;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv;
import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.region.SurfaceRegion;
import g2elab.mipse.meshCore.region.VolumeRegion;
import g2elab.mipse.tools.files.Ecriture;
import g2elab.mipse.tools.multiThreads.GestionnaireTaches;

import java.io.File;
import java.io.IOException;

import static g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_ConducEquiv.eps0;

/**
 *
 * @author jsiau
 */
public class LoopAntenna_WithDielec_Zf {

    public static void main(String[] args) throws IOException {
        
        int nbF = 50;
        double f_dep = 1e9;
        double f_fin = 5e9;
        double pas = (f_fin - f_dep) / (nbF - 1);
        System.out.println("pas= " + pas);
        double f[] = new double[nbF];
        for (int i = 0; i < nbF; i++) {
            f[i] = f_dep + i * pas;
        }
        
        double z[] = new double[nbF];
        Ecriture save = new Ecriture("D:/Zf_LAWD_rlmpc1.txt");
        
        for (int i = 0; i < nbF; i++) {
            System.out.println(i+" / "+nbF+"\t f= "+f[i]);
            z[i] = doIt(f[i]);
            for (int j = 0; j <= i; j++) {
                System.out.println(f[j]+" "+z[j]);
            }
            save.ecrire(f[i] + " " + z[i] + "\n");
        }
        save.close();
        GestionnaireTaches.getGestionnaireTaches().stop();
    }
    
    public static double doIt(double f){
        File fil = new File("");
        String path = fil.getAbsolutePath();
        /*
         Nombre de regions importees : 4
         Region 0, Nom : LOOP_ANTENNA, type : 3, 1305 elements
         Region 1, Nom : DIELECTRIC, type : 3, 2354 elements
         Region 2, Nom : POSITIVE, type : 2, 4 elements
         Region 3, Nom : NEGATIVE, type : 2, 4 elements
         */
//        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/LOOPANTENNA_SMALL.DEC";
        String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/U_PEEC_DIELECTRIC/LOOPANTENNA_DIELEC_VOL4.DEC";
        
        ImportFlux Imp = new ImportFlux(file);
        int nDof[] = new int[4];
        for (int i = 0; i < nDof.length; i++) {
            nDof[i] = Imp.getRegion(i).getElementSet().getNbElement();
        }
        VolumeRegion rV[] = new VolumeRegion[2];
        rV[0] = (VolumeRegion) Imp.getRegion(0);// Region conductrice
        rV[1] = (VolumeRegion) Imp.getRegion(1);// Region dielectrique

        /*
         CONSTRUCTION DU SUPPORT DES CONDUCTIVITES EQUIVALENTES
         */
        double Sigma[][] = new double[2][nDof[0] + nDof[1]];// Size = [2][nbElmttotal]
        // Conducteur
        double Re = 1/ 1.72e-8;
        double Im = eps0;
        for (int i = 0; i < nDof[0]; i++) {
            Sigma[0][i] = Re;// Conductivite du cuivre
            Sigma[1][i] = Im;// Epsilon0
        }
        // Dielectrique
        for (int i = nDof[0]; i < nDof[0] + nDof[1]; i++) {
            Sigma[0][i] = 0;// 0 pour le dielec
            Sigma[1][i] = 5*eps0;// Espilon
        }
        // Create the solver
        PEEC_RLMPC_ConducEquiv solP = new PEEC_RLMPC_ConducEquiv(rV, Sigma, (SurfaceRegion) Imp.getRegion(2), (SurfaceRegion) Imp.getRegion(3), true);
        solP.setPtsDeGaussInductifs(15, 4);
        solP.setPtsDeGaussCapacitifs(4, 4);
        solP.setVerbose(false);

        BlocElectriqueBasique circuitPur = new BlocElectriqueBasique();

        FaceDeg1 Omega = solP.getFunctionSpace();
        System.out.println("Activedofcount= " + Omega.getActiveDofCount());

        /*
        CREATION DES BORNES
        */

        int nbBranches = solP.getNbLignes();
        int C1 = Omega.getNbElement() + 1;
        int C2 = C1 + 1;

        circuitPur.addSourceISimple(nbBranches, C1, C2, "Source I", 1, 0);

        circuitPur.finSaisie();
        solP.setCircuitElectrique(circuitPur);

        // Resolution
        double ib[][] = solP.resolutionDirecte(f);
//        Matrix res = new Matrix(2, Omega.getActiveDofCount());
//        for (int i = 0; i < Omega.getActiveDofCount(); i++) {
//            res.setElement(0, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i]);
//            res.setElement(1, solP.getSolveurCircuit().getNumIdBranche(i), ib[0][2 * i + 1]);
//        }

        System.out.println("I= "+ib[0][2*nbBranches]+" + j* "+ib[0][2*nbBranches+1]);
        System.out.println("U= "+ib[1][2*nbBranches]+" + j* "+ib[1][2*nbBranches+1]);
        return Math.hypot(ib[1][2 * nbBranches], ib[1][2 * nbBranches + 1]);
        
//        // On reorganise les dof   
//        RealFaceQuantity Jreal = new RealFaceQuantity(Omega, res.row(0));
//        RealFaceQuantity Jimag = new RealFaceQuantity(Omega, res.row(1));
//
//        ExportGmshHdiv exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/WireRE_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
//        exportJreal.addQuantity(Jreal, "Jreal");
//
//        ExportGmshHdiv exportJimag = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/WireIM_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
//        exportJimag.addQuantity(Jimag, "Jimag");
//
//        ComplexFaceQuantity J = new ComplexFaceQuantity(Omega, res);
//        exportJreal = new ExportGmshHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/WireJmod" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
//        exportJreal.addQuantityExportMod(J, "Jmod");
//        
//        ExportVtkHdiv tmp = new ExportVtkHdiv(Omega, "D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/Jmod.vtk");
//        tmp.addQuantity(J, "J");
//
//        File out_file = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/WireJmod" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
//        File out_file1 = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/WireRE_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
//        File out_file2 = new File("D:/jsiau/_Backup_Sources/Resultats/Dielectric/loopAntenna/WireIM_f" + f + "_nDof" + Omega.getActiveDofCount() + ".msh");
//        try {
//            Runtime.getRuntime().exec("C:/gmsh-2.8.4-Windows/gmsh.exe " + out_file + " " + out_file1 + " " + out_file2);
//
//        } catch (IOException ex) {
//            Logger.getLogger(PEEC_RLMPC_ConducEquiv.class.getName()).log(Level.SEVERE, null, ex);
//        }
//
//        double x[] = new double[2*solP.getNbLignes()];
//        for (int i = 0; i < x.length; i++) {
//            x[i] = Math.random();
//        }
//        
//        double zb[][] = new double[solP.getNbLignes()][2*solP.getNbLignes()];
//        solP.getZbFull(zb, 0, solP.getNbLignes()-1, 0, solP.getNbLignes()-1, f);
//        Basic2D Z = new Basic2D(zb);
//        double vf[] = Z.product(x, new double[2*solP.getNbLignes()]);
//        ColumnVector rvf = new ColumnVector(vf);
//        
//        double vp[] = solP.produit(x, new double[2*solP.getNbLignes()], f);
//        ColumnVector rvp = new ColumnVector(vp);
//        
//        rvp.sub(rvf);
//        System.out.println("Erreur abs= "+rvp.norm());
//        System.out.println("erreur relative= "+rvp.norm()/rvf.norm());
//        GestionnaireTaches.getGestionnaireTaches().stop();

    }

}
