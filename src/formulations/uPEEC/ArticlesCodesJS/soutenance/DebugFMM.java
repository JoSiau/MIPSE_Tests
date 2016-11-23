/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.ArticlesCodesJS.soutenance;

import static formulations.uPEEC.ArticlesCodesJS.soutenance.Comparison_CompTechs.getSolver;
import g2elab.mipse.formulationInProgress.magnetodynamic.U_PEEC_DIELECTRIC.Cond_n_Dielec_Volumic.PEEC_RLMPC_Volume;
import g2elab.mipse.mipseCore.matrixCompression.Compression;
import g2elab.mipse.numericalTools.matrix.complex.AbstractMatrixComplex;
import g2elab.mipse.numericalTools.matrix.real.AbstractMatrixReal;
import g2elab.mipse.numericalTools.vector.full.VectorFullComplex;
import g2elab.mipse.numericalTools.vector.full.VectorFullReal;

/**
 *
 * @author jsiau
 */
public class DebugFMM {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        double f = 1e6;
        PEEC_RLMPC_Volume solFMM = getSolver(Compression.HCA, f);
        PEEC_RLMPC_Volume solFULL = getSolver(Compression.No, f);

        AbstractMatrixReal[][] lfmm = solFMM.getMatL();
        AbstractMatrixReal[][] l = solFULL.getMatL();
        //
        //
        System.out.println("Erreur par bloc");
        for (int i = 0; i < solFULL.getFinalMatrix().getNbBlocRows(); i++) {
            for (int j = 0; j < solFULL.getFinalMatrix().getNbBlocColumns(); j++) {
                if (solFULL.getFinalMatrix().getBlock(i, j) != null) {
                    System.out.println("Erreur bloc[" + i + "][" + j + "] = " + getErr(solFULL.getFinalMatrix().getBlock(i, j), solFMM.getFinalMatrix().getBlock(i, j)));
                }
            }
        }
        //
        //
        System.out.println("Erreur L reels");
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                System.out.println("Erreur L[" + i + "][" + j + "] = " + getErr(l[i][j], lfmm[i][j]));
            }
        }
        //
        //
        System.out.println("Erreur L complexe");
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                System.out.println("Erreur L[" + i + "][" + j + "] = " + getErr(solFULL.getFinalMatrix().getBlock(i, j), solFMM.getFinalMatrix().getBlock(i, j)));
            }
        }
        System.out.println("Erreur global = " + getErr(solFULL.getFinalMatrix(), solFMM.getFinalMatrix()));
    }

    public static double getErr(AbstractMatrixReal m, AbstractMatrixReal fmm) {
        VectorFullReal x = VectorFullReal.generateRandom(m.getColumns());
        VectorFullReal b = new VectorFullReal(m.product(x.getValues(false), new double[m.getRows()]));
        VectorFullReal bfmm = new VectorFullReal(fmm.product(x.getValues(false), new double[m.getRows()]));
        return bfmm.sub(b, null).norm2() / b.norm2();
    }

    public static double getErr(AbstractMatrixComplex ref, AbstractMatrixComplex fmm) {
        VectorFullComplex x = VectorFullComplex.generateRandom(ref.getColumns());
        VectorFullComplex b = new VectorFullComplex(ref.product(x.getValues(false), new double[2 * ref.getRows()]));
        VectorFullComplex bfmm = new VectorFullComplex(fmm.product(x.getValues(false), new double[2 * ref.getRows()]));
        return bfmm.sub(b, null).norm2() / b.norm2();
    }

}
