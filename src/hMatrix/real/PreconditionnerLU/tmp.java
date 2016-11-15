/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.PreconditionnerLU;

import g2elab.mipse.numericalTools.matrix.real.dense.basic2D.DecompositionLUBasic2D;
import g2elab.mipse.tools.files.Lecture;
import got.matrix.LuDecomposition;
import got.matrix.Matrix;
import got.matrix.RankDeficientException;

import java.io.IOException;

/**
 * @author jsiau
 */
public class tmp {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {

        Lecture fic = new Lecture("D:/A.out");

        int n = fic.getNbLignes();
        double m[] = new double[n * n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(fic.lireLigneDoubleVecteur(' '), 0, m, i * n, n);
        }

        Matrix M = new Matrix(n, n, m);
        System.out.println("M:\n" + M);

        System.out.println("");
        LuDecomposition LUfull = new LuDecomposition();
        ;
        LUfull.setEpsilon(1e-120);
        try {
            LUfull.decompose(M, false);
        } catch (RankDeficientException ex) {
            System.out.println(ex);
            System.out.println("Pivot simple ne fonctionne pas !");
        }

        int pr[] = new int[n];
        LUfull.getRowPermutations(pr);
        int pc[] = new int[n];
        LUfull.getColumnPermutations(pc);

        /*
         * Assemble L and U into one matrice
         */
        Matrix tmpL = new Matrix(n, n);
        LUfull.getG(tmpL);
        Matrix tmpU = new Matrix(n, n);
        LUfull.getH(tmpU);


        Matrix lu = new Matrix(n, n);
        lu.mul(tmpL, tmpU);
        lu.permuteRows(pr, new int[n]);
        lu.sub(M);
        System.out.println("Erreur relative= " + lu.norm() / M.norm());
        
        /*
        
        */
        try {
            LUfull.decompose(M, true);
        } catch (RankDeficientException ex) {
            System.out.println(ex);
            System.out.println("Pivot complet ne fonctionne pas !");
        }
        pr = new int[n];
        LUfull.getRowPermutations(pr);
        pc = new int[n];
        LUfull.getColumnPermutations(pc);

        /*
         * Assemble L and U into one matrice
         */
        tmpL = new Matrix(n, n);
        LUfull.getG(tmpL);
        tmpU = new Matrix(n, n);
        LUfull.getH(tmpU);


        lu = new Matrix(n, n);
        lu.mul(tmpL, tmpU);
        lu.permuteRows(pr, new int[n]);
        lu.permuteColumns(pc, new int[n]);
        lu.sub(M);
        System.out.println("Erreur relative= " + lu.norm() / M.norm());


        /*
        
         */
        double A[][] = new double[n][n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(m, n * i, A[i], 0, n);
        }

        DecompositionLUBasic2D lu2d = new DecompositionLUBasic2D(A, 1e-120, true, true, false);
        double L2d[][] = lu2d.getL();
        double U2d[][] = lu2d.getU();
        int pL[] = lu2d.getPermLigne();

        double L1[] = new double[n * n];
        double U1[] = new double[n * n];
        for (int i = 0; i < n; i++) {
            System.arraycopy(L2d[i], 0, L1, n * i, n);
            System.arraycopy(U2d[i], 0, U1, n * i, n);
        }

        Matrix L = new Matrix(n, n, L1);
        Matrix U = new Matrix(n, n, U1);
        Matrix P = new Matrix(n, n);
        for (int i = 0; i < n; i++) {
            P.setElement(i, i, pL[i]);
        }
        Matrix LU = new Matrix(n, n);
        LU.mul(L, U);
        Matrix PA = new Matrix(n, n);
        PA.mul(P, M);

        LU.sub(PA);
        System.out.println("");

    }

}
