/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.PreconditionnerLU;

import got.matrix.ColumnVector;
import got.matrix.LuDecomposition;
import got.matrix.Matrix;
import got.matrix.RowVector;

/**
 * @author siau
 */
public class TestPermutationsLU {

    /**
     * @param args
     */
    public static void main(String[] args) {
        double deb, fin;
        int n = 650;
        Matrix A = new Matrix(n, n);
        A.setZero();
        for (int i = 0; i < A.getRowCount(); i++)
            for (int j = 0; j < A.getColumnCount(); j++)
                A.setElement(i, j, Math.random() * 100);
//        System.out.println("A\n"+A);

        LuDecomposition dLU = new LuDecomposition(A, true);
        Matrix L = new Matrix(n, n);
        dLU.getG(L);
        Matrix U = new Matrix(n, n);
        dLU.getH(U);
        //*
        L = L.copy();
        U = U.copy();
        //*/
        Matrix mlu = new Matrix(n, n);
        mlu.mul(L, U);
        mlu.sub(A);
        System.out.println("Erreur exacte LxU= " + mlu.norm());

        deb = System.nanoTime();
        int pr[] = new int[n];
        dLU.getRowPermutations(pr);
        Matrix permr = new Matrix(n, n);
        for (int i = 0; i < permr.getRowCount(); i++) {
            permr.setElement(pr[i], i, 1);
        }
        Matrix tmpp = new Matrix(n, n);
        tmpp.mul(permr, L);
        fin = System.nanoTime();
        System.out.println("Time P1= " + (fin - deb) * 1e-9);

        // P^{-1}.L        
        deb = System.nanoTime();
        Matrix tmp0 = new Matrix(n, n);
        tmp0 = L.copy();
        tmp0.permuteRows(pr, new int[n]);
        fin = System.nanoTime();
        System.out.println("Time P2= " + (fin - deb) * 1e-9);
        // tmp0.mul(permr, L);


        deb = System.nanoTime();
        Matrix Am = new Matrix(n, n);
        Matrix tmp = L;
        //a
        RowVector vr = new RowVector(Am.getColumnCount());
        for (int i = 0; i < Am.getRowCount(); i++) {
            tmp.getRow(i, vr);
            Am.setRow(pr[i], vr);
        }
        fin = System.nanoTime();
        System.out.println("Time P3= " + (fin - deb) * 1e-9);
        
        /*
        COLUMN PERM
        */
        int pc[] = new int[n];
        dLU.getColumnPermutations(pc);
        Matrix permc = new Matrix(n, n);
        permc.setIdentity();
        permc.permuteColumns(pc, new int[n]);
        // U.Q
        Matrix tmp1 = new Matrix(n, n);
        tmp1 = U.copy();
        tmp1.permuteColumns(pc, new int[n]);


        Matrix LU = new Matrix(n, n);
        LU.mul(tmp0, U);

        LU.sub(A);
        System.out.println("Erreur rel P-1xLxU= " + LU.norm() / A.norm());

        LU.mul(tmp0, tmp1);
        LU.sub(A);
        System.out.println("Erreur rel P-1xLxUxQ= " + LU.norm() / A.norm());


        Matrix L1 = new Matrix(n, n);
        RowVector r = new RowVector(n);
        for (int i = 0; i < n; i++) {
            L.getRow(i, r);
            L1.setRow(pr[i], r);
        }
        LU.mul(L1, U);
        LU.sub(A);
        System.out.println("Erreur rel P-1xLxU= " + LU.norm() / A.norm());


        ColumnVector b = new ColumnVector(n);
        for (int i = 0; i < n; i++)
            b.setElement(i, Math.random() * 1000);


        ColumnVector xT = testSolve(dLU, n, b);

        ColumnVector xLU = new ColumnVector(n);
        dLU.solve(b, xLU);

        xT.sub(xLU);
        System.out.println("Error abs= " + xT.norm());
        System.out.println("Error rel= " + xT.norm() / xLU.norm());
    }

    /**
     * @param dLU
     * @param n
     * @param b
     * @return
     */
    public static ColumnVector testSolve(LuDecomposition dLU, int n, ColumnVector b) {

        Matrix L = new Matrix(n, n);
        dLU.getG(L);
        Matrix U = new Matrix(n, n);
        dLU.getH(U);

        int pr[] = new int[n];
        dLU.getRowPermutations(pr);
        int pc[] = new int[n];
        dLU.getColumnPermutations(pc);
        
        
        /*
         * 1ere etape
         */
        ColumnVector c = new ColumnVector(b.getRowCount());
        for (int i = 0; i < n; i++) {
//                c.setElement(this.Pr[L.getRowStart()+i], b.getElement(i));
            c.setElement(i, b.getElement(pr[i]));
        }
        ColumnVector y = new ColumnVector(n);
        LaDescente(c, L, y);  
        
        /*
         * 2eme etape
         */
        ColumnVector z = new ColumnVector(n);
        LaRemontee(y, U, z);

        ColumnVector x = new ColumnVector(n);
        for (int i = 0; i < n; i++)
            x.setElement(pc[i], z.getElement(i));


        return x;
    }

    /**
     * @param A
     * @param L
     * @param M
     */
    protected static void LaDescente(Matrix A, Matrix L, Matrix M) {
        // Set the 1st row
//       double pivot = L.getElement(0, 0);// Normally pivot==1
        for (int i = 0; i < M.getColumnCount(); i++)
            M.setElement(0, i, A.getElement(0, i));

        // Do all the work
        double res;
        for (int i = 1; i < M.getRowCount(); i++) {
//           pivot = L.getElement(i, i);
            for (int j = 0; j < M.getColumnCount(); j++) {
                res = A.getElement(i, j);
                for (int k = 0; k < i; k++)
                    res -= L.getElement(i, k) * M.getElement(k, j);

//               M.setElement(i, j, res/pivot);    
                M.setElement(i, j, res);
            }
        }
    }

    /**
     * @param b
     * @param U
     * @param x
     * @throws IllegalArgumentException
     */
    protected static void LaRemontee(ColumnVector b, Matrix U, ColumnVector x) throws IllegalArgumentException {
        if (U.getColumnCount() != x.getRowCount() || U.getRowCount() != b.getRowCount())
            throw new IllegalArgumentException("Error in LaRemontee(): Dimensions Dismatches");

//        int n = x.getRowCount()-1;
//        x.setElement(n, b.getElement(n)/U.getElement(n, n));

        for (int i = x.getRowCount() - 1; i >= 0; i--) {
//            double pivot = U.getElement(i, i);
            double res = b.getElement(i);

            for (int j = i + 1; j < x.getRowCount(); j++)
                res -= U.getElement(i, j) * x.getElement(j);

            x.setElement(i, res / U.getElement(i, i));
        }
    }

}
