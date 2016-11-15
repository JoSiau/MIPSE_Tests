/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real;

import g2elab.mipse.mipseCore.matrixCompression.aca.ACACompatible;
import got.matrix.ColumnVector;
import got.matrix.Matrix;
import got.matrix.RowVector;

/**
 * @author siau
 */
public class CompColACA implements ACACompatible {

    Matrix X, Y;

    /**
     * @param A
     * @param B
     */
    public CompColACA(Matrix A, Matrix B) {
        X = A;
        Y = B;
    }

    @Override
    public int getACAColumnCount() {
        return Y.getRowCount();
    }

    @Override
    public int getACARowCount() {
        return Y.getRowCount();
    }

    @Override
    public ColumnVector getACAColumn(int j) {
        return CompCol(X, Y, j, 0);
    }

    @Override
    public RowVector getACARow(int i) {
        return CompCol(X, Y, i, 1).transpose();
    }

    private ColumnVector CompCol(Matrix Xi, Matrix Yi, int ind, int rorc) {
        double alpha = 0.1;
        ColumnVector v = null;
        Matrix S;

        if (rorc == 0) {
            RowVector Coli = new RowVector(Yi.getColumnCount());
            Yi.getRow(ind, Coli);

            S = new Matrix(Xi.getColumnCount(), Xi.getRowCount());
            // S.sub(Xi.transpose(), kronProdVect(one,Coli.transpose()));
            S.sub(Xi.transpose(), this.MakeKronRes(Xi.getRowCount(), Coli.transpose()));

            // Compute the : "S.^2"
            for (int i = 0; i < S.getRowCount(); i++) {
                for (int j = 0; j < S.getColumnCount(); j++) {
                    S.setElement(i, j, S.getElement(i, j) * S.getElement(i, j));
                }
            }
            // Additionne les lignes de S en stockant les résultats dans un Columnvector !
            v = new ColumnVector(S.getColumnCount());

            double tmp;
            for (int i = 0; i < S.getColumnCount(); i++) {
                // v.setElement(i, 0);
                tmp = 0;
                for (int j = 0; j < S.getRowCount(); j++) {
                    //v.setElement(i,v.getElement(i)+S.getElement(j, i));
                    tmp += S.getElement(j, i);
                }

                v.setElement(i, 1 / (alpha + tmp));
            }
        } else {
            RowVector Coli = new RowVector(Xi.getColumnCount());
            Xi.getRow(ind, Coli);

            S = new Matrix(Yi.getColumnCount(), Yi.getRowCount());
            // S.sub(Yi.transpose(), kronProdVect(one,Coli.transpose()));
            S.sub(Yi.transpose(), this.MakeKronRes(Yi.getRowCount(), Coli.transpose()));

            // Compute the : "S.^2"
            for (int i = 0; i < S.getRowCount(); i++) {
                for (int j = 0; j < S.getColumnCount(); j++) {
                    S.setElement(i, j, S.getElement(i, j) * S.getElement(i, j));
                }
            }
            // Additionne les lignes de S en stockant les résultats dans un Columnvector !
            v = new ColumnVector(S.getColumnCount());

            double tmp;
            for (int i = 0; i < S.getColumnCount(); i++) {
                // v.setElement(i, 0);
                tmp = 0;
                for (int j = 0; j < S.getRowCount(); j++) {
                    //v.setElement(i,v.getElement(i)+S.getElement(j, i));
                    tmp += S.getElement(j, i);
                }

                v.setElement(i, 1 / (alpha + tmp));
            }
        }
        return v;
    }

    private Matrix kronProdVect(RowVector V, ColumnVector W) {
        Matrix K = new Matrix(W.getRowCount(), V.getColumnCount());
        for (int i = 0; i < K.getRowCount(); i++)
            for (int j = 0; j < K.getColumnCount(); j++)
                K.setElement(i, j, V.getElement(j) * W.getElement(i));

        return K;
    }

    private Matrix MakeKronRes(int n, ColumnVector W) {
        Matrix K = new Matrix(W.getRowCount(), n);
        for (int i = 0; i < n; i++)
            K.setColumn(i, W);

        return K;
    }

}
