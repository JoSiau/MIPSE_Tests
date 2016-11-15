/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.complex.imp;


import g2elab.mipse.numericalTools.matrix.complex.dense.basic2D.DecompositionLUBasic2D;

import java.util.Arrays;

/**
 * @author jsiau
 */
public class TestLU {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        int n = 10;
        double m[][] = new double[n][2 * n];
        for (double[] m1 : m) {
            for (int j = 0; j < m1.length; j++) {
                m1[j] = Math.random() * 10;
            }
        }
        DecompositionLUBasic2D lu = new DecompositionLUBasic2D(m, 1e-20, false, true, true);
        System.out.println("Pr = " + Arrays.toString(lu.getPermLigne()));
        System.out.println("Pc = " + Arrays.toString(lu.getPermColonne()));


    }

}
