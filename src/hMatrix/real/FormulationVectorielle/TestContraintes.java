/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.FormulationVectorielle;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.contraints.InactiveBorderNode;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.functionSpace.CurlSEdgeDeg1;
import g2elab.mipse.meshCore.functionSpace.GradNodalDeg1;
import g2elab.mipse.meshCore.functionSpace.Hgrad;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.CrossDGmulAlpha;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;

/**
 * TEST LASSEMBLAGE AVEC CONTRAINTE DE DIRICHLET NULLE
 *
 * @author jsiau
 */
public class TestContraintes {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File f = new File("");
        String path = f.getAbsolutePath();
        String file = path + "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/PLAQUE2000.DEC";

        ImportFlux ImF = new ImportFlux(file);
        Region reg = ImF.getRegion(0);
        ElementSurfSetHomogene mesh = (ElementSurfSetHomogene) reg.getElementSet();

        System.out.println("Nombre elements " + mesh.getNbElement());
        System.out.println("Nombre de noeuds " + mesh.getNbNoeud());

        double mu0 = 4 * Math.PI * Math.pow(10, -7);
        double sigma = 6 * Math.pow(10, 7);
        double frequence = 0.01;
        double epaisseur = 50e-3;

        System.out.println(" Conductivité électrique = " + sigma + " S/m");
        System.out.println(" Perméabilité = " + mu0 + " H/m");
        System.out.println(" Fréquence = " + frequence + " Hz");
        System.out.println(" ");

        // Constrainte de dirichlet null sur la frontiere
        InactiveBorderNode Dir0border = new InactiveBorderNode();

        // Espace nodal gradHgrad
        GradNodalDeg1 gradalpha = new GradNodalDeg1(mesh, Dir0border);
        // Espace Fonctionnel pour T
        CurlSEdgeDeg1 rotSalpha = new CurlSEdgeDeg1(gradalpha);

        // Espace fonction fonction de forme
        NodalDeg1 alpha = new NodalDeg1(mesh, Dir0border);
        // Espace fcontionel alpha.N
        Hgrad alphaN = alpha.createProjOnNormal(1);

        int ndof = alpha.getActiveDofCount();

        int renum[] = alpha.getRenum();
//        Node nFixd[] = new Node[mesh.getNbNoeud() - ndof];
//        System.out.println("nFixd.length= " + nFixd.length);
//        System.arraycopy(alpha.getNumberedNodesList(), ndof, nFixd, 0, nFixd.length);
//
//        Node[] N = Dir0border.getRegion(alphaN);
//        System.out.println("Nlength= " + N.length);
//
//        int cpt = 0;
//        for (int i = 0; i < N.length; i++) {
//            boolean b = true;
//            for (int j = 0; (j < nFixd.length) && b; j++) {
//                if( Arrays.equals(N[i].getCoord(), nFixd[j].getCoord()) ){
//                    System.out.println("("+i+" , "+j+")");
//                    b = false;
//                    cpt++;
//                }                    
//            }
//        }
//        System.out.println("cpt= "+cpt);


        System.out.println("Assemblage Plein");
        GalerkinIntegralFormulationFull IFp = new GalerkinIntegralFormulationFull(alphaN, rotSalpha, new CrossDGmulAlpha(2 * Math.PI * frequence * mu0 * epaisseur), new SelfElementFixedGauss(3, new Cancel()), 3);
        long deb = System.nanoTime();
        IFp.assembly();
        long fin = System.nanoTime();
        System.out.println("fin = " + (fin - deb) / 1e9 + "secs");

        System.out.println("Assemblage Hmatrix");
        double eps = 1e-7;
        int kmax = 50, nmin = 30, order = 6;
        GalerkinIntegralFormulationHCA IF = new GalerkinIntegralFormulationHCA(alphaN, rotSalpha, new CrossDGmulAlpha(2 * Math.PI * frequence * mu0 * epaisseur), new SelfElementFixedGauss(3, new Cancel()), 3, 3,
                eps, kmax, nmin, order);
        deb = System.nanoTime();
        IF.assembly();
        fin = System.nanoTime();
        System.out.println("fin = " + (fin - deb) / 1e9 + "secs");

        ColumnVector v = new ColumnVector(ndof);
        for (int i = 0; i < ndof; i++) {
            v.setElement(i, Math.random() * 1000);
        }

        System.out.println("CALCUL DU PMV Hmat");
        StorageHmatrix H = (StorageHmatrix) IF.getStore();
        ColumnVector vH = H.prod(v);

        System.out.println("CALCUL DU PMV Plein");
        Matrix M = ((StorageFull) IFp.getStore()).getMatrix();
        ColumnVector vP = new ColumnVector(ndof);
        vP.mul(M, v);

        vH.sub(vP);
        double errA = vH.norm();
        System.out.println("Erreur absolue= " + errA + "\t Erreur relative= " + errA / M.norm());

        /*
         CHECK THE ERROR BLOCK-WISE
         */
        Matrix M2 = new Matrix(ndof, ndof);
        int ind[] = H.getDof2Idx();
        for (int i = 0; i < M.getRowCount(); i++) {
            for (int j = 0; j < M.getColumnCount(); j++) {
                M2.setElement(i, j, M.getElement(ind[i], ind[j]));
            }
        }
        H.CheckError(M2, eps);
    }

}
