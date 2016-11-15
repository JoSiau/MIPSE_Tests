/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Magnetodynamic;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.elements.Node;
import g2elab.mipse.meshCore.functionSpace.CurlSEdgeDeg1;
import g2elab.mipse.meshCore.functionSpace.GradNodalDeg1;
import g2elab.mipse.meshCore.functionSpace.Hgrad;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.NoCorrection;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.integrationstrategies.FixedGaussSource;
import g2elab.mipse.mipseCore.integralIntegration.kernel.CrossDGmulAlpha;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.integralIntegration.kernel.NegDotDG;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.BinaryTree.BlockCluster;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.BinaryTree.BlockClusterTree;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.HCA.CollocReduit;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.StorageBlock.Admissible.RkMatrix;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.clustering.Cluster;
import g2elab.mipse.mipseCore.numericalMethods.CollocationIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.File;
import java.io.IOException;

/**
 * @author jsiau
 */
public class tmp {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {

        double mu0 = 4 * Math.PI * Math.pow(10, -7);
        double sigma = 6 * Math.pow(10, 7);
        double frequence = 0.01;
        double epaisseur = 50e-6;

        File f = new File("");
        String path = f.getAbsolutePath();
        String file = path + "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/PLAQUE2000.DEC";
        //String file = path + "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere/SPHERE_1884.DEC";
        ImportFlux ImF = new ImportFlux(file);
        Region reg = ImF.getRegion(0);
        ElementSurfSetHomogene mesh = (ElementSurfSetHomogene) reg.getElementSet();

        System.out.println("Nombre elements " + mesh.getNbElement());
        System.out.println("Nombre de noeuds " + mesh.getNbNoeud());

        NodalDeg1 alpha = new NodalDeg1(mesh);//, Dir0border);
        Hgrad alphaN = alpha.createProjOnNormal(1);

        GradNodalDeg1 gradalpha = new GradNodalDeg1(mesh);//, Dir0border);
        CurlSEdgeDeg1 rotSalpha = new CurlSEdgeDeg1(gradalpha);

        double eps = 1e-5;
        int kmax = 50, nmin = 30, order = 4, nbGauss = 7;
        GalerkinIntegralFormulationHCA IFcross = new GalerkinIntegralFormulationHCA(alphaN, rotSalpha, new CrossDGmulAlpha(1), new SelfElementFixedGauss(3, new Cancel()), nbGauss, nbGauss,
                eps, kmax, nmin, order);

        GalerkinIntegralFormulationHCA IFdot = new GalerkinIntegralFormulationHCA(alpha, gradalpha, new NegDotDG(), new SelfElementFixedGauss(nbGauss, new Cancel()), nbGauss, nbGauss,
                eps, kmax, nmin, order);

        BlockClusterTree bct = new BlockClusterTree(alpha, 3, nmin, 2.0);
// 377 33 1033 34
        BlockCluster bc = bct.getBC(377, 1033);
        int ind[] = bct.getDof2IdxSrc();

        Cluster row = bc.getRowCluster();
        Cluster col = bc.getColCluster();

        RkMatrix Rcross = IFcross.getCompressedBlock(row, col);
        RkMatrix Rdot = IFdot.getCompressedBlock(row, col);
        System.out.println("Rcross.rank= " + Rcross.getRank() + "\t Rdot.rank= " + Rdot.getRank());

        Matrix tmp = Rcross.getFullMatrix();
        Matrix M = Rdot.getFullMatrix();
        tmp.sub(M);
        double tmpNorm = tmp.norm(), mNorm = M.norm();
        System.out.println("Erreur relative= " + tmpNorm / mNorm + "\t Erreur absolue= " + tmpNorm + "\t M.norm= " + mNorm);

        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        System.out.println("++++++++++++++++++++++++++++++++++++++++++++++++++++++");

        Matrix Xtarget = new Matrix(3, mesh.getNbNoeud());
        Node[] N = mesh.getLocalNodes();// rotSalpha.getNumberedNodesList();
        for (int i = 0; i < N.length; i++) {
            Xtarget.setColumn(i, (new ColumnVector(N[i].getCoord())).copy());
        }

        mesh.computeVois();

        Matrix Mr = new Matrix(mesh.getNbNoeud(), mesh.getNbNoeud());
        CollocReduit Cr = new CollocReduit(alphaN, new MultGvect());
        Mr = Cr.AssemblyGeneral(Xtarget, nbGauss, Mr, null);

        System.out.println("Mr=\n" + Mr.submatrix(0, 0, 15, 45));

        Matrix Mc = new Matrix(mesh.getNbNoeud(), mesh.getNbNoeud());
        CollocationIntegralFormulation Cc = new CollocationIntegralFormulation(alphaN, new MultGvect(), new NoCorrection(new FixedGaussSource(nbGauss)));
        Cc.matAssembly(rotSalpha, Mc);
        Mc = ((StorageFull) Cc.getStore()).getMatrix();

        Mr.sub(Mc);
        System.out.println("Erreur relative des collocations= " + Mr.norm() / Mc.norm() + "\t erreur abs= " + Mr.norm());

        if (mesh.getLocalNodes() != mesh.getGlobalNodes()) {
            System.out.println("Local != Global !");
        }

        if (rotSalpha.getNumberedNodesList() != alphaN.getNumberedNodesList()) {
            System.out.println("rotAlpha != alphaN !!!");
        }

        if (rotSalpha.getNumberedNodesList() != mesh.getLocalNodes()) {
            System.out.println("rotAlpha != Local !!!!!");
        }

        if (rotSalpha.getElementSet() != alphaN.getElementSet()) {
            System.err.println("WTF");
        }

    }

}
