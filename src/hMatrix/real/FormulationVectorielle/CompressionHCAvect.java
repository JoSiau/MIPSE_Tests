/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.FormulationVectorielle;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ImportGmshMeshRegion;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.functionSpace.CurlSEdgeDeg1;
import g2elab.mipse.meshCore.functionSpace.GradNodalDeg1;
import g2elab.mipse.meshCore.functionSpace.Hgrad;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.CrossDGmulAlpha;
import g2elab.mipse.mipseCore.matrixCompression.hMatrix.Real.TruncationControl;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;

import java.io.IOException;

/**
 * FICHIERS POUR OBTENIR LES COMPRESSIONS DE LA HCA EN FORMULATION VECTORIELLE
 *
 * @author jsiau
 */
public class CompressionHCAvect {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        int nbGauss = 7;
        double eps = 1e-5;
        TruncationControl tol = new TruncationControl("rel", eps);
        int kmax = 50, nmin = 30, order = (int) -Math.log10(eps) - 1;
        System.out.println("ordre= " + order);

        String file = new java.io.File(".").getCanonicalPath();
        file += "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/Sphere_gmsh/sphere_";
        //*        
        int nbElm[] = new int[]{6481, 12529, 34449, 62713, 132241};
         /*/
        String file = "D:/Meshs/Cube/CUBE_";
        // 1568 , 6152 , 18928 , 38444 , 75320 , 153312
        int nbElm[] = new int[]{1568, 6152, 18928, 38444, 75320, 153312};
        //*/

        int nTest = nbElm.length;

        double stor[][] = new double[nTest][2];
        double time[][] = new double[nTest][2];
        double storAgglo[][] = new double[nTest][2];
        double timeAgglo[][] = new double[nTest][2];


        GalerkinIntegralFormulationHCA IF2;
        ImportFlux ImF;
        ImportGmshMeshRegion mesh1;
        ElementSurfSetHomogene mesh;
        StorageHmatrix H;
        long deb;
        for (int i = 0; i < nTest; i++) {
            /*
            ImF = new ImportFlux( file + nbElm[i] + ".DEC");
            mesh = (ElementSurfSetHomogene) ImF.getRegions(0).getElementSet();
            /*/
            mesh1 = new ImportGmshMeshRegion(0, file + nbElm[i] + ".msh");
            mesh = (ElementSurfSetHomogene) mesh1.createHomogeneESet();
            //*/
            System.out.println("Nombre elements " + mesh.getNbElement());
            System.out.println("Nombre de noeuds " + mesh.getNbNoeud());

            /*
             CONSTRUCTION FUNCTION SPACE
             */
            NodalDeg1 alpha = new NodalDeg1(mesh);
            Hgrad alphaN = alpha.createProjOnNormal(1);
            GradNodalDeg1 gradalpha = new GradNodalDeg1(mesh);
            CurlSEdgeDeg1 rotSalpha = new CurlSEdgeDeg1(gradalpha);

            // Test sans recompression
            IF2 = new GalerkinIntegralFormulationHCA(alphaN, rotSalpha, new CrossDGmulAlpha(/*frequence * mu0 * epaisseur*/1), new SelfElementFixedGauss(3, new Cancel()), nbGauss, nbGauss,
                    eps, kmax, nmin, order, 2.0, false);
            deb = System.currentTimeMillis();
            IF2.assembly();
            time[i][0] = (System.currentTimeMillis() - deb) * 1e-3;
            H = (StorageHmatrix) IF2.getStore();
            stor[i][0] = H.getStorageSizePerDofs();

            // Agglomeration
            deb = System.currentTimeMillis();
            H.Agglomerate(tol);
            timeAgglo[i][0] = (System.currentTimeMillis() - deb) * 1e-3;
            storAgglo[i][0] = H.getStorageSizePerDofs();


            // Test avec recompression
            IF2 = new GalerkinIntegralFormulationHCA(alphaN, rotSalpha, new CrossDGmulAlpha(/*frequence * mu0 * epaisseur*/1), new SelfElementFixedGauss(3, new Cancel()), nbGauss, nbGauss,
                    eps, kmax, nmin, order, 2.0, true);
            deb = System.currentTimeMillis();
            IF2.assembly();
            time[i][1] = (System.currentTimeMillis() - deb) * 1e-3;
            H = (StorageHmatrix) IF2.getStore();
            stor[i][1] = H.getStorageSizePerDofs();

            // Agglomeration
            deb = System.currentTimeMillis();
            H.Agglomerate(tol);
            timeAgglo[i][1] = (System.currentTimeMillis() - deb) * 1e-3;
            storAgglo[i][1] = H.getStorageSizePerDofs();
        }

        for (int i = 0; i < nTest; i++) {
            System.out.println(time[i][0] + " \t , " + stor[i][0] + " \t , " + timeAgglo[i][0] + " \t , " + storAgglo[i][0] + " \t , " + time[i][1] + " \t , " + stor[i][1] + " \t , " + timeAgglo[i][1] + " \t , " + storAgglo[i][1]);
        }
    }

}
/*
170.374 	 , 78.55691712350854 	 , 88.474 	 , 8.010641728474685 	 , 246.485 	 , 79.53789100290228 	 , 14.214 	 , 7.912931312479845
361.577 	 , 94.88239169823753 	 , 148.19 	 , 8.722780431559874 	 , 508.837 	 , 95.96244440784055 	 , 16.445 	 , 8.705485093065393
1077.292 	 , 128.39201419698315 	 , 475.308 	 , 10.762259686483288 	 , 1551.489 	 , 129.46909198461992 	 , 46.477 	 , 10.755042886719906
2152.262 	 , 146.14245060311094 	 , 929.669 	 , 12.339908805743297 	 , 3089.485 	 , 147.23303043042395 	 , 86.885 	 , 12.341655078743978
5739.05 	 , 168.10285696834583 	 , 9020.579 	 , 14.305257371466963 	 , 7721.034 	 , 169.11113316740216 	 , 3313.878 	 , 14.31400693246194
*/