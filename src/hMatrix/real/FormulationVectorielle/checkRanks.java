/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.FormulationVectorielle;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.functionSpace.CurlSEdgeDeg1;
import g2elab.mipse.meshCore.functionSpace.GradNodalDeg1;
import g2elab.mipse.meshCore.functionSpace.Hgrad;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.CrossDGmulAlpha;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;

import java.io.File;

/**
 * Compare les rangs des Hmatrices avec et sans recompression.
 *
 * @author jsiau
 */
public class checkRanks {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {

        File f = new File("");
        String path = f.getAbsolutePath();
        String file = path + "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/";
        //*
        //        ImportFlux ImF = new ImportFlux(file+"PLAQUE2000.DEC");
        ImportFlux ImF = new ImportFlux(file + "/sphere/SPHERE_3538.DEC");
        Region reg = ImF.getRegion(0);
        ElementSurfSetHomogene mesh = (ElementSurfSetHomogene) reg.getElementSet();
         /*/
        ImportGmshMeshRegion ImF = new ImportGmshMeshRegion(0, "D:/Meshs/plaqueInclinee/plaqueIncline_2400e.msh");
//        ImportGmshMeshRegion ImF = new ImportGmshMeshRegion(0,"D:/Meshs/plaqueNonCoPlanaire/plaqueNCP_2800e.msh");
        ElementSurfSetHomogene mesh = (ElementSurfSetHomogene) ImF.createHomogeneESet();
        //*/
        System.out.println("Nombre elements " + mesh.getNbElement());
        System.out.println("Nombre de noeuds " + mesh.getNbNoeud());
        // Espace nodal gradHgrad
        GradNodalDeg1 gradalpha = new GradNodalDeg1(mesh);
        // Espace Fonctionnel pour T
        CurlSEdgeDeg1 rotSalpha = new CurlSEdgeDeg1(gradalpha);

        // Espace fonction fonction de forme
        NodalDeg1 alpha = new NodalDeg1(mesh);
        // Espace fcontionel alpha.N
        Hgrad alphaN = alpha.createProjOnNormal(1);

        System.out.println("Assemblage Hmatrix avec recomp.");
        double eps = 1e-5;
        int kmax = 50, nmin = 30, order = 6;
        GalerkinIntegralFormulationHCA IF = new GalerkinIntegralFormulationHCA(alphaN, rotSalpha, new CrossDGmulAlpha(1), new SelfElementFixedGauss(3, new Cancel()), 3, 3,
                eps, kmax, nmin, order, 2.0, true);
        long deb = System.nanoTime();
        IF.assembly();
        long fin = System.nanoTime();
        System.out.println("fin = " + (fin - deb) / 1e9 + "secs");
        StorageHmatrix H1 = (StorageHmatrix) IF.getStore();

        System.out.println("Assemblage Hmatrix sans recomp.");
        IF = new GalerkinIntegralFormulationHCA(alphaN, rotSalpha, new CrossDGmulAlpha(1), new SelfElementFixedGauss(3, new Cancel()), 3, 3,
                eps, kmax, nmin, order, 2.0, false);
        deb = System.nanoTime();
        IF.assembly();
        fin = System.nanoTime();
        System.out.println("fin = " + (fin - deb) / 1e9 + "secs");
        StorageHmatrix H2 = (StorageHmatrix) IF.getStore();

        System.out.println("Check Ranks: H1 vs H2");
        H1.checkRanks(H2);
    }
}
