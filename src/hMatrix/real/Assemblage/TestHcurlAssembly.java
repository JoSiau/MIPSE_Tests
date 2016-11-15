/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Assemblage;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.EdgeDeg1;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.ColumnVector;
import got.matrix.Matrix;

import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author jsiau
 */
public class TestHcurlAssembly {

    /**
     * @param args the command line arguments
     * @throws java.io.IOException
     */
    public static void main(String[] args) throws IOException {
        String meshDir = new java.io.File(".").getCanonicalPath();
        meshDir += "\\src\\g2elab\\mipse\\formulationInProgress\\magnetodynamic\\U_PEEC_DIELECTRIC\\";
        ImportFlux importM = new ImportFlux(meshDir + "COND_VOL16.DEC");
        //
        Region r = importM.getRegion(1);
        EdgeDeg1 fs = new EdgeDeg1((ElementSetHomogene) r.getElementSet());

        int nbGauss = 8, nbGaussCor = 27;

        GalerkinIntegralFormulationHCA as = new GalerkinIntegralFormulationHCA(fs, fs, new MultGvect(), new SelfElementFixedGauss(nbGauss, new InCreasedPGSourceNumber(nbGaussCor)), nbGauss, nbGaussCor,
                1e-4, 50, 50, 3);
        as.assembly();
        StorageHmatrix h = (StorageHmatrix) as.getStore();

        GalerkinIntegralFormulationFull fu = new GalerkinIntegralFormulationFull(fs, fs, new MultGvect(), new SelfElementFixedGauss(nbGauss, new InCreasedPGSourceNumber(nbGaussCor)), nbGauss);
        fu.assembly();
        Matrix m = ((StorageFull) fu.getStore()).getMatrix();

        double b[] = new double[fs.getActiveDofCount()];
        Arrays.fill(b, 1.0);

        ColumnVector xh = h.prod(new ColumnVector(b));

        ColumnVector xm = new ColumnVector(fs.getActiveDofCount());
        xm.mul(m, new ColumnVector(b));

        xh.sub(xm);
        System.out.println("Erreur = " + xh.norm() / xm.norm());

    }

}
