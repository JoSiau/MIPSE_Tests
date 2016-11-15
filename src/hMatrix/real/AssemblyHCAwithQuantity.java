/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.elements.ElementSetHomogene;
import g2elab.mipse.meshCore.functionSpace.Cell;
import g2elab.mipse.meshCore.functionSpace.FaceDeg1;
import g2elab.mipse.meshCore.functionSpace.NodalDeg1;
import g2elab.mipse.meshCore.quantity.RealScalarCellQuantity;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.AnalyticalCorrection;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.InCreasedPGSourceNumber;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultG;
import g2elab.mipse.mipseCore.integralIntegration.kernel.MultGvect;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationHCA;
import g2elab.mipse.mipseCore.storage.StorageFull;
import g2elab.mipse.mipseCore.storage.StorageHmatrix;
import got.matrix.Matrix;
import got.matrix.RowVector;

import java.io.IOException;

/**
 * @author jsiau
 */
public class AssemblyHCAwithQuantity {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        String meshDir = new java.io.File(".").getCanonicalPath();
        //*
        meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
        //        ImportFlux mesh = new ImportFlux(meshDir + "/sphere/SPHERE_1884.DEC");
        ImportFlux mesh = new ImportFlux(meshDir + "/PLAQUE2000.DEC");
        //        ImportFlux mesh = new ImportFlux("VariateurDEC/ATV71_105963.DEC");
        ElementSetHomogene ES = (ElementSetHomogene) mesh.getRegion(0).getElementSet();
        /*/
         meshDir += "/src/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh";
         //         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/plaque/plaque_4548.msh");   
         //         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0,meshDir+"/cube.msh");   
         ImportGmshMeshRegion mesh1 = new ImportGmshMeshRegion(0, meshDir + "/sphere/SPHERE_4893.msh");
         ElementSetHomogene ES = mesh1.createHomogeneESet();
         //*/
        int d = 3;
        int n = ES.getNbNoeud();

        System.out.println("nbNoeuds= " + ES.getNbNoeud());
        System.out.println("nbElmts= " + ES.getNbElement());

        RowVector q = new RowVector(ES.getNbElement());
        q.setAllElements(100);
        for (int i = 0; i < 10; i++) {
            q.setElement(i, 27);
        }

        System.out.println("////////////////////////////////////////////////////");
        System.out.println("////////////////        CELL        ////////////////");
        System.out.println("////////////////////////////////////////////////////");
        Cell C = new Cell(ES);

        RealScalarCellQuantity quant = new RealScalarCellQuantity(new Cell(ES), q);

        int epsPow = 5;
        GalerkinIntegralFormulationHCA f = new GalerkinIntegralFormulationHCA(C, C, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7, 7,
                Math.pow(10, -epsPow), 50, 30, epsPow - 1);
        f.assembly(quant);
        StorageHmatrix H = (StorageHmatrix) f.getStore();

        GalerkinIntegralFormulationFull IF = new GalerkinIntegralFormulationFull(C, C, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7);
        IF.assembly(quant);
        Matrix M = ((StorageFull) IF.getStore()).getMatrix();

        M = H.reOrder_Full2HM(M);
        H.CheckError(M, Math.pow(10, -epsPow));


        System.out.println("////////////////////////////////////////////////////");
        System.out.println("////////////////        HGRAD        ///////////////");
        System.out.println("////////////////////////////////////////////////////");
        NodalDeg1 ND1 = new NodalDeg1(ES);
        f = new GalerkinIntegralFormulationHCA(ND1, ND1, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7, 7,
                Math.pow(10, -epsPow), 50, 30, epsPow - 1);
        f.assembly(quant);
        H = (StorageHmatrix) f.getStore();

        IF = new GalerkinIntegralFormulationFull(ND1, ND1, new MultG(), new SelfElementFixedGauss(7, new AnalyticalCorrection()), 7);
        IF.assembly(quant);
        M = ((StorageFull) IF.getStore()).getMatrix();

        M = H.reOrder_Full2HM(M);
        H.CheckError(M, Math.pow(10, -epsPow));

        System.out.println("////////////////////////////////////////////////////");
        System.out.println("////////////////        HDIV        ////////////////");
        System.out.println("////////////////////////////////////////////////////");
        FaceDeg1 FD1 = new FaceDeg1(ES);

        f = new GalerkinIntegralFormulationHCA(FD1, FD1, new MultGvect(), new SelfElementFixedGauss(7, new InCreasedPGSourceNumber(3)), 7, 7,
                Math.pow(10, -epsPow), 50, 30, epsPow - 1);
        f.assembly(quant);
        H = (StorageHmatrix) f.getStore();

        IF = new GalerkinIntegralFormulationFull(FD1, FD1, new MultGvect(), new SelfElementFixedGauss(7, new InCreasedPGSourceNumber(3)), 7);
        IF.assembly(quant);
        M = ((StorageFull) IF.getStore()).getMatrix();

        M = H.reOrder_Full2HM(M);
        H.CheckError(M, Math.pow(10, -epsPow));
    }

}
