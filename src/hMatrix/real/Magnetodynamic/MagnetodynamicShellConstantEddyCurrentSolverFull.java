/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package hMatrix.real.Magnetodynamic;

import g2elab.mipse.meshCore.IO.flux.ImportFlux;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshCell;
import g2elab.mipse.meshCore.IO.gmsh.ExportGmshHgrad;
import g2elab.mipse.meshCore.IO.paraview.ExportVtkCell;
import g2elab.mipse.meshCore.IO.paraview.ExportVtkHgrad;
import g2elab.mipse.meshCore.contraints.InactiveBorderNode;
import g2elab.mipse.meshCore.elements.ElementSurfSetHomogene;
import g2elab.mipse.meshCore.functionSpace.*;
import g2elab.mipse.meshCore.quantity.*;
import g2elab.mipse.meshCore.region.Region;
import g2elab.mipse.mipseCore.integralIntegration.correctionstrategies.typecorrection.Cancel;
import g2elab.mipse.mipseCore.integralIntegration.integrationcorrectionstrategies.SelfElementFixedGauss;
import g2elab.mipse.mipseCore.integralIntegration.kernel.CrossDGmulAlpha;
import g2elab.mipse.mipseCore.numericalMethods.FiniteElementFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulation;
import g2elab.mipse.mipseCore.numericalMethods.GalerkinIntegralFormulationFull;
import g2elab.mipse.mipseCore.storage.Storage;
import g2elab.mipse.mipseCore.storage.StorageSparse;
import g2elab.mipse.numericalTools.iterativeSolver.ProductComplex;
import g2elab.mipse.numericalTools.iterativeSolver.complex.FGMResComplex;
import g2elab.mipse.numericalTools.matrix.complex.MatrixComplexPartReIm;
import g2elab.mipse.numericalTools.matrix.real.sparse.rcReal.SparseMatrixRCReal;
import g2elab.mipse.numericalTools.preconditioner.ComplexPrecond;
import got.matrix.Matrix;
import got.matrix.RowVector;

import java.io.File;

/**
 * @author chadebec
 */
public class MagnetodynamicShellConstantEddyCurrentSolverFull implements ProductComplex, ComplexPrecond {

    private SparseMatrixRCReal zbRe;
    private Storage zbIm;
    private MatrixComplexPartReIm zb;
    private ComplexNodalQuantity Res;
    /*
     * Second membre 
     */
    private double[] secondMembre;
    /*
     * nombre de degres de liberte
     */
    /*
     * espace fonctionnel des fonctions de forme
     */
    private double frequence;
    private double sigma;
    private NodalDeg1 alpha;
    private int nbdof;
    private double epaisseur;
    static double mu0 = 4 * Math.PI * 1e-7;

    public MagnetodynamicShellConstantEddyCurrentSolverFull(ElementSurfSetHomogene mesh, RealVectorCellQuantity H0cell, double frequence, double sigma, double epaisseur) {

        this.frequence = frequence;
        this.sigma = sigma;
        this.epaisseur = epaisseur;

        // Constrainte de dirichlet null sur la frontiere
        InactiveBorderNode Dir0border = new InactiveBorderNode();

        // Espace nodal gradHgrad
        GradNodalDeg1 gradalpha = new GradNodalDeg1(mesh, Dir0border);//,new InactiveLastNode());//, Dir0border);
        // Espace Fonctionnel pour T
        CurlSEdgeDeg1 rotSalpha = new CurlSEdgeDeg1(gradalpha);
        //new CurlSHcurl(mesh, 1, Dir0border);

        // Espace fonction fonction de forme
        alpha = new NodalDeg1(mesh, Dir0border);//,new InactiveLastNode());//, Dir0border);
        nbdof = alpha.getActiveDofCount() * 2;

        // Espace fcontionel alpha.N
        Hgrad alphaN = alpha.createProjOnNormal(1);

        System.out.println("debut integration FEM");
        // support pour les conductivites
        Matrix resVal = new Matrix(1, this.alpha.getElementSet().getNbElement());
        resVal.setAllElements(1 / this.sigma);
        Cell supportConductivite = new Cell(alpha.getElementSet());
        RealScalarCellQuantity invSigma = new RealScalarCellQuantity(supportConductivite, resVal);
        FiniteElementFormulation FEM = new FiniteElementFormulation(gradalpha);
        FEM.assembly(1, invSigma);
        StorageSparse SFE = (StorageSparse) FEM.getStore();
        this.zbRe = SFE.getMatrixPrecond(null).getSparseMatrixRCReal();
        System.out.println("FEM integres");

        System.out.println("debut integral IV");
        GalerkinIntegralFormulation IF = new GalerkinIntegralFormulationFull(alphaN, rotSalpha, new CrossDGmulAlpha(2 * Math.PI * this.frequence * mu0 * this.epaisseur), new SelfElementFixedGauss(3, new Cancel()), 3);
        long deb = System.nanoTime();
        IF.assembly();
        long fin = System.nanoTime();

        System.out.println("fin = " + (fin - deb) / 1e9 + "secs");
        //IF.FMMAssembly(alphaN, 3,3, new RepartitionElemNiveauConstant(0), 1,4,2,1);
        this.zbIm = IF.getStore();
        System.out.println("IV integres");

//        System.out.println(((StorageFull)IF.getStore()).getMatrix());

        // Second membre
        RealScalarCellQuantity H0Ncell = H0cell.normalProjection();
        RowVector h0n = H0Ncell.projWithExplicitDof(alpha, 3);
        h0n.scale(-1 * 2 * Math.PI * frequence * mu0);
        double[] bim = new double[alpha.getActiveDofCount()];
        h0n.get(bim);
        secondMembre = new double[2 * alpha.getActiveDofCount()];
        for (int i = 0; i < alpha.getActiveDofCount(); i++) {
            secondMembre[2 * i + 1] = bim[i];
        }

    }

    public ComplexNodalQuantity solve() {
        FGMResComplex gmres = new FGMResComplex(this, this);
        gmres.setInfoResolution(new double[]{100, 1e-12, 1, -50});
        double[] x = new double[nbdof];
        gmres.solve(x, secondMembre);
        Matrix res = new Matrix(alpha.getActiveDofCount(), 2, x);
        RealNodalQuantity Treal = new RealNodalQuantity(alpha, res.transpose().row(0));
        RealNodalQuantity Timag = new RealNodalQuantity(alpha, res.transpose().row(1));
        ComplexNodalQuantity T = new ComplexNodalQuantity(alpha);
        T.setReal(Treal);
        T.setImag(Timag);
        Res = T;
        return T;
    }

    @Override
    public double[] product(double[] x, double[] mult) {
        zb = new MatrixComplexPartReIm();
        zb.setPartReIm(zbRe, this.zbIm);
        mult = zb.product(x, new double[x.length]);
        return mult;
    }

    @Override
    public double[] precond(double[] x, double[] mult) {
        if (mult == null) {
            mult = new double[x.length];
        }
        System.arraycopy(x, 0, mult, 0, x.length);
        return mult;
    }

    @Override
    public long getMemoryUsed() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public static void main(String[] args) throws Exception {
        String className = "";
        if (args.length != 0) {
            className = args[0];
        }

        File f = new File("");
        String path = f.getAbsolutePath();
        String file = path + "/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/PLAQUE2000.DEC";
//        String file = path + "./../mipsecore/src/trunk/g2elab/mipse/mipseCore/matrixCompression/hMatrix/mesh/sphere/SPHERE_3538.DEC";
        //String file = path + "/src/g2elab/mipse/formulationInProgress/magnetodynamic/PLAQUESANSTROU_68.DEC";
        ImportFlux IF = new ImportFlux(file);
        Region reg = IF.getRegion(0);
        ElementSurfSetHomogene mesh = (ElementSurfSetHomogene) reg.getElementSet();

        System.out.println("Nombre elements " + mesh.getNbElement());
        System.out.println("Nombre de noeuds " + mesh.getNbNoeud());

        double mu = 4 * Math.PI * Math.pow(10, -7);
        double sigma = 6 * Math.pow(10, 7);
        double freq = 0.01;
        double epaisseur = 50e-3;

        System.out.println(" Conductivité électrique = " + sigma + " S/m");
        System.out.println(" Perméabilité = " + mu + " H/m");
        System.out.println(" Fréquence = " + freq + " Hz");
        System.out.println(" ");


        Matrix Hext = new Matrix(3, mesh.getNbElement());
        for (int i = 0; i < mesh.getNbElement(); i++) {
            Hext.setElement(2, i, Math.sqrt(2));
        }
        RealVectorCellQuantity H0cell = new RealVectorCellQuantity(new Cell(mesh));
        H0cell.setFullMatrixValues(Hext);
        // Projection sur les normales du maillage

        MagnetodynamicShellConstantEddyCurrentSolverFull formulation = new MagnetodynamicShellConstantEddyCurrentSolverFull(mesh, H0cell, freq, sigma, epaisseur);
        ComplexNodalQuantity T = formulation.solve();

        System.out.println(T.getMatrixValues());

        NodalDeg1 FS = (NodalDeg1) T.getFS();

        RealNodalQuantity Treal = T.getReal();
        RealNodalQuantity Timag = T.getImag();

        RealVectorCellQuantity Jreal = Treal.curlS();
        RealVectorCellQuantity Jimag = Timag.curlS();

        ComplexVectorCellQuantity J = new ComplexVectorCellQuantity((Cell) Jreal.getFS());
        J.setReal(Jreal);
        J.setImag(Jimag);

        ExportVtkHgrad exportT = new ExportVtkHgrad(FS, "Y:/problemT.vtk");
        exportT.addQuantity(T, "T");

        ExportVtkCell exportJ = new ExportVtkCell((Cell) J.getFS(), "Y:/problemJ.vtk");
        exportJ.addQuantity(J, "J");

        ExportGmshHgrad exportTm = new ExportGmshHgrad(FS, "Y:/problemT.msh");
        exportTm.addQuantity(T, "T");

        ExportGmshCell exportJm = new ExportGmshCell((Cell) J.getFS(), "Y:/problemJ.msh");
        exportJm.addQuantity(J, "J");

        ExportGmshHgrad test = new ExportGmshHgrad(FS, "Y:/problemTcell.msh");
        test.vizualiseCellValue(Timag);


        g2elab.mipse.tools.multiThreads.GestionnaireTaches.getGestionnaireTaches().stop();
    }

    public void printRes() {
        NodalDeg1 FS = (NodalDeg1) Res.getFS();

        RealNodalQuantity Treal = Res.getReal();
        RealNodalQuantity Timag = Res.getImag();

        RealVectorCellQuantity Jreal = Treal.curlS();
        RealVectorCellQuantity Jimag = Timag.curlS();

        ComplexVectorCellQuantity J = new ComplexVectorCellQuantity((Cell) Jreal.getFS());
        J.setReal(Jreal);
        J.setImag(Jimag);

        ExportVtkHgrad exportT = new ExportVtkHgrad(FS, "Y:/problemT.vtk");
        exportT.addQuantity(Res, "T");

        ExportVtkCell exportJ = new ExportVtkCell((Cell) J.getFS(), "Y:/problemJ.vtk");
        exportJ.addQuantity(J, "J");

        ExportGmshHgrad exportTm = new ExportGmshHgrad(FS, "Y:/problemT.msh");
        exportTm.addQuantity(Res, "T");

        ExportGmshCell exportJm = new ExportGmshCell((Cell) J.getFS(), "Y:/problemJ.msh");
        exportJm.addQuantity(J, "J");

        ExportGmshHgrad test = new ExportGmshHgrad(FS, "Y:/problemTcell.msh");
        test.vizualiseCellValue(Timag);
    }

    @Override
    public void configuration(double[] config) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setMatrix(Object mat, double[] configuration) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
