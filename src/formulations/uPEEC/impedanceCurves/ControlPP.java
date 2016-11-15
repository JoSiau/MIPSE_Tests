/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package formulations.uPEEC.impedanceCurves;

import java.io.IOException;
import java.io.InputStream;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author jsiau
 */
public class ControlPP {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {
        ControlPP c = new ControlPP();

    }

    public ControlPP() throws Exception {
//        computeMain c[] = new computeMain[4];
//        c[0] = new computeMain(new String[]{"1e6", "4e7", "20"});
//        c[1] = new computeMain(new String[]{"4.001e7", "1e8", "30"});
//        c[2] = new computeMain(new String[]{"1.001e8", "5e8", "30"});
//        c[3] = new computeMain(new String[]{"5.001e8", "1e9", "30"});
//        for (int i = 0; i < 4; i++) {
//            ExecutorService s = Executors.newSingleThreadExecutor();
//            s.submit(c[i]);            
//        }

        String s[][] = new String[5][];
        s[0] = new String[]{"1e6", "4e7", "20"};
        s[1] = new String[]{"4.001e7", "1e8", "30"};
        s[2] = new String[]{"1.001e8", "2e8", "30"};
        s[3] = new String[]{"2.001e8", "5e8", "30"};
        s[4] = new String[]{"5.001e8", "1e9", "30"};

        Process process[] = new Process[s.length];
        final InputStream stdOut[] = new InputStream[s.length];
        for (int i = 0; i < s.length; i++) {

            String separator = System.getProperty("file.separator");
            String classpath = System.getProperty("java.class.path");
            String path = System.getProperty("java.home")
                    + separator + "bin" + separator + "java";
            ProcessBuilder processBuilder
                    = new ProcessBuilder(path, "-cp",
                            classpath,
                            Zf_PEEC_RLMPC_SURF_PP.class.getName(), s[i][0], s[i][1], s[i][2]);
            process[i] = processBuilder.start();
            stdOut[i] = process[i].getInputStream();
//            process.waitFor();
//            process.destroy();
        }

        for (int ii = 0; ii < s.length; ii++) {
            final int i = ii;
            new Thread(new Runnable() {
                public void run() {
                    byte[] buffer = new byte[8192];
                    int len = -1;
                    try {
                        while ((len = stdOut[i].read(buffer)) > 0) {
                            System.out.write(buffer, 0, len);
                        }
                    } catch (IOException ex) {
                        Logger.getLogger(ControlPP.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }).start();
        }
        
//        for (int i = 0; i < s.length; i++) {
//            process[i].waitFor();
//            System.out.println("stop waitFor: "+i);
//        }
//        for (int i = 0; i < s.length; i++) {
//            process[i].destroy();
//            System.out.println("Destroy: "+i);
//        }
    }

    
}
