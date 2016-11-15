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
public class Control_BobZach {
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {
        Control_BobZach c = new Control_BobZach();
    }

    public Control_BobZach() throws Exception {
         String s[][] = new String[5][];
         
         s[0] = new String[]{"1e6", "5e6", "25"};
         s[1] = new String[]{"5.001e6", "1e7", "25"};
         s[2] = new String[]{"1.001e7", "5e7", "25"};
         s[3] = new String[]{"5.001e7", "1e8", "25"};
         
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
                            Zf_BobZacharie.class.getName(), s[i][0], s[i][1], s[i][2]);
            process[i] = processBuilder.start();
            stdOut[i] = process[i].getInputStream();
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
                        Logger.getLogger(Control_BobZach.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            }).start();
        }
    }

    
}
