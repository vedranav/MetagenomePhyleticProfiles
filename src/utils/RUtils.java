package utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 *
 * @author Vedrana Vidulin [vedrana.vidulin@ijs.si]
 */
public class RUtils
{
    public static void executeRScript(String rscriptFileName) throws IOException, InterruptedException
    {
        Runtime rt = Runtime.getRuntime();
        String[] rscript = {"Rscript", rscriptFileName};

        Process p = rt.exec(rscript);
        p.waitFor();

        InputStream stdout = p.getInputStream();
        InputStream stderr = p.getErrorStream();   			    

        //Clean up if any output in stdout
        BufferedReader brCleanUp = new BufferedReader(new InputStreamReader(stdout));

        String line;
        
        while ((line = brCleanUp.readLine ()) != null)
            System.out.println ("[Stdout] " + line);

        brCleanUp.close();

        //Clean up if any output in stderr
        brCleanUp = new BufferedReader(new InputStreamReader(stderr));
        while ((line = brCleanUp.readLine ()) != null)
            System.out.println ("[Stderr] " + line);

        brCleanUp.close();

        p.destroy();
    }
}
