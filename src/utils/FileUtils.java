package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Enumeration;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

/**
 *
 * @author Vedrana Vidulin [vedrana.vidulin@ijs.si]
 */
public class FileUtils
{
    public static BufferedReader findReaderType(File inFile) throws IOException
    {
        BufferedReader br = null;
        
        if (inFile.getAbsolutePath().endsWith("gz"))
        {
            GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(inFile));       
            br = new BufferedReader(new InputStreamReader(gzip));
        }
        else if (inFile.getAbsolutePath().endsWith(".zip"))
        {
            ZipFile zipFile = new ZipFile(inFile);

            Enumeration<? extends ZipEntry> entries = zipFile.entries();

            while(entries.hasMoreElements()) //considers only the first file in the zip archive
            {
                ZipEntry entry = entries.nextElement();
                InputStream stream = zipFile.getInputStream(entry);
            
                br = new BufferedReader(new InputStreamReader(stream));
                
                break;
            }
        }
        else
            br = new BufferedReader(new FileReader(inFile));
        
        return br;
    }
}
