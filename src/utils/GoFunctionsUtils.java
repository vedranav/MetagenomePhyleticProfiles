package utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import static utils.FileUtils.findReaderType;

/**
 *
 * @author Vedrana Vidulin [vedrana.vidulin@ijs.si]
 */
public class GoFunctionsUtils
{
    public static Set<Integer> extractOnlyProkaryoticGOs(Set<Integer> goFunctions, File geneOntologyFile) throws IOException
    {
        Set<Integer> prokaryoticGOs = new TreeSet<>();
        
        BufferedReader br = findReaderType(geneOntologyFile);
        
        Integer func = null;
        
        String line;
        
        while ((line = br.readLine()) != null)
        {
            line = line.trim();

            if (line.startsWith("<id>GO:"))
                func = Integer.valueOf(line.substring(line.indexOf(":") + 1, line.lastIndexOf("<")));
            else if (line.startsWith("<subset>") && func != null)
            {
                if (line.substring(line.indexOf(">") + 1, line.lastIndexOf("<")).trim().toLowerCase().equals("gosubset_prok") &&
                    goFunctions.contains(func))
                    prokaryoticGOs.add(func);
            }
        }
        
        return prokaryoticGOs;
    }
    
    public static Map<Integer, Double> loadGOFunctionFrequencies(File inFileWithUniprotFrequencies) throws IOException
    {
        Map<Integer, Double> func2Freq = new HashMap<>();
        
        BufferedReader br = findReaderType(inFileWithUniprotFrequencies);
        
        String line;

        while ((line=br.readLine()) != null)
            if (!line.startsWith("#"))
            {
                String[] parts = line.split("\t");

                int go = Integer.parseInt(parts[0]);
                double frequency = Double.parseDouble(parts[1]);

                func2Freq.put(go, frequency);                
            }
        
        return func2Freq;
    }
}
