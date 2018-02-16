package utils;

import com.google.common.collect.ArrayTable;
import com.google.common.collect.Table;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import static utils.FileUtils.findReaderType;

/**
 *
 * @author Vedrana Vidulin [vedrana.vidulin@ijs.si]
 */
public class ResourceLoaders
{
    public static Table<Integer, Integer, Double> loadPrScoresTable(File inFileWithPrScores, Set<Integer> loadOnlyFunctions,
                                                                    Set<Integer> loadOnlyOGs) throws IOException
    {
        List<Integer> ogs = new ArrayList<>();
        List<Integer> funcs = new ArrayList<>();        
        Table<Integer, Integer, Double> prScoresTable = null;
        
        //Collect gene families to initialize the table with Pr scores
        BufferedReader br = findReaderType(inFileWithPrScores);
        
        String line;
        
        while((line = br.readLine()) != null)
            if (!line.startsWith("OG"))
                ogs.add(ogStrToOgInt(line.substring(0, line.indexOf("\t"))));
        
        System.out.println("Gene families: " + ogs.size());
        
        //Collect functions for initializing the table, initialize and load Pr scores
        List<Integer> allFuncsWithPrList = new ArrayList<>();
        
        br = findReaderType(inFileWithPrScores);
        
        while((line = br.readLine()) != null)
        {
            String[] parts = line.split("\t");
            
            if (line.startsWith("OG"))
            {
                for (int i = 1; i < parts.length; i++)
                    allFuncsWithPrList.add(Integer.parseInt(parts[i]));
                
                Set<Integer> allFuncsWithPr = new TreeSet<>();
                allFuncsWithPr.addAll(allFuncsWithPrList);
                
                if (loadOnlyFunctions != null)
                    allFuncsWithPr.retainAll(loadOnlyFunctions);
                
                funcs.addAll(allFuncsWithPr);
                
                System.out.println("Functions: " + funcs.size());
                
                prScoresTable = ArrayTable.create(ogs, funcs);
                
                continue;
            }
            
            
            int og = ogStrToOgInt(line.substring(0, line.indexOf("\t")));
            
            if (loadOnlyOGs != null && !loadOnlyOGs.contains(og))
                continue;
            
            for (int i = 1; i < parts.length; i++)
            {
                int func = allFuncsWithPrList.get(i-1);
                
                if (funcs.contains(func))
                {
                    String val = parts[i].trim();
                    
                    if (!val.equals("NA"))
                        prScoresTable.put(og, func, Double.parseDouble(val));
                }
            }
        }
        
        return prScoresTable;
    }
    
    public static Integer ogStrToOgInt(String ogStr)
    {
        Pattern p = Pattern.compile("([a-z]+)([0-9]+)");
        Matcher m = p.matcher(ogStr.toLowerCase());
        
        if (m.find())
        {
            String ogType = m.group(1);
            int ogId = Integer.valueOf(m.group(2));
            return (ogType.equals("nog") ? -ogId : ogId);
        }
        
        return null;
    }
    
    public static Map<Integer, Set<Integer>> extractOg2FuncsFromPrTable(Table<Integer, Integer, Double> prScoresTable, double prThreshold)
    {
        Map<Integer, Set<Integer>> og2funcs = new TreeMap<>();
        
        for (int og : prScoresTable.rowKeySet())
        {
            Set<Integer> funcs = new TreeSet<>();
            for (int func : prScoresTable.columnKeySet())
            {
                Double val = prScoresTable.get(og, func);
                
                if (val != null && val >= prThreshold)
                    funcs.add(func);
            }
            
            og2funcs.put(og, funcs);
        }
        
        return og2funcs;
    }
    
    public static Map<Integer, Set<Integer>> loadOg2FunctionsFromFile(File inOg2FunctionsFile) throws IOException
    {
        Map<Integer, Set<Integer>> og2functions = new TreeMap<>();
        
        Set<Integer> allFunctions = new HashSet<>(); 
        
        BufferedReader br = findReaderType(inOg2FunctionsFile);
        
        String line;
        
        while ((line = br.readLine()) != null)
        {
            if (line.startsWith("@") || line.isEmpty() || line.startsWith("#"))
                continue;
            
            int og;
            String functionsStr;
            
            if (inOg2FunctionsFile.getName().endsWith(".arff.zip") || inOg2FunctionsFile.getName().endsWith(".arff"))
            {
                og = ogStrToOgInt(line.substring(0, line.indexOf(",")));
                functionsStr = line.substring(line.lastIndexOf(",") + 1).trim();
            }
            else
            {
                og = Integer.parseInt(line.substring(0, line.indexOf("\t")));
                functionsStr = line.substring(line.indexOf("\t") + 1);
            }
            
            String[] parts = functionsStr.split("@");
                
            Set<Integer> functions = new TreeSet<>();
            for (String part : parts)
                if (!part.isEmpty())
                    functions.add(Integer.parseInt(part));

            if (og != 0 && !functions.isEmpty())
            {
                og2functions.put(og, functions);
                allFunctions.addAll(functions);
            }
        }
        
        System.out.println("Gene families with known functions: " + og2functions.size());
        System.out.println("Known functions: " + allFunctions.size());
        
        return og2functions;
    }
}
