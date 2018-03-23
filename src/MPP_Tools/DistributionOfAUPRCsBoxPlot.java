package MPP_Tools;

import com.google.common.collect.ArrayTable;
import com.google.common.collect.Table;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import static org.apache.commons.math3.util.Precision.round;
import static utils.FileUtils.findReaderType;
import static utils.GoFunctionsUtils.extractOnlyProkaryoticGOs;
import static utils.GoFunctionsUtils.loadGOFunctionFrequencies;
import static utils.RUtils.executeRScript;


/**
 *
 * @author Vedrana Vidulin [vedrana.vidulin@ijs.si]
 */
public class DistributionOfAUPRCsBoxPlot
{
    /**
     * Draws box-plots, each showing a distribution of GO function-level AUPRCs for a specific classification model and
     * a specific level of GO functions generality. Generality is defined as information content (IC), where high numbers
     * indicate that a GO function is specific and, therefore, contributes more information about gene's function.
     * This code assumes that you have R on your computer.
     * 
     * @param inFilesWithAUPRCsAndPredictionsStats      An array of tab delimited files, one for each classifier. 
     *                                                  Tab delimited file lists for each GO function an AUPRC for a specific 
     *                                                  classifier and an information about the number of gene families for 
     *                                                  which the function was predicted at Pr>=50% by the classifier.
     * @param classifierNames                           An array with classifier names.
     * @param colors                                    An array of colors that will be used to color box plots of different
     *                                                  classifiers. Use R specific color definitions.
     * @param inGeneOntologyFile                        Gene ontology in obo-xml.gz format.
     * @param inFileWithFunctionFrequencies             File with frequencies of GO functions appearance in Uniprot-GOA.
     * @param outFolder                                 Folder that will contain graph and statistics files.
     * 
     * @throws IOException
     * @throws InterruptedException 
     */
    public static void drawDistributionOfAUPRCsBoxPlot(File[] inFilesWithAUPRCsAndPredictionsStats, String[] classifierNames,
                                                       String[] colors, File inGeneOntologyFile, File inFileWithFunctionFrequencies,
                                                       File outFolder) throws IOException, InterruptedException
    {
        if (!outFolder.exists())
            outFolder.mkdirs();
        
        //Extract only learnable prokaryotic GO functions from the GO subset gosubset_prok
        Map<Integer, List<Double>> funcAUPRCs = new TreeMap<>();
        Map<Integer, List<Integer>> funcStats = new TreeMap<>();
        
        for (File f : inFilesWithAUPRCsAndPredictionsStats)
        {
            BufferedReader br = findReaderType(f);
            
            String line;
            
            while((line=br.readLine()) != null)
                if (!line.startsWith("#"))
                {
                    String[] parts = line.split("\t");
                    
                    int func = Integer.parseInt(parts[0]);
                    double auprc = Double.parseDouble(parts[1]);
                    int numPredictions = Integer.parseInt(parts[2]);
                    
                    List<Double> auprcs = new ArrayList<>();
                    if (funcAUPRCs.containsKey(func))
                        auprcs = funcAUPRCs.get(func);
                    
                    auprcs.add(auprc);
                    
                    funcAUPRCs.put(func, auprcs);
                    
                    
                    List<Integer> stats = new ArrayList<>();
                    if (funcStats.containsKey(func))
                        stats = funcStats.get(func);
                    
                    stats.add(numPredictions);
                    
                    funcStats.put(func, stats);
                }
        }
        
        Set<Integer> prokFuncs = extractOnlyProkaryoticGOs(funcAUPRCs.keySet(), inGeneOntologyFile);
        funcAUPRCs.keySet().retainAll(prokFuncs); //Remove non-prokaryotic GO functions
        
        Set<Integer> nonLearnableFuncs = new HashSet<>();
        for (int func : funcAUPRCs.keySet())
        {
            List<Integer> stats = funcStats.get(func);
            
            int sum = 0;
            for (int stat : stats)
                sum += stat;
            
            if (sum == 0)
                nonLearnableFuncs.add(func);
        }
        
        funcAUPRCs.keySet().removeAll(nonLearnableFuncs); //Remove non-learnable GO functions
        
        
        //Divide functions in categories according to their level of generality
        Map<Integer, Double> funcFreq = loadGOFunctionFrequencies(inFileWithFunctionFrequencies);
        
        String[] generalityLevels = {"S", "M", "G"};
        List<String> methodGenerality = new ArrayList<>();
        for (String generalityLevel : generalityLevels)
            for (String method : classifierNames)            
                methodGenerality.add(method + "-" + generalityLevel);
        
        Table<Integer, String, Double> function_MethodGenerality_AUPRCs = ArrayTable.create(funcAUPRCs.keySet(), methodGenerality);
        
        for (int func : funcAUPRCs.keySet())
        {
            List<Double> auprcs = funcAUPRCs.get(func);
            
            for (int i = 0; i < auprcs.size(); i++)
            {
                double IC = -(Math.log(funcFreq.get(func)) / Math.log(2));
                
                String generalityCategory;
                if (IC > 8)
                    generalityCategory = "S";
                else if (IC >= 4)
                    generalityCategory = "M";
                else
                    generalityCategory = "G";
                
                function_MethodGenerality_AUPRCs.put(func, classifierNames[i] + "-" + generalityCategory, auprcs.get(i));
            }
        }
        
        
        //Write the table to file
        String dataFile = outFolder.getAbsolutePath() + "/AUPRCs_data.txt";
        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(dataFile)));
        
        bw.write("Function");
        for (String mg : function_MethodGenerality_AUPRCs.columnKeySet())
            bw.write("\t" + mg);
        bw.write("\n");
        
        for (int func : function_MethodGenerality_AUPRCs.rowKeySet())
        {
            bw.write(String.valueOf(func));
            
            for (String mg : function_MethodGenerality_AUPRCs.columnKeySet())
            {
                Double auprc = function_MethodGenerality_AUPRCs.get(func, mg);
                bw.write("\t" + (auprc == null ? "NA" : String.valueOf(round(auprc, 4))));
            }
            
            bw.write("\n");
        }
        
        bw.close();
        
        
        //Write file with statistics and legend
        bw = new BufferedWriter(new FileWriter(new File(outFolder + "/Graph_statistics+Legend.txt")));
        
        bw.write(function_MethodGenerality_AUPRCs.rowKeySet().size() + " learnable GO functions with at least one prediction at Pr>=50%\n");
                
        bw.write("\nLEGEND:\n");
        bw.write("S - specific GO functions with IC>8\n");
        bw.write("M - medium specific GO functions with 4<=IC<=8\n");
        bw.write("G - general GO functions with IC<4\n");
        
        bw.close();
        
        
        //Generate R script
        String generalityRScript = outFolder.getAbsolutePath() + "/Rscript.r";
        bw = new BufferedWriter(new FileWriter(new File(generalityRScript)));
        
        bw.write("pdf('" + outFolder.getAbsolutePath().replace('\\', '/') + "/BoxPlot.pdf')\n");
        bw.write("data = read.table('" + dataFile.replace('\\', '/') + "', header = T, sep = '\t', na.strings = 'NA', check.names = F)\n");
        bw.write("par(cex.axis=0.9)\n");
        bw.write("par(cex.lab=1.2)\n");
        bw.write("boxplot(data[c(2:ncol(data))], las = 2, at = c(");
        String delimiter = "";
        for (int i = 1; i <= (classifierNames.length * 3) + 2; i++)
            if (i % (classifierNames.length + 1) != 0)
            {
                bw.write(delimiter + i);
                delimiter = ",";
            }
        
        bw.write("), col = c(");
        delimiter = "";
        for (int i = 0; i < 3; i++)
        {
            bw.write(delimiter);
            
            String del = "";
            for (int j = 0; j < classifierNames.length; j++)
            {
                bw.write(del + "'" + colors[j] + "'");
                del = ",";
            }
            
            delimiter = ",";
        }
        
        bw.write("), main = '', ylab = 'AUPRC', notch = T, varwidth = T, outline = T)\n");
        
        bw.write("dev.off()\n");
        
        bw.close();
        
        executeRScript(generalityRScript);
        
        new File(generalityRScript).delete();
    }
}
