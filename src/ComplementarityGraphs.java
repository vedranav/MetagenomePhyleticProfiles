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
import java.util.TreeSet;
import static org.apache.commons.math3.util.Precision.round;
import static utils.FileUtils.findReaderType;
import static utils.GoFunctionsUtils.extractOnlyProkaryoticGOs;
import static utils.RUtils.executeRScript;
import static utils.ResourceLoaders.extractOg2FuncsFromPrTable;
import static utils.ResourceLoaders.loadOg2FunctionsFromFile;
import static utils.ResourceLoaders.loadPrScoresTable;
import static utils.ResourceLoaders.ogToStr;

/**
 *
 * @author Vedrana Vidulin [vedrana.vidulin@ijs.si]
 */
public class ComplementarityGraphs
{
    public static void main(String[] args) throws Exception
    {
        //IMPORTANT!!!!
        //This code assumes that you have R already installed on your computer.
        //If you get an error "Cannot run program "Rscript": error=2, No such file or directory", you should
        //set the path to Rscript in "utils.RUtils" by adding the path to the line
        //"String[] rscript = {"/ADD PATH HERE/Rscript", rscriptFileName};".
        //You can find out the path by running "type -a Rscript" in a terminal (tested on Ubuntu and MacOS).
        
        //-----------------------------------
        //SELECT THE EXAMPLE FROM THE PAPER
        boolean Fig1_d_e = true;
        boolean Fig1_fgh_A = true;
        boolean Fig1_fgh_B = true;
        boolean Fig1_fgh_C = true;
        //-----------------------------------
        
        //-----------------------------------
        //SET PATHS
        String dataDir = ".../MetagenomePhyleticProfiles/src/data/";
        String outDir = "...";
        //-----------------------------------
        

        double[] prThresholds = {0.5, 0.7, 0.9};
        
        if (Fig1_d_e)
            for (double prThreshold : prThresholds)
            {
                File inFirstMethodPrecisionFile = new File(dataDir + "MPP-H_predictions.tsv.gz");
                File inSecondMethodPrecisionFile = new File(dataDir + "MPP-O_predictions.tsv.gz");
                String firstMethodName = "MPP-H";
                String secondMethodName = "MPP-O";
                File inOg2FunctionsFile = new File(dataDir + "og2funcs-eggNOG_3.txt.gz");
                File inGeneOntologyFile = new File(dataDir + "go_201401-termdb.obo-xml.gz");
                File outFolderForFig1d = new File(outDir + "/Fig1d/PR-" + prThreshold);
                                                                                    
                drawFunctionBasedComplementarityGraphForTwoMethods(inFirstMethodPrecisionFile, inSecondMethodPrecisionFile,
                    firstMethodName, secondMethodName, "'yellow','green','blue'", prThreshold, inOg2FunctionsFile,
                    inGeneOntologyFile, outFolderForFig1d);
                
                drawGeneFamilyBasedComplementarityGraphForTwoMethods(inFirstMethodPrecisionFile, inSecondMethodPrecisionFile,
                    new File(outFolderForFig1d + "/Functions_predicted_by_" + firstMethodName + "_and_" + secondMethodName + ".txt"),
                    firstMethodName, secondMethodName, ".7, .17", "yellow, blue", prThreshold, inOg2FunctionsFile,
                    inGeneOntologyFile, new File(outDir + "/Fig1e/PR-" + prThreshold));
            }
        
        if (Fig1_fgh_A)
            for (double prThreshold : prThresholds)
                drawFunctionBasedComplementarityGraphForTwoMethods(new File(dataDir + "MPP-H_predictions.tsv.gz"),
                    new File(dataDir + "PP-H_predictions.tsv.gz"), "MPP-H", "PP", "'red','green','blue'", prThreshold,
                    new File(dataDir + "og2funcs-eggNOG_3.txt.gz"), new File(dataDir + "go_201401-termdb.obo-xml.gz"),
                    new File(outDir + "/Fig1_fgh_A/PR-" + prThreshold));
        
        if (Fig1_fgh_B)
            for (double prThreshold : prThresholds)
                drawFunctionBasedComplementarityGraphForTwoMethods(new File(dataDir + "MPP-O_predictions.tsv.gz"),
                    new File(dataDir + "PP-O_predictions.tsv.gz"), "MPP-O", "PP", "'red','green','blue'", prThreshold,
                    new File(dataDir + "og2funcs-eggNOG_3.txt.gz"), new File(dataDir + "go_201401-termdb.obo-xml.gz"),
                    new File(outDir + "/Fig1_fgh_B/PR-" + prThreshold));
        
        if (Fig1_fgh_C)
            for (double prThreshold : prThresholds)
                drawFunctionBasedComplementarityGraphForTwoMethods(new File(dataDir + "MPP-I_predictions.tsv.gz"),
                    new File(dataDir + "PP-I_predictions.tsv.gz"), "MPP-I", "PP", "'red','green','blue'", prThreshold,
                    new File(dataDir + "og2funcs-eggNOG_4.txt.gz"), new File(dataDir + "go_201401-termdb.obo-xml.gz"),
                    new File(outDir + "/Fig1_fgh_C/PR-" + prThreshold));
    }
    
    /**
     * Draws graph that shows the level of complementarity between two classification models in terms of functions
     * they are able to predict at a specific level of precision (Pr).
     * This code assumes that you have R on your computer.
     * 
     * @param inFirstMethodPrecisionFile    The table with Pr scores outputted by the first classification model.
     * @param inSecondMethodPrecisionFile   The table with Pr scored outputted by the second classification model.
     * @param firstMethodName               The name of the first method, i.e., classification model.
     * @param secondMethodName              The name of the second method.
     * @param colors                        Color scheme, composed of three colors. Accepts color definitions from the R color pallete.
     * @param prThreshold                   Predictions with Pr >= prThreshold will be considered as positive.
     * @param inOg2FunctionsFile            File with known functions from Uniprot-GOA that are assigned to gene families.
     * @param inGeneOntologyFile            Gene ontology in obo-xml.gz format.
     * @param outFolder                     Folder that will contain graph and statistics files.
     * 
     * @throws IOException
     * @throws InterruptedException 
     */
    public static void drawFunctionBasedComplementarityGraphForTwoMethods(File inFirstMethodPrecisionFile,
                                                                          File inSecondMethodPrecisionFile,
                                                                          String firstMethodName, String secondMethodName,
                                                                          String colors, double prThreshold,
                                                                          File inOg2FunctionsFile, File inGeneOntologyFile,
                                                                          File outFolder) throws IOException, InterruptedException
    {
        if (!outFolder.exists())
            outFolder.mkdirs();
        
        //Load known functions
        Map<Integer, Set<Integer>> og2known_funcs = loadOg2FunctionsFromFile(inOg2FunctionsFile);
        
        //Consider only GO functions from prokaryotic GO subset gosubset_prok
        Set<Integer> knownProkFuncs = new TreeSet<>();
        for (Set<Integer> known_funcs : og2known_funcs.values())
            knownProkFuncs.addAll(known_funcs);
        
        knownProkFuncs = extractOnlyProkaryoticGOs(knownProkFuncs, inGeneOntologyFile);
        
        //Load predictions
        Map<Integer, Set<Integer>> og2predicted_funcs_First = extractOg2FuncsFromPrTable(
                loadPrScoresTable(inFirstMethodPrecisionFile, knownProkFuncs, og2known_funcs.keySet()), prThreshold);
        
        Map<Integer, Set<Integer>> og2predicted_funcs_Second = extractOg2FuncsFromPrTable(
                loadPrScoresTable(inSecondMethodPrecisionFile, knownProkFuncs, og2known_funcs.keySet()), prThreshold);
        
        //Extract set of correctly predicted functions
        Set<Integer> correctlyPredictedFuncs = new TreeSet<>();
        
        for (int og : og2known_funcs.keySet())
        {
            if (og2predicted_funcs_First.containsKey(og))
            {
                og2predicted_funcs_First.get(og).retainAll(og2known_funcs.get(og));
                correctlyPredictedFuncs.addAll(og2predicted_funcs_First.get(og));
            }
            
            if (og2predicted_funcs_Second.containsKey(og))
            {
                og2predicted_funcs_Second.get(og).retainAll(og2known_funcs.get(og));
                correctlyPredictedFuncs.addAll(og2predicted_funcs_Second.get(og));
            }
        }
        
        List<String> methodCombinations = new ArrayList<>();
        methodCombinations.add(firstMethodName);
        methodCombinations.add(firstMethodName + "+" + secondMethodName);
        methodCombinations.add(secondMethodName);
        
        //Count overlaps between methods in correctly predicted functions
        Table<Integer, String, Integer> funcCombCount = ArrayTable.create(correctlyPredictedFuncs, methodCombinations);
        
        for (int func : funcCombCount.rowKeySet())
            for (String comb : funcCombCount.columnKeySet())
                funcCombCount.put(func, comb, 0);
        
        Set<Integer> commonOgs = new HashSet<>();
        commonOgs.addAll(og2predicted_funcs_First.keySet());
        commonOgs.addAll(og2predicted_funcs_Second.keySet());
        
        for (int og : commonOgs)
        {
            Set<Integer> funcs = new HashSet<>();
            
            Set<Integer> firstFuncs = new HashSet<>();
            Set<Integer> secondFuncs = new HashSet<>();
            
            if (og2predicted_funcs_First.containsKey(og))
                firstFuncs = og2predicted_funcs_First.get(og);
            
            if (og2predicted_funcs_Second.containsKey(og))
                secondFuncs = og2predicted_funcs_Second.get(og);
            
            funcs.addAll(firstFuncs);
            funcs.addAll(secondFuncs);
            
            for (int func : funcs)
                if (firstFuncs.contains(func) && !secondFuncs.contains(func))
                {
                    int cnt = funcCombCount.get(func, firstMethodName);
                    cnt++;
                    funcCombCount.put(func, firstMethodName, cnt);
                }
                else if (!firstFuncs.contains(func) && secondFuncs.contains(func))
                {
                    int cnt = funcCombCount.get(func, secondMethodName);
                    cnt++;
                    funcCombCount.put(func, secondMethodName, cnt);
                }
                else
                {
                    int cnt = funcCombCount.get(func, firstMethodName + "+" + secondMethodName);
                    cnt++;
                    funcCombCount.put(func, firstMethodName + "+" + secondMethodName, cnt);
                }
        }
        
        BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outFolder + "/Number_of_times_a_classifier_predicted_a_function.txt")));
        
        bw.write("Function");
        for (String comb : funcCombCount.columnKeySet())
            bw.write("\t" + comb);
        bw.write("\n");
        
        for (int func : funcCombCount.rowKeySet())
        {
            bw.write(String.valueOf(func));
            for (String comb : funcCombCount.columnKeySet())
                bw.write("\t" + funcCombCount.get(func, comb));
            bw.write("\n");
        }
        
        bw.close();
        
        //Create data file for drawing histogram
        String dataFilePath = outFolder.getAbsolutePath() + "/Histogram_data.txt";
        bw = new BufferedWriter(new FileWriter(new File(dataFilePath)));
        
        bw.write("COUNTS\tCATEGORY\n");
        
        Set<Integer> firstMethodSpecificFunctions = new TreeSet<>();
        Set<Integer> secondMethodSpecificFunctions = new TreeSet<>();
        Set<Integer> bothMethodsFunctions = new TreeSet<>();
        
        for (int func : funcCombCount.rowKeySet())
        {
            boolean first = funcCombCount.get(func, firstMethodName) > 0;
            boolean second = funcCombCount.get(func, secondMethodName) > 0;
            boolean both = funcCombCount.get(func, firstMethodName + "+" + secondMethodName) > 0;
            
            if (first && !second && !both)
            {
                bw.write("1\t1\n"); //First
                firstMethodSpecificFunctions.add(func);
            }
            else if (!first && second && !both)
            {
                bw.write("3\t1\n"); //Second
                secondMethodSpecificFunctions.add(func);
            }
            else
            {
                bw.write("2\t1\n"); //Both
                bothMethodsFunctions.add(func);
            }
        }
        
        bw.close();
        
        bw = new BufferedWriter(new FileWriter(new File(outFolder + "/Graph_statistics+Legend.txt")));
        
        int total = firstMethodSpecificFunctions.size() + secondMethodSpecificFunctions.size() + bothMethodsFunctions.size();
        
        bw.write("# functions: " + total + "\n\n");
        bw.write(firstMethodName + ": " + firstMethodSpecificFunctions.size() + " >> " + round((double)firstMethodSpecificFunctions.size()/(double)total*(double)100,2) + "%\n");
        bw.write(secondMethodName + ": " + secondMethodSpecificFunctions.size() + " >> " + round((double)secondMethodSpecificFunctions.size()/(double)total*(double)100,2) + "%\n");
        bw.write(firstMethodName + "+" + secondMethodName + ": " + bothMethodsFunctions.size() + " >> " + round((double)bothMethodsFunctions.size()/(double)total*(double)100,2) + "%\n");
        
        bw.write("\nLEGEND:\n");
        String[] color = colors.split(",");
        bw.write(color[0].replace("'", "") + ": " + firstMethodName + "\n");
        bw.write(color[1].replace("'", "") + ": " + firstMethodName + " and " + secondMethodName + " overlap\n");
        bw.write(color[2].replace("'", "") + ": " + secondMethodName + "\n");
        
        bw.close();
        
        //Write to file lists of GO functions that are predicted only by a single method and by both methods
        bw = new BufferedWriter(new FileWriter(new File(outFolder + "/Functions_predicted_by_" + firstMethodName + ".txt")));
        for (int func : firstMethodSpecificFunctions)
            bw.write(func + "\n");
        bw.close();
        
        bw = new BufferedWriter(new FileWriter(new File(outFolder + "/Functions_predicted_by_" + secondMethodName + ".txt")));
        for (int func : secondMethodSpecificFunctions)
            bw.write(func + "\n");
        bw.close();
        
        bw = new BufferedWriter(new FileWriter(new File(outFolder + "/Functions_predicted_by_" + firstMethodName + "_and_" + secondMethodName + ".txt")));
        for (int func : bothMethodsFunctions)
            bw.write(func + "\n");
        bw.close();
        
        //Generate and run R script
        String rScriptFilePath = outFolder.getAbsolutePath() + "/Rscript.r";
        bw = new BufferedWriter(new FileWriter(rScriptFilePath));
        
        Double prThresholdPerc = prThreshold * 100;
        
        bw.write("jpeg('" + (outFolder.getAbsolutePath() + "/Histogram.jpg").replace('\\', '/') + "',width=500,height=180)\n");  
        bw.write("data = read.table('" + dataFilePath.replaceAll("\\\\", "/") + "', header = T, sep = '\t', na.strings = '?')\n");
        bw.write("counts <- table(data$COUNTS, data$CATEGORY)\n");
        bw.write("barplot(counts,horiz=T,col=c(" + colors + "),cex.axis=1.8,cex.names=1.8,cex.lab=1.4,names.arg = c(\"Pr \\u2265 " +
                 prThresholdPerc.intValue() + "%\")," + "xlab=\"GO functions predicted by: " + firstMethodName + ", both, " +
                 secondMethodName + "\")\n");
        
        bw.write("dev.off()\n");
                
        bw.close();
        
        executeRScript(rScriptFilePath);
        

        new File(dataFilePath).delete();
        new File(rScriptFilePath).delete();
    }
    
    
    /**
     * Draws graph that shows the level of complementarity between two classification models in terms of
     * gene families to which they are able to assign GO functions at a specific level of precision (Pr).
     * This code assumes that you have R on your computer and R package ‘venneuler’.
     * 
     * @param inFirstMethodPrecisionFile            The table with Pr scores outputted by the first classification model. 
     * @param inSecondMethodPrecisionFile           The table with Pr scored outputted by the second classification model.
     * @param inFileWithListOfConsideredFunctions   Supply file with a list of functions if you want to observe complementarity
     *                                              patterns only for the specific subset of functions. Set to null to consider
     *                                              all available functions.
     * @param firstMethodName                       The name of the first method, i.e., classification model.
     * @param secondMethodName                      The name of the second method.
     * @param colorCodes                            Color scheme, composed of two colors separated by comma. Accepts color
     *                                              definitions typical for R package 'venneuler': .X, where X is a number.
     * @param colorNames                            Two color names separated by comma that will be used to generate graph's legend.
     * @param prThreshold                           Predictions with Pr >= prThreshold will be considered as positive.
     * @param inOg2FunctionsFile                    File with known functions from Uniprot-GOA that are assigned to gene families.
     * @param inGeneOntologyFile                    Gene ontology in obo-xml.gz format.
     * @param outFolder                             Folder that will contain graph and statistics files.
     * 
     * @throws IOException
     * @throws InterruptedException 
     */
    public static void drawGeneFamilyBasedComplementarityGraphForTwoMethods(File inFirstMethodPrecisionFile,
                                                                            File inSecondMethodPrecisionFile,
                                                                            File inFileWithListOfConsideredFunctions,
                                                                            String firstMethodName, String secondMethodName,
                                                                            String colorCodes, String colorNames,
                                                                            double prThreshold, File inOg2FunctionsFile,
                                                                            File inGeneOntologyFile, File outFolder) throws IOException, InterruptedException
    {
        if (!outFolder.exists())
            outFolder.mkdirs();
        
        //Load known functions
        Map<Integer, Set<Integer>> og2known_funcs = loadOg2FunctionsFromFile(inOg2FunctionsFile);
        
         //Consider only GO functions from prokaryotic GO subset gosubset_prok
        Set<Integer> knownProkFuncs = new TreeSet<>();
        for (Set<Integer> known_funcs : og2known_funcs.values())
            knownProkFuncs.addAll(known_funcs);
        
        knownProkFuncs = extractOnlyProkaryoticGOs(knownProkFuncs, inGeneOntologyFile);
        
        //Load a subset of relevant functions when supplied, and consider only those prokaryotic functions that
        //overlap with the relevant functions
        if (inFileWithListOfConsideredFunctions != null)
        {
            Set<Integer> relevantFuncs = new HashSet<>();
            
            BufferedReader br = findReaderType(inFileWithListOfConsideredFunctions);
            
            String line;
            
            while((line=br.readLine()) != null)
                relevantFuncs.add(Integer.parseInt(line));
            
            knownProkFuncs.retainAll(relevantFuncs);
        }
        
        //Load predicted functions
        Map<Integer, Set<Integer>> og2predicted_funcs_First = extractOg2FuncsFromPrTable(
                loadPrScoresTable(inFirstMethodPrecisionFile, knownProkFuncs, og2known_funcs.keySet()), prThreshold);
        
        Map<Integer, Set<Integer>> og2predicted_funcs_Second = extractOg2FuncsFromPrTable(
                loadPrScoresTable(inSecondMethodPrecisionFile, knownProkFuncs, og2known_funcs.keySet()), prThreshold);
        
        //Count predicted GO functions per gene families
        Set<Integer> commonOgs = new HashSet<>();
        commonOgs.addAll(og2predicted_funcs_First.keySet());
        commonOgs.addAll(og2predicted_funcs_Second.keySet());
        
        List<String> methodCombinations = new ArrayList<>();
        methodCombinations.add(firstMethodName);
        methodCombinations.add(firstMethodName + "+" + secondMethodName);
        methodCombinations.add(secondMethodName);
        
        Table<Integer, String, Integer> ogCombCount = ArrayTable.create(commonOgs, methodCombinations);
        
        for (int og : ogCombCount.rowKeySet())
            for (String comb : ogCombCount.columnKeySet())
                ogCombCount.put(og, comb, 0);
        
        for (int og : commonOgs)
        {
            Set<Integer> predsFirst = new HashSet<>();
            if (og2predicted_funcs_First.containsKey(og))
                predsFirst = og2predicted_funcs_First.get(og);
            
            Set<Integer> predsSecond = new HashSet<>();
            if (og2predicted_funcs_Second.containsKey(og))
                predsSecond = og2predicted_funcs_Second.get(og);
            
            Set<Integer> known = og2known_funcs.get(og);
            
            predsFirst.retainAll(known);
            predsSecond.retainAll(known);
            
            Set<Integer> overlap = new HashSet<>(predsFirst);
            overlap.retainAll(predsSecond);
            
            Set<Integer> firstOnly = new HashSet<>(predsFirst);
            firstOnly.removeAll(predsSecond);
            
            Set<Integer> secondOnly = new HashSet<>(predsSecond);
            secondOnly.removeAll(predsFirst);
            
            ogCombCount.put(og, firstMethodName, firstOnly.size());
            ogCombCount.put(og, firstMethodName + "+" + secondMethodName, overlap.size());
            ogCombCount.put(og, secondMethodName, secondOnly.size());
        }
        
        //Compute ratio of predictions
        int firstOnlyCount = 0;
        int secondOnlyCount = 0;
        int overlapCount = 0;
        
        for (int og : ogCombCount.rowKeySet())
        {
            firstOnlyCount += ogCombCount.get(og, firstMethodName);
            secondOnlyCount += ogCombCount.get(og, secondMethodName);
            overlapCount += ogCombCount.get(og, firstMethodName + "+" + secondMethodName);
        }
        
        //Write data file with counts
        FileWriter fw = new FileWriter(new File(outFolder + "/Number_of_predicted_functions.txt"));
        BufferedWriter bw = new BufferedWriter(fw);
        
        bw.write("Gene family\t" + firstMethodName + "\t" + firstMethodName + "+" + secondMethodName + "\t" + secondMethodName + "\n");
        
        int cntOGsWithPredictions = 0;
        
        for (int og : ogCombCount.rowKeySet())
            if ((ogCombCount.get(og, firstMethodName) + ogCombCount.get(og, firstMethodName + "+" + secondMethodName) +
                ogCombCount.get(og, secondMethodName)) > 0)
            {
                bw.write(ogToStr(og) + "\t" + ogCombCount.get(og, firstMethodName) + "\t" +
                         ogCombCount.get(og, firstMethodName + "+" + secondMethodName) + "\t" +
                         ogCombCount.get(og, secondMethodName) + "\n");
                
                cntOGsWithPredictions++;
            }
        
        bw.close();
        
        fw = new FileWriter(new File(outFolder + "/Graph_statistics+Legend.txt"));
        bw = new BufferedWriter(fw);
        
        int total = firstOnlyCount + overlapCount + secondOnlyCount;
        
        bw.write(cntOGsWithPredictions + " gene families received " + total + " predictions assigned by:\n\n");
        bw.write(firstMethodName + ": " + firstOnlyCount + " >> " + round((double)firstOnlyCount/(double)total*(double)100,2) + "%\n");
        bw.write(secondMethodName + ": " + secondOnlyCount + " >> " + round((double)secondOnlyCount/(double)total*(double)100,2) + "%\n");
        bw.write(firstMethodName + "+" + secondMethodName + ": " + overlapCount + " >> " + round((double)overlapCount/(double)total*(double)100,2) + "%\n");
        
        bw.write("\nLEGEND:\n");
        String[] color = colorNames.split(",");
        bw.write(color[0].trim() + ": " + firstMethodName + "\n");
        bw.write(color[1].trim() + ": " + secondMethodName + "\n");
        
        bw.close();
        
        //Write R script
        String rScriptFilePath = outFolder.getAbsolutePath() + "/Rscript.r";
        fw = new FileWriter(new File(rScriptFilePath));
        bw = new BufferedWriter(fw);
        
        bw.write("library(venneuler)\n");
        bw.write("pdf('" + outFolder.getAbsolutePath().replace("\\", "/") + "/Venn_diagram.pdf')\n");
        bw.write("vd <- venneuler(c('" + secondMethodName + "'=" + secondOnlyCount + ",'" + firstMethodName + "&" + secondMethodName +
                 "'=" + overlapCount + ",'" + firstMethodName + "'=" + firstOnlyCount + "))\n");
        bw.write("vd$labels <- rep(\"\", length(vd$labels))\n");
        bw.write("vd$colors <- c(" + colorCodes + ")\n");
        bw.write("plot(vd)\n");
        bw.write("dev.off()\n");
        
        bw.close();
        
        executeRScript(rScriptFilePath);
        
        new File(rScriptFilePath).delete();
    }
}
