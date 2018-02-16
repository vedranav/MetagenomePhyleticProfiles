import com.google.common.collect.ArrayTable;
import com.google.common.collect.Table;
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
import static utils.GoFunctionsUtils.extractOnlyProkaryoticGOs;
import static utils.RUtils.executeRScript;
import static utils.ResourceLoaders.extractOg2FuncsFromPrTable;
import static utils.ResourceLoaders.loadOg2FunctionsFromFile;
import static utils.ResourceLoaders.loadPrScoresTable;

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
        boolean Fig1dExample = true;
        boolean Fig1_fgh_A_Example = true;
        boolean Fig1_fgh_B_Example = true;
        boolean Fig1_fgh_C_Example = true;
        //-----------------------------------
        
        //-----------------------------------
        //SET PATHS
        String dataDir = ".../MetagenomePhyleticProfiles/src/data/";
        String outDir = "...";
        //-----------------------------------
        

        double[] prThresholds = {0.5, 0.7, 0.9};
        
        if (Fig1dExample)
            for (double prThreshold : prThresholds)
                drawFunctionBasedComplementarityGraphForTwoMethods(new File(dataDir + "MPP-H_predictions.tsv.gz"),
                    new File(dataDir + "MPP-O_predictions.tsv.gz"), "MPP-H", "MPP-O", "'yellow','green','blue'", prThreshold,
                    new File(dataDir + "og2funcs-eggNOG_3.txt.gz"), new File(dataDir + "go_201401-termdb.obo-xml.gz"),
                    new File(outDir + "/Fig1d/PR-" + prThreshold));
        
        if (Fig1_fgh_A_Example)
            for (double prThreshold : prThresholds)
                drawFunctionBasedComplementarityGraphForTwoMethods(new File(dataDir + "MPP-H_predictions.tsv.gz"),
                    new File(dataDir + "PP-H_predictions.tsv.gz"), "MPP-H", "PP", "'red','green','blue'", prThreshold,
                    new File(dataDir + "og2funcs-eggNOG_3.txt.gz"), new File(dataDir + "go_201401-termdb.obo-xml.gz"),
                    new File(outDir + "/Fig1_fgh_A/PR-" + prThreshold));
        
        if (Fig1_fgh_B_Example)
            for (double prThreshold : prThresholds)
                drawFunctionBasedComplementarityGraphForTwoMethods(new File(dataDir + "MPP-O_predictions.tsv.gz"),
                    new File(dataDir + "PP-O_predictions.tsv.gz"), "MPP-O", "PP", "'red','green','blue'", prThreshold,
                    new File(dataDir + "og2funcs-eggNOG_3.txt.gz"), new File(dataDir + "go_201401-termdb.obo-xml.gz"),
                    new File(outDir + "/Fig1_fgh_B/PR-" + prThreshold));
        
        if (Fig1_fgh_C_Example)
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
        
        FileWriter fw = new FileWriter(new File(outFolder + "/Functions_counts.txt"));
        BufferedWriter bw = new BufferedWriter(fw);
        
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
        fw = new FileWriter(new File(dataFilePath));
        bw = new BufferedWriter(fw);
        
        bw.write("COUNTS\tCATEGORY\n");
        
        int firstCnt = 0;
        int secondCnt = 0;
        int bothCnt = 0;
        
        for (int func : funcCombCount.rowKeySet())
        {
            boolean first = funcCombCount.get(func, firstMethodName) > 0;
            boolean second = funcCombCount.get(func, secondMethodName) > 0;
            boolean both = funcCombCount.get(func, firstMethodName + "+" + secondMethodName) > 0;
            
            if (first && !second && !both)
            {
                bw.write("1\t1\n"); //First
                firstCnt++;
            }
            else if (!first && second && !both)
            {
                bw.write("3\t1\n"); //Second
                secondCnt++;
            }
            else
            {
                bw.write("2\t1\n"); //Both
                bothCnt++;
            }
        }
        
        bw.close();
        
        fw = new FileWriter(new File(outFolder + "/Counts+Legend.txt"));
        bw = new BufferedWriter(fw);
        
        int total = firstCnt + secondCnt + bothCnt;
        
        bw.write("# functions: " + total + "\n\n");
        bw.write(firstMethodName + ": " + firstCnt + " >> " + round((double)firstCnt/(double)total*(double)100,2) + "%\n");
        bw.write(secondMethodName + ": " + secondCnt + " >> " + round((double)secondCnt/(double)total*(double)100,2) + "%\n");
        bw.write(firstMethodName + "+" + secondMethodName + ": " + bothCnt + " >> " + round((double)bothCnt/(double)total*(double)100,2) + "%\n");
        
        bw.write("\nLEGEND:\n");
        String[] color = colors.split(",");
        bw.write(color[0].replace("'", "") + ": " + firstMethodName + "\n");
        bw.write(color[1].replace("'", "") + ": " + firstMethodName + " and " + secondMethodName + " overlap\n");
        bw.write(color[2].replace("'", "") + ": " + secondMethodName + "\n");
        
        bw.close();
        
        //Generate and run R script
        String rScriptFilePath = outFolder.getAbsolutePath() + "/Rscript.r";
        fw = new FileWriter(rScriptFilePath);
        bw = new BufferedWriter(fw);
        
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
}
