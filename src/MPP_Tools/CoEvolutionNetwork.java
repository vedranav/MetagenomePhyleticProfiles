package MPP_Tools;

import com.google.common.collect.ArrayTable;
import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import static org.apache.commons.math3.util.Precision.round;
import static utils.ArffUtils.selectSubsetOfLabelsFromHMCArff;
import static utils.FileUtils.findReaderType;
import static utils.RUtils.executeRScript;
import weka.core.Instances;
import weka.core.converters.CSVSaver;
import weka.core.converters.ConverterUtils;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Remove;
import weka.filters.unsupervised.instance.RemoveRange;

/**
 *
 * @author Vedrana Vidulin [vedrana.vidulin@ijs.si]
 * 
 */
public class CoEvolutionNetwork
{
    /**
     * Performs Random Forest-based feature selection by keeping the features with positive values of Gini Index.
     * This procedure depends on R and 'randomForest' R package.
     * 
     * @param inMPPArffFile         Input file with metagenome phyletic profiles (MPP) data set in ARFF format.
     * @param inPPArffFile          Input file with phyletic profiles (PP) data set in ARFF format matched with MPP.
     *                              It is assumed that both data sets have the same set of gene families as instances,
     *                              that they are both in CLUS HMC ARFF format (see examples in "data" folder) and that
     *                              the first feature is string type instance ID feature (name of a gene family).
     * @param mppRelatedFunction    GO function that is predicted by MPP better than by PP
     *                              (see file "MPP-I_and_PP-I_AUPRCs.tsv" in "data").
     * @param ppRelatedFunction     GO function that is predicted by PP better than by MPP.
     *                              (see file "MPP-I_and_PP-I_AUPRCs.tsv" in "data").
     * @param outFolder             MPP and PP data sets with the selected features will be outputted to this folder.
     * 
     * @throws Exception 
     */
    public static void randomForestFeatureSelection(File inMPPArffFile, File inPPArffFile,
                                                    int mppRelatedFunction, int ppRelatedFunction,
                                                    File outFolder) throws Exception
    {
        if (!outFolder.exists())
            outFolder.mkdirs();
        
        //----------------------------------------------------------------
        //SELECT GENE FAMILIES ANNOTATED WITH THE mppRelatedFunction AND
        //ppRelatedFunction GO FUNCTIONS
        //----------------------------------------------------------------
        Set<String> selectedOGs = new TreeSet<>();
        
        BufferedReader br = findReaderType(inMPPArffFile); //becuase MPP and PP data sets are paired, it is enough to read gene families from one of them
        
        String line;
                
        while((line=br.readLine()) != null)
            if (!line.startsWith("@") && !line.isEmpty())
            {
                String og = line.substring(0, line.indexOf(","));
                String funcsStr = line.substring(line.lastIndexOf(",") + 1).trim();
                
                String[] parts = funcsStr.split("@");
                
                Set<Integer> funcs = new HashSet<>();
                for (String part : parts)
                    funcs.add(Integer.parseInt(part));
                
                if (funcs.contains(mppRelatedFunction) || funcs.contains(ppRelatedFunction))
                    selectedOGs.add(og);
            }
        
        System.out.println(selectedOGs.size() + " gene families are annotated with the selected GO functions: " + mppRelatedFunction + " and " + ppRelatedFunction);
        
        
        //---------------------------
        //PERFORM FEATURE SELECTION
        //---------------------------
        File[] arffFiles = {inMPPArffFile, inPPArffFile};
        
        for (File f : arffFiles)
        {
            System.out.println("Performing feature selection for the data set: " + f.getName());
            
            System.out.println("\tPreparing data for feature selection");
            
            //Create binary class ARFF file for the method-specific selected GO function 
            List<String> selectedFunction = new ArrayList<>();

            if (f == inMPPArffFile)
                selectedFunction.add(String.valueOf(mppRelatedFunction));
            else
                selectedFunction.add(String.valueOf(ppRelatedFunction));

            String singleLabelArffFile = outFolder + "/" + f.getName().substring(0, f.getName().lastIndexOf(".arff")) + ".arff";

            selectSubsetOfLabelsFromHMCArff(selectedFunction, f, new File(singleLabelArffFile));

            
            //Remove ID attribute and save the data set into CSV format
            ConverterUtils.DataSource source = new ConverterUtils.DataSource(singleLabelArffFile);
            Instances data = source.getDataSet();
            
            Remove remove = new Remove();
            remove.setAttributeIndices("1");
            remove.setInvertSelection(false);
            remove.setInputFormat(data);
            data = Filter.useFilter(data, remove);
            
            String goFunc = data.attribute(data.numAttributes() - 1).name();
            
            String singleLabelCsvFile = outFolder + "/" + f.getName().substring(0, f.getName().lastIndexOf(".arff")) + ".csv";

            CSVSaver saver = new CSVSaver();
            saver.setInstances(data);
            saver.setFile(new File(singleLabelCsvFile));
            saver.writeBatch();
            

            //Create R script that will compute feature importances
            System.out.println("\tGenerating R script");
            
            String singleLabelRScriptFile = outFolder + "/" + f.getName().substring(0, f.getName().lastIndexOf(".arff")) + ".r";

            String singleLabelImportanceMatrixFile = outFolder + "/" + f.getName().substring(0, f.getName().lastIndexOf(".arff")) + "-FIs.txt";
            
            
            BufferedWriter bwS = new BufferedWriter(new FileWriter(singleLabelRScriptFile));

            bwS.write("library(randomForest)\n");
            bwS.write("data <- read.csv(\"" + singleLabelCsvFile.replace("\\", "/") + "\")\n");

            if (f == inMPPArffFile)
                bwS.write("data$X" + mppRelatedFunction + " <- factor(data$X" + mppRelatedFunction + ")\n");
            else
            {
                bwS.write("cols <- colnames(data)\n");
                bwS.write("data[cols] <- lapply(data[cols], factor)\n");
            }

            bwS.write("set.seed(1)\n");
            bwS.write("data.rf <- randomForest(data$X" + goFunc + " ~ ., data=data, ntree=200, importance=TRUE)\n");
            bwS.write("imp <- round(importance(data.rf), 2)\n");
            bwS.write("impdf <- as.data.frame(imp)\n");
            bwS.write("write.table(impdf, file=\"" + singleLabelImportanceMatrixFile.replace("\\", "/") +
                      "\", sep = \"\\t\")\n");

            bwS.close();

            System.out.println("\tExecuting R script");
            executeRScript(singleLabelRScriptFile);
            
            
            //Extract indices of the selected features
            List<Integer> indicesOfTheSelectedFeatures = new ArrayList<>();
            indicesOfTheSelectedFeatures.add(0); //keep ID

            br = findReaderType(new File(singleLabelImportanceMatrixFile));

            int featureIndex = 1;

            while((line=br.readLine()) != null)
            {
                if (line.startsWith("\"0\"\t\"1\"\t\"MeanDecreaseAccuracy\"\t\"MeanDecreaseGini\"")) //skip header
                    continue;

                String[] parts = line.split("\t");

                if (Double.parseDouble(parts[2]) > 0)
                    indicesOfTheSelectedFeatures.add(featureIndex);
                
                featureIndex++;
            }
            
            System.out.println((indicesOfTheSelectedFeatures.size() - 1) + " features selected");
            
            
            //Create data set with the selected features
            source = new ConverterUtils.DataSource(singleLabelArffFile);
            data = source.getDataSet();
            
            indicesOfTheSelectedFeatures.add(data.numAttributes() - 1);
            
            int[] indicesOfTheSelectedFeaturesArray = new int[indicesOfTheSelectedFeatures.size()];
            for (int i = 0; i < indicesOfTheSelectedFeatures.size(); i++)
                indicesOfTheSelectedFeaturesArray[i] = indicesOfTheSelectedFeatures.get(i);

            remove = new Remove();
            remove.setAttributeIndicesArray(indicesOfTheSelectedFeaturesArray);
            remove.setInvertSelection(true);
            remove.setInputFormat(data);
            data = Filter.useFilter(data, remove);
            
            
            //Keep only gene families that are annotated with the selected GO functions 
            String indicesOfOGsThatAreIrrelevant = "";
            for (int i = 0; i < data.numInstances(); i++)
            {
                String og = data.instance(i).stringValue(0);
                
                if (!selectedOGs.contains(og))
                    indicesOfOGsThatAreIrrelevant += "," + (i + 1);
            }
            
            indicesOfOGsThatAreIrrelevant = indicesOfOGsThatAreIrrelevant.substring(1);
            
            RemoveRange removeInstances = new RemoveRange();
            removeInstances.setInstancesIndices(indicesOfOGsThatAreIrrelevant);
            removeInstances.setInvertSelection(false);
            removeInstances.setInputFormat(data);
            data = Filter.useFilter(data, removeInstances);
            
            saver = new CSVSaver();
            saver.setInstances(data);
            saver.setFile(new File(singleLabelCsvFile));
            saver.writeBatch();
            
            
            //Remove reduntant files
            new File(singleLabelArffFile).delete();
            new File(singleLabelRScriptFile).delete();
            new File(singleLabelImportanceMatrixFile).delete();
        }
    }
    
    /**
     * Computes similarities between gene family profiles. Similarity is measured as Pearson correlation coefficient.
     * Please note that the direction of correlation is irrelevant in this application.
     * Consequently, we considered only the strength of correlation. 
     * 
     * @param inMPPDatasetWithSelectedFeatures  The MPP data set with selected features outputted by "randomForestFeatureSelection" procedure.
     * @param inPPDatasetWithSelectedFeatures   The PP data set with selected features outputted by "randomForestFeatureSelection" procedure.
     * @param outFileWithPearsonCorrelationCoefficientsForOGPairs   Output file that will contain Pearson correlation coefficients.
     * 
     * @throws IOException 
     */
    public static void computeSimilaritiesBetweenOGs(File inMPPDatasetWithSelectedFeatures, File inPPDatasetWithSelectedFeatures,
                                                     File outFileWithPearsonCorrelationCoefficientsForOGPairs) throws IOException
    {
        //Extract the selected gene families
        Set<String> selectedOGs = new TreeSet<>();
        
        BufferedReader br = findReaderType(inMPPDatasetWithSelectedFeatures);
        
        String line;
        
        boolean header = true;
                
        while((line=br.readLine()) != null)
            if (header)
                header = false;
            else
                selectedOGs.add(line.substring(0, line.indexOf(",")));
        
        System.out.println("Pairwise similarity will be computed for " + selectedOGs.size() + " gene families");

        //Construct all possible pairs of gene families
        Set<Set<String>> ogPairs = new HashSet<>();

        for (String ogFirst : selectedOGs)
            for (String ogSecond : selectedOGs)
            {
                Set<String> pair = new TreeSet<>();
                pair.add(ogFirst);
                pair.add(ogSecond);

                if (pair.size() == 2)
                    ogPairs.add(pair);
            }
        
        List<String> methods = new ArrayList<>();
        methods.add("MPP");
        methods.add("PP");
        
        Table<Set<String>, String, Double> ogPair_method_pearsonCC = ArrayTable.create(ogPairs, methods);

        System.out.println("Number of gene family pairs: " + ogPairs.size());
        
        
        File[] csvFiles = {inMPPDatasetWithSelectedFeatures, inPPDatasetWithSelectedFeatures};
        
        for (File f : csvFiles)
        {
            //Load instances
            Map<String, double[]> ogRepresentingVectors = new TreeMap<>();
            
            br = findReaderType(f);
            
            header = true;
            
            while((line=br.readLine()) != null)
                if (header)
                    header = false;
                else
                {
                    String parts[] = line.substring(line.indexOf(",") + 1).trim().split(",");
                    
                    double[] vals = new double[parts.length - 1];
                    for (int i = 0; i < parts.length - 1; i++)
                        vals[i] = Double.parseDouble(parts[i]);
                    
                    ogRepresentingVectors.put(line.substring(0, line.indexOf(",")), vals);
                }
            
            //Compute similarities between gene families in MPP and PP
            for (Set<String> pair : ogPairs)
            {
                List<String> pairList = new ArrayList<>(pair);
                
                double pearsonCC = Math.abs(round(
                    new PearsonsCorrelation().correlation(
                        ogRepresentingVectors.get(pairList.get(0)),
                        ogRepresentingVectors.get(pairList.get(1))), 4));

                ogPair_method_pearsonCC.put(pair, (f.equals(inMPPDatasetWithSelectedFeatures) ? "MPP" : "PP"), pearsonCC);
            }
        }
    
        //Save computed correlation coefficients
        BufferedWriter bw = new BufferedWriter(new FileWriter(outFileWithPearsonCorrelationCoefficientsForOGPairs));
        
        //Header
        bw.write("Gene family pair");
        for (String method : ogPair_method_pearsonCC.columnKeySet())
            bw.write("\t" + method);
        bw.write("\n");
        
        for (Set<String> ogPair : ogPair_method_pearsonCC.rowKeySet())
        {
            List<String> pairList = new ArrayList<>(ogPair);
            
            bw.write(pairList.get(0) + "-" + pairList.get(1));
            for (String method : ogPair_method_pearsonCC.columnKeySet())
                bw.write("\t" + ogPair_method_pearsonCC.get(ogPair, method));
            bw.write("\n");
        }
        
        bw.close();
    }
    
    /**
     * Draws histograms with distributions of Pearson corelation coefficients for MPP and PP.
     * Histograms can help to determine the threshold that will be used to filter out less relevant edges in a graph.
     * This procedure depends on R.
     * 
     * @param inFileWithPCC File that contains Pearson correlation coefficients outputted by "computeSimilaritiesBetweenOGs".
     * @param outFolder     Histograms in PNG format for MPP and PP will be outputted to this folder.
     * 
     * @throws IOException
     * @throws InterruptedException 
     */
    public static void drawPearsonCorrelationCoefficientDistributionGraph(File inFileWithPCC, File outFolder) throws IOException, InterruptedException
    {
        if (!outFolder.exists())
            outFolder.mkdirs();
        
        //Load Pearson correlation coefficients from file
        Table<String, String, Double> ogPair_method_PCC = TreeBasedTable.create();
        
        List<String> methods = new ArrayList<>();
        
        BufferedReader br = findReaderType(inFileWithPCC);
        
        String line;
        
        boolean header = true;
                
        while((line=br.readLine()) != null)
            if (header)
            {
                header = false;
                
                String[] parts = line.split("\t");
                
                for (int i = 1; i < parts.length; i++)
                    methods.add(parts[i]);
            }
            else
            {
                String[] parts = line.split("\t");
                
                for (int i = 1; i < parts.length; i++)
                    ogPair_method_PCC.put(parts[0], methods.get(i - 1), Double.parseDouble(parts[i]));
            }
        
        
        //Draw histograms with distributions of Pearson correlation coefficients for MPP and PP
        for (String method : methods)
        {
            String histogramDataFilePath = outFolder + "/" + method + "-histogram.csv";
            
            BufferedWriter bw = new BufferedWriter(new FileWriter(histogramDataFilePath));
            
            for (String ogPair : ogPair_method_PCC.rowKeySet())
                bw.write(String.valueOf(ogPair_method_PCC.get(ogPair, method)) + "\n");
            
            bw.close();
            
            String histogramRScriptFilePath = outFolder + "/" + method + "-histogram.r";

            bw = new BufferedWriter(new FileWriter(histogramRScriptFilePath));

            bw.write("png('" + histogramDataFilePath.substring(0, histogramRScriptFilePath.lastIndexOf(".")).replace("\\", "/") + ".png')\n");
            bw.write("data <- read.csv('" + histogramDataFilePath.replace("\\", "/") + "', header=F)\n");
            bw.write("histogram <- hist(data$V1, plot=F)\n");
            bw.write("plot(histogram, ylim=c(0, max(histogram$counts) + 5), main='', xlab='Pearson correlation coefficient')\n");
            bw.write("text(histogram$mids, histogram$counts + 60, histogram$counts, cex=0.75)\n");
            bw.write("dev.off()\n");

            bw.close();

            executeRScript(histogramRScriptFilePath);
            
            //Remove reduntant files
            new File(histogramDataFilePath).delete();
            new File(histogramRScriptFilePath).delete();
        }
    }
    
    /**
     * Generates gexf file with the description of coevolution network, which can be visualized in Gephi (http://gephi.org).
     * 
     * @param inMPPDatasetWithSelectedFeatures  The MPP data set with selected features outputted by "randomForestFeatureSelection" procedure.
     * @param inPPDatasetWithSelectedFeatures   The PP data set with selected features outputted by "randomForestFeatureSelection" procedure.
     * @param pccThreshold  The network will include only edges with the absolute value of Pearson correlation coefficient > pccThreshold
     *                      and interconnected nodes.
     * @param inFileWithPCC File that contains Pearson correlation coefficients outputted by "computeSimilaritiesBetweenOGs".
     * @param outGephiFile  Gexf file with the description of network.
     * 
     * @throws IOException 
     */
    public static void composeGephiFileWithNetwork(File inMPPDatasetWithSelectedFeatures, File inPPDatasetWithSelectedFeatures,
                                                   double pccThreshold, File inFileWithPCC, File outGephiFile) throws IOException
    {
        String mppSpecificFunction = "";
        String ppSpecificFunction = "";
        
        //---------------------
        //LOAD NODES AND EDGES
        //---------------------
        //Load gene families and determine whether they are associated with MPP- or PP-specific GO function
        Map<String, String> ogMethod = new HashMap<>();
        
        BufferedReader br = findReaderType(inMPPDatasetWithSelectedFeatures);
        
        String line;
        
        boolean header = true;
        
        while((line=br.readLine()) != null)
            if (header)
            {
                mppSpecificFunction = line.substring(line.lastIndexOf(",") + 1);
                header = false;
            }
            else if (line.substring(line.lastIndexOf(",") + 1).equals("1"))
               ogMethod.put(line.substring(0, line.indexOf(",")), "MPP");
        
        
        br = findReaderType(inPPDatasetWithSelectedFeatures);
        
        header = true;
        
        while((line=br.readLine()) != null)
            if (header)
            {
                ppSpecificFunction = line.substring(line.lastIndexOf(",") + 1);
                header = false;
            }
            else if (line.substring(line.lastIndexOf(",") + 1).equals("1"))
            {
                String og = line.substring(0, line.indexOf(","));
                
                if (ogMethod.containsKey(og))
                    ogMethod.put(og, "both");
                else
                    ogMethod.put(og, "PP");
            }
        
        
        //Load correlations that are > threshold
        Map<String, Double> ogPairPCC_MPP = new HashMap<>();
        Map<String, Double> ogPairPCC_PP = new HashMap<>();
        
        br = findReaderType(inFileWithPCC);
        
        header = true;
        
        while((line=br.readLine()) != null)
            if (header)
                header = false;
            else
            {
                String[] parts = line.split("\t");
                
                Double pccMPP = Double.parseDouble(parts[1]);
                Double pccPP = Double.parseDouble(parts[2]);
                
                if (pccMPP > pccThreshold)
                    ogPairPCC_MPP.put(parts[0], pccMPP);
                
                if (pccPP > pccThreshold)
                    ogPairPCC_PP.put(parts[0], pccPP);
            }
        
        
        //--------------------------------
        //COMPUTE WEIGHTs OF SHARED EDGES
        //--------------------------------
        Map<String, Double> ogPairPCC_MPP_PP = new HashMap<>();
        
        Set<String> sharedPairs = new HashSet<>(ogPairPCC_MPP.keySet());
        sharedPairs.retainAll(ogPairPCC_PP.keySet());

        for (String sharedPair : sharedPairs)
        {
            double mppCC = ogPairPCC_MPP.get(sharedPair);
            double ppCC = ogPairPCC_PP.get(sharedPair);
            
            ogPairPCC_MPP_PP.put(sharedPair, round((mppCC + ppCC) / (double)2, 4));
            
            ogPairPCC_MPP.remove(sharedPair);
            ogPairPCC_PP.remove(sharedPair);
        }

        
        //--------------------
        //GENERATE GEPHI FILE
        //--------------------
        FileWriter fw = new FileWriter(outGephiFile);
        BufferedWriter bw = new BufferedWriter(fw);

        bw.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        bw.write("<gexf xmlns=\"http://www.gexf.net/1.3draft\" xmlns:viz=\"http://www.gexf.net/1.1draft/viz\" version=\"1.3\">\n");
        bw.write("\t<graph defaultedgetype=\"undirected\" idtype=\"string\" type=\"static\">\n");

        bw.write("\t\t<nodes>\n");
        
        //Include in graph only interconnected nodes
        Set<String> interconnectedNodes = new TreeSet<>();

        for (String ogPair : ogPairPCC_MPP.keySet())
        {
            interconnectedNodes.add(ogPair.substring(0, ogPair.indexOf("-")));
            interconnectedNodes.add(ogPair.substring(ogPair.indexOf("-") + 1));
        }

        for (String ogPair : ogPairPCC_PP.keySet())
        {
            interconnectedNodes.add(ogPair.substring(0, ogPair.indexOf("-")));
            interconnectedNodes.add(ogPair.substring(ogPair.indexOf("-") + 1));
        }
        
        for (String ogPair : ogPairPCC_MPP_PP.keySet())
        {
            interconnectedNodes.add(ogPair.substring(0, ogPair.indexOf("-")));
            interconnectedNodes.add(ogPair.substring(ogPair.indexOf("-") + 1));
        }

        List<String> ogs = new ArrayList<>(interconnectedNodes);

        
        //Enumerate interconnected nodes and define their properties
        Set<String> ogsAnnotatedWith_MPP_SpecificFunction = new HashSet<>();
        Set<String> ogsAnnotatedWith_PP_SpecificFunction = new HashSet<>();
        Set<String> ogsAnnotatedWith_BOTH_Functions = new HashSet<>();
        
        for (int i = 0; i < ogs.size(); i++)
        {
            bw.write("\t\t\t<node id=\"" + i + "\" label=\"" + ogs.get(i) + "\">\n");

            String method = ogMethod.get(ogs.get(i));

            if (method.equals("MPP")) //red
            {
                bw.write("\t\t\t\t<viz:color r=\"255\" g=\"0\" b=\"0\" a=\"0.6\"/>\n");
                ogsAnnotatedWith_MPP_SpecificFunction.add(ogs.get(i));
            }
            else if (method.equals("PP")) //blue
            {
                bw.write("\t\t\t\t<viz:color r=\"0\" g=\"0\" b=\"255\" a=\"0.6\"/>\n");
                ogsAnnotatedWith_PP_SpecificFunction.add(ogs.get(i));
            }
            else if (method.equals("both")) //green
            {
                bw.write("\t\t\t\t<viz:color r=\"0\" g=\"255\" b=\"0\" a=\"0.6\"/>\n");
                ogsAnnotatedWith_BOTH_Functions.add(ogs.get(i));
            }

            bw.write("\t\t\t</node>\n");
        }
        bw.write("\t\t</nodes>\n");

        
        //Enumerate edges and define their properties
        bw.write("\t\t<edges>\n");

        Set<String> MPP_specificFunctionOGsConnectedIn_MPP_Network = new TreeSet<>();
        Set<String> MPP_specificFunctionOGsConnectedIn_PP_Network = new TreeSet<>();
        Set<String> PP_specificFunctionOGsConnectedIn_MPP_Network = new TreeSet<>();
        Set<String> PP_specificFunctionOGsConnectedIn_PP_Network = new TreeSet<>();

        int cnt = 0;
        
        //MPP layer
        for (String ogPair : ogPairPCC_MPP.keySet())
        {
            String firstOG = ogPair.substring(0, ogPair.indexOf("-"));
            String secondOG = ogPair.substring(ogPair.indexOf("-") + 1);
            
            bw.write("\t\t\t<edge id=\"" + cnt + "\" source=\"" +
                     ogs.indexOf(firstOG) + "\" target=\"" +
                     ogs.indexOf(secondOG) + "\" weight=\"" +
                     ogPairPCC_MPP.get(ogPair) + "\">\n");
            bw.write("\t\t\t\t<viz:color r=\"255\" g=\"0\" b=\"0\" a=\"1\"/>\n");
            bw.write("\t\t\t</edge>\n");
            
            cnt++;

            if (ogsAnnotatedWith_MPP_SpecificFunction.contains(firstOG) || ogsAnnotatedWith_BOTH_Functions.contains(firstOG))
                MPP_specificFunctionOGsConnectedIn_MPP_Network.add(firstOG);
            
            if (ogsAnnotatedWith_MPP_SpecificFunction.contains(secondOG) || ogsAnnotatedWith_BOTH_Functions.contains(secondOG))
                MPP_specificFunctionOGsConnectedIn_MPP_Network.add(secondOG);

            if (ogsAnnotatedWith_PP_SpecificFunction.contains(firstOG) || ogsAnnotatedWith_BOTH_Functions.contains(firstOG))
                PP_specificFunctionOGsConnectedIn_MPP_Network.add(firstOG);
            
            if (ogsAnnotatedWith_PP_SpecificFunction.contains(secondOG) || ogsAnnotatedWith_BOTH_Functions.contains(secondOG))
                PP_specificFunctionOGsConnectedIn_MPP_Network.add(secondOG);
        }

        //PP layer
        for (String ogPair : ogPairPCC_PP.keySet())
        {
            String firstOG = ogPair.substring(0, ogPair.indexOf("-"));
            String secondOG = ogPair.substring(ogPair.indexOf("-") + 1);

            bw.write("\t\t\t<edge id=\"" + cnt + "\" source=\"" +
                     ogs.indexOf(firstOG) + "\" target=\"" +
                     ogs.indexOf(secondOG) + "\" weight=\"" +
                     ogPairPCC_PP.get(ogPair) + "\">\n");
            bw.write("\t\t\t\t<viz:color r=\"0\" g=\"0\" b=\"255\" a=\"1\"/>\n");
            bw.write("\t\t\t</edge>\n");
            
            cnt++;

            if (ogsAnnotatedWith_MPP_SpecificFunction.contains(firstOG) || ogsAnnotatedWith_BOTH_Functions.contains(firstOG))
                MPP_specificFunctionOGsConnectedIn_PP_Network.add(firstOG);
            
            if (ogsAnnotatedWith_MPP_SpecificFunction.contains(secondOG) || ogsAnnotatedWith_BOTH_Functions.contains(secondOG))
                MPP_specificFunctionOGsConnectedIn_PP_Network.add(secondOG);

            if (ogsAnnotatedWith_PP_SpecificFunction.contains(firstOG) || ogsAnnotatedWith_BOTH_Functions.contains(firstOG))
                PP_specificFunctionOGsConnectedIn_PP_Network.add(firstOG);
            
            if (ogsAnnotatedWith_PP_SpecificFunction.contains(secondOG) || ogsAnnotatedWith_BOTH_Functions.contains(secondOG))
                PP_specificFunctionOGsConnectedIn_PP_Network.add(secondOG);
        }

        //Connections in both layers
        for (String ogPair : ogPairPCC_MPP_PP.keySet())
        {
            String firstOG = ogPair.substring(0, ogPair.indexOf("-"));
            String secondOG = ogPair.substring(ogPair.indexOf("-") + 1);

            bw.write("\t\t\t<edge id=\"" + cnt + "\" source=\"" +
                     ogs.indexOf(firstOG) + "\" target=\"" +
                     ogs.indexOf(secondOG) + "\" weight=\"" +
                     ogPairPCC_MPP_PP.get(ogPair) + "\">\n");
            bw.write("\t\t\t\t<viz:color r=\"0\" g=\"255\" b=\"0\" a=\"1\"/>\n");
            bw.write("\t\t\t</edge>\n");
            
            cnt++;
            
            if (ogsAnnotatedWith_MPP_SpecificFunction.contains(firstOG) || ogsAnnotatedWith_BOTH_Functions.contains(firstOG))
            {
                MPP_specificFunctionOGsConnectedIn_MPP_Network.add(firstOG);
                MPP_specificFunctionOGsConnectedIn_PP_Network.add(firstOG);
            }
            
            if (ogsAnnotatedWith_MPP_SpecificFunction.contains(secondOG) || ogsAnnotatedWith_BOTH_Functions.contains(secondOG))
            {
                MPP_specificFunctionOGsConnectedIn_MPP_Network.add(secondOG);
                MPP_specificFunctionOGsConnectedIn_PP_Network.add(secondOG);
            }

            if (ogsAnnotatedWith_PP_SpecificFunction.contains(firstOG) || ogsAnnotatedWith_BOTH_Functions.contains(firstOG))
            {
                PP_specificFunctionOGsConnectedIn_MPP_Network.add(firstOG);
                PP_specificFunctionOGsConnectedIn_PP_Network.add(firstOG);
            }
            
            if (ogsAnnotatedWith_PP_SpecificFunction.contains(secondOG) || ogsAnnotatedWith_BOTH_Functions.contains(secondOG))
            {
                PP_specificFunctionOGsConnectedIn_MPP_Network.add(secondOG);
                PP_specificFunctionOGsConnectedIn_PP_Network.add(secondOG);
            }
        }

        bw.write("\t\t</edges>\n");

        bw.write("\t</graph>\n");

        bw.write("</gexf>\n");

        bw.close();

        //Statistics
        System.out.println("-------------------");
        System.out.println("NETWORK STATISTICS");
        System.out.println("-------------------");
        System.out.println(interconnectedNodes.size() + " gene families have similar profiles in MPP and/or PP, and are included in the network, of which:");
        
        System.out.println("\t" + (ogsAnnotatedWith_MPP_SpecificFunction.size() + ogsAnnotatedWith_BOTH_Functions.size()) + " are annotated in Uniprot-GOA with a GO function, which MPP predicted better than PP (MPP-specific)");
        System.out.println("\t" + (ogsAnnotatedWith_PP_SpecificFunction.size() + ogsAnnotatedWith_BOTH_Functions.size()) + " are annotated in Uniprot-GOA with a GO function, which PP predicted better than MPP (PP-specific)");
        System.out.println("\t" + ogsAnnotatedWith_BOTH_Functions.size() + " of the " + interconnectedNodes.size() + " gene families are annotated in Uniprot-GOA with both functions");
        
        System.out.println("\n" + MPP_specificFunctionOGsConnectedIn_MPP_Network.size() + " gene families annotated with an MPP-specific GO function are connected in the MPP profiles-based similarity network");
        System.out.println(PP_specificFunctionOGsConnectedIn_MPP_Network.size() + " gene families annotated with a PP-specific GO function are connected in the MPP profiles-based similarity network");
        System.out.println(MPP_specificFunctionOGsConnectedIn_PP_Network.size() + " gene families annotated with an MPP-specific GO function are connected in the PP profiles-based similarity network");
        System.out.println(PP_specificFunctionOGsConnectedIn_PP_Network.size() + " gene falimies annotated with a PP-specific GO function are connected in the PP profiles-based similarity network");
        
        //Instructions on how to draw the network
        System.out.println("\n---------------------------------------------");
        System.out.println("INSTRUCTIONS ON HOW TO VISUALIZE THE NETWORK");
        System.out.println("---------------------------------------------");
        System.out.println("This script outputs a gexf file that can be visualized in Gephi (http://gephi.org).");
        System.out.println("After instaling Gephi, select File > Open and choose the gexf file. On the \"Import report\"");
        System.out.println("window select \"Ok\" and you will see the raw network. The next step is to select a layout.");
        System.out.println("We recommend Fruchterman Reingold, which can be selected in the panel \"Layout\" on the left side.");
        System.out.println("Press the \"Run\" button to apply the layout. Please note that layout algorithms are stohastic");
        System.out.println("and each time will rearrange the network in a sligtly different manner.");
        
        //Network legend
        System.out.println("\n-------");
        System.out.println("LEGEND");
        System.out.println("-------");
        System.out.println("Nodes - Gene families annotated in Uniprot-GOA with a GO function - predicted better by:");
        System.out.println("\tRed: GO:" + mppSpecificFunction + " - MPP");
        System.out.println("\tBlue: GO:" + ppSpecificFunction + " - PP");
        System.out.println("\tGreen: both GOs - both MPP and PP");
        
        System.out.println("\nEdges - Similarities (measured as Pearson correlation coefficient) between gene family profiles in:");
        System.out.println("\tRed: MPP");
        System.out.println("\tBlue: PP");
        System.out.println("\tGreen: both");
        System.out.println("\nWidth of the edges represents the level of similarity, where thicker edges represent higher similarity.");
        System.out.println("The edges represent the absolute values of correlation coefficients > " + pccThreshold);
    }
}
