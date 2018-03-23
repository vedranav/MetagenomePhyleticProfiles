import static MPP_Tools.CoEvolutionNetwork.composeGephiFileWithNetwork;
import static MPP_Tools.CoEvolutionNetwork.computeSimilaritiesBetweenOGs;
import static MPP_Tools.CoEvolutionNetwork.drawPearsonCorrelationCoefficientDistributionGraph;
import static MPP_Tools.CoEvolutionNetwork.randomForestFeatureSelection;
import static MPP_Tools.ComplementarityGraphs.drawFunctionBasedComplementarityGraphForMultipleModels;
import static MPP_Tools.ComplementarityGraphs.drawFunctionBasedComplementarityGraphForTwoMethods;
import static MPP_Tools.ComplementarityGraphs.drawGeneFamilyBasedComplementarityGraphForTwoMethods;
import static MPP_Tools.DistributionOfAUPRCsBoxPlot.drawDistributionOfAUPRCsBoxPlot;
import java.io.File;

/**
 *
 * @author Vedrana Vidulin [vedrana.vidulin@ijs.si]
 * 
 * If you find this code useful, please cite the paper Vidulin, V., Smuc, T., Dzeroski, S. and Supek, F. (2018)
 * Automated gene function prediction using metagenome data. Under review in Microbiome.
 * 
 */
public class ExperimentsFromThePaper
{
    public static void main(String[] args) throws Exception
    {
        //IMPORTANT!!!!
        //Please begin with README file
        
        //-------------------------------------
        //SELECT THE EXAMPLE(S) FROM THE PAPER
        boolean Fig1a_MPP_H = true;
        boolean Fig1a_MPP_O = true;
        boolean Fig1_d_e = true;
        boolean Fig1_fgh_A = true;
        boolean Fig1_fgh_B = true;
        boolean Fig1_fgh_C = true;
        boolean Fig2_a = true;
        boolean Fig3a = true;
        boolean Fig3c = true;
        boolean Fig4a = true;
        boolean Fig4b = true;
        boolean Fig4c = true;
        //-------------------------------------
        
        //-----------------------------------
        //SET PATHS
        String dataDir = "/MetagenomePhyleticProfiles/src/data/";
        String outDir = "";
        //-----------------------------------
        
        
        double[] prThresholds = {0.5, 0.7, 0.9};
        
        File inGeneOntologyFile = new File(dataDir + "go_201401-termdb.obo-xml.gz");
        File inFileWithFunctionFrequencies = new File(dataDir + "Uniprot-freqs-2071_organisms-Uniprot-GOA-2013-12-10.txt");
        
        if (Fig1a_MPP_H)
        {
            File[] inFilesWithAUPRCsAndPredictionsStats = {new File(dataDir + "MPP-H-AUPRCs+prediction_statistics.txt"),
                                                           new File(dataDir + "MPP-H-Baseline-AUPRCs+prediction_statistics.txt")};
            String[] classifierNames = {"MPP-H", "Baseline"};
            String[] colors = {"red", "black"};
            
            drawDistributionOfAUPRCsBoxPlot(inFilesWithAUPRCsAndPredictionsStats, classifierNames, colors, inGeneOntologyFile,
                                            inFileWithFunctionFrequencies, new File(outDir + "/Fig1a-MPP-H/"));
        }
        
        if (Fig1a_MPP_O)
        {
            File[] inFilesWithAUPRCsAndPredictionsStats = {new File(dataDir + "MPP-O-AUPRCs+prediction_statistics.txt"),
                                                           new File(dataDir + "MPP-O-Baseline-AUPRCs+prediction_statistics.txt")};
            String[] classifierNames = {"MPP-O", "Baseline"};
            String[] colors = {"red", "black"};
            
            drawDistributionOfAUPRCsBoxPlot(inFilesWithAUPRCsAndPredictionsStats, classifierNames, colors, inGeneOntologyFile,
                                            inFileWithFunctionFrequencies, new File(outDir + "/Fig1a-MPP-O/"));
        }
        
        if (Fig1_d_e)
            for (double prThreshold : prThresholds)
            {
                File inFirstMethodPrecisionFile = new File(dataDir + "MPP-H_predictions.tsv.gz");
                File inSecondMethodPrecisionFile = new File(dataDir + "MPP-O_predictions.tsv.gz");
                String firstMethodName = "MPP-H";
                String secondMethodName = "MPP-O";
                File inOg2FunctionsFile = new File(dataDir + "og2funcs-eggNOG_3.txt.gz");
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
                    new File(dataDir + "og2funcs-eggNOG_3.txt.gz"), inGeneOntologyFile,
                    new File(outDir + "/Fig1_fgh_A/PR-" + prThreshold));
        
        if (Fig1_fgh_B)
            for (double prThreshold : prThresholds)
                drawFunctionBasedComplementarityGraphForTwoMethods(new File(dataDir + "MPP-O_predictions.tsv.gz"),
                    new File(dataDir + "PP-O_predictions.tsv.gz"), "MPP-O", "PP", "'red','green','blue'", prThreshold,
                    new File(dataDir + "og2funcs-eggNOG_3.txt.gz"), inGeneOntologyFile,
                    new File(outDir + "/Fig1_fgh_B/PR-" + prThreshold));
        
        if (Fig1_fgh_C)
            for (double prThreshold : prThresholds)
                drawFunctionBasedComplementarityGraphForTwoMethods(new File(dataDir + "MPP-I_predictions.tsv.gz"),
                    new File(dataDir + "PP-I_predictions.tsv.gz"), "MPP-I", "PP", "'red','green','blue'", prThreshold,
                    new File(dataDir + "og2funcs-eggNOG_4.txt.gz"), inGeneOntologyFile,
                    new File(outDir + "/Fig1_fgh_C/PR-" + prThreshold));
        
        if (Fig2_a)
        {
            File[] inPrecisionFiles = {new File(dataDir + "Freshwater_predictions.tsv.gz"), new File(dataDir + "Marine_predictions.tsv.gz"),
                                       new File(dataDir + "Thermal_springs_predictions.tsv.gz"), new File(dataDir + "Soil_predictions.tsv.gz"),
                                       new File(dataDir + "Engineered_predictions.tsv.gz"), new File(dataDir + "Human_predictions.tsv.gz"),
                                       new File(dataDir + "Plants_predictions.tsv.gz")};
            
            String[] methodsNames = {"Freshwater", "Marine", "Thermal springs", "Soil", "Engineered", "Human", "Plants"};
            
            for (double prThreshold : prThresholds)
                drawFunctionBasedComplementarityGraphForMultipleModels(inPrecisionFiles, methodsNames, prThreshold,
                                                                       new File(dataDir + "og2funcs-eggNOG_4.txt.gz"),
                                                                       inGeneOntologyFile,
                                                                       new File(outDir + "/Fig2a/PR-" + prThreshold));
        }
        
        if (Fig3a)
        {
            String outDir3 = outDir + "Fig3a/";
            
            int goFunctionOnWhich_MPP_PerformsBetter = 51540;
            int goFunctionOnWhich_PP_PerformsBetter = 3954;
            
            randomForestFeatureSelection(new File(dataDir + "MPP-I.arff.zip"), new File(dataDir + "PP-I.arff.zip"),
                                         goFunctionOnWhich_MPP_PerformsBetter, goFunctionOnWhich_PP_PerformsBetter,
                                         new File(outDir3));

            computeSimilaritiesBetweenOGs(new File(outDir3 + "MPP-I.csv"), new File(outDir3 + "PP-I.csv"),
                                          new File(outDir3 + "Pearson_correlation_coefficients.txt"));

            drawPearsonCorrelationCoefficientDistributionGraph(new File(outDir3 + "Pearson_correlation_coefficients.txt"),
                                                               new File(outDir3));

            composeGephiFileWithNetwork(new File(outDir3 + "MPP-I.csv"), new File(outDir3 + "PP-I.csv"), 0.7,
                                        new File(outDir3 + "Pearson_correlation_coefficients.txt"),
                                        new File(outDir3 + "CoEvolution_network.gexf"));
        }
        
        if (Fig3c)
        {
            String outDir3 = outDir + "Fig3c/";
            
            int goFunctionOnWhich_MPP_PerformsBetter = 4812;
            int goFunctionOnWhich_PP_PerformsBetter = 6520;
            
            randomForestFeatureSelection(new File(dataDir + "MPP-I.arff.zip"), new File(dataDir + "PP-I.arff.zip"),
                                         goFunctionOnWhich_MPP_PerformsBetter, goFunctionOnWhich_PP_PerformsBetter,
                                         new File(outDir3));

            computeSimilaritiesBetweenOGs(new File(outDir3 + "MPP-I.csv"), new File(outDir3 + "PP-I.csv"),
                                          new File(outDir3 + "Pearson_correlation_coefficients.txt"));

            drawPearsonCorrelationCoefficientDistributionGraph(new File(outDir3 + "Pearson_correlation_coefficients.txt"),
                                                               new File(outDir3));

            composeGephiFileWithNetwork(new File(outDir3 + "MPP-I.csv"), new File(outDir3 + "PP-I.csv"), 0.7,
                                        new File(outDir3 + "Pearson_correlation_coefficients.txt"),
                                        new File(outDir3 + "CoEvolution_network.gexf"));
        }
        
        if (Fig4a)
        {
            File[] inFilesWithAUPRCsAndPredictionsStats = {new File(dataDir + "MPP-I-AUPRCs+prediction_statistics.txt"),
                                                           new File(dataDir + "MPP-16S-5k-AUPRCs+prediction_statistics.txt"),
                                                           new File(dataDir + "MPP-16S-Baseline-AUPRCs+prediction_statistics.txt")};
            String[] classifierNames = {"I", "16S", "Baseline"};
            String[] colors = {"red", "black", "grey"};
            
            drawDistributionOfAUPRCsBoxPlot(inFilesWithAUPRCsAndPredictionsStats, classifierNames, colors, inGeneOntologyFile,
                                            inFileWithFunctionFrequencies, new File(outDir + "/Fig4a/"));
        }
        
        if (Fig4b)
        {
            File[] inFilesWithAUPRCsAndPredictionsStats = {new File(dataDir + "MPP-I-ABUNDANT-AUPRCs+prediction_statistics.txt"),
                                                           new File(dataDir + "MPP-16S-ABUNDANT-AUPRCs+prediction_statistics.txt")};
            String[] classifierNames = {"I", "16S"};
            String[] colors = {"red", "black"};
            
            drawDistributionOfAUPRCsBoxPlot(inFilesWithAUPRCsAndPredictionsStats, classifierNames, colors, inGeneOntologyFile,
                                            inFileWithFunctionFrequencies, new File(outDir + "/Fig4b/"));
        }
        
        if (Fig4c)
        {
            File[] inFilesWithAUPRCsAndPredictionsStats = {new File(dataDir + "MPP-I-RARE-AUPRCs+prediction_statistics.txt"),
                                                           new File(dataDir + "MPP-16S-RARE-AUPRCs+prediction_statistics.txt")};
            String[] classifierNames = {"I", "16S"};
            String[] colors = {"red", "black"};
            
            drawDistributionOfAUPRCsBoxPlot(inFilesWithAUPRCsAndPredictionsStats, classifierNames, colors, inGeneOntologyFile,
                                            inFileWithFunctionFrequencies, new File(outDir + "/Fig4c/"));
        }
    }
}
