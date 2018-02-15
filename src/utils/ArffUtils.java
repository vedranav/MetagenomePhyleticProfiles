package utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import static utils.FileUtils.findReaderType;

/**
 *
 * @author Vedrana Vidulin [vedrana.vidulin@ijs.si]
 */
public class ArffUtils
{
    public static void selectSubsetOfLabelsFromHMCArff(List<String> labels, File inArffFile, File outArffFile) throws IOException
    {
        FileWriter fw = new FileWriter(outArffFile);
        BufferedWriter bw = new BufferedWriter(fw);
        
        BufferedReader br = findReaderType(inArffFile);
        
        String line;
        
        while((line=br.readLine()) != null)
            if (line.startsWith("@ATTRIBUTE class hierarchical"))
            {
                line = line.substring(line.indexOf("hierarchical") + 12).trim();
                
                Set<String> labelsInArff = new HashSet<>();
                
                String[] parts = line.split(",");
                
                for (String part : parts)
                {
                    String firstLabel = part.substring(0, part.indexOf("/"));
                    String secondLabel = part.substring(part.indexOf("/") + 1);
                    
                    if (!firstLabel.equals("root"))
                        labelsInArff.add(firstLabel);
                    
                    if (!secondLabel.equals("root"))
                        labelsInArff.add(secondLabel);
                }
                
                labels.retainAll(labelsInArff);
                
                for (String label : labels)
                    bw.write("@ATTRIBUTE " + label + "\t{0,1}\n");
            }
            else if (line.startsWith("@") || line.isEmpty())
                bw.write(line + "\n");
            else
            {                
                bw.write(line.substring(0, line.lastIndexOf(",")));
                
                String labelsStr = line.substring(line.lastIndexOf(",") + 1).trim();
                
                String[] parts = labelsStr.split("@");
                List<String> instanceLabels = Arrays.asList(parts);
                
                for (String label : labels)
                    if (instanceLabels.contains(label))
                        bw.write(", 1");
                    else
                        bw.write(", 0");
                
                bw.write("\n");
            }
        
        bw.close();
    }
}
