package motifcatcher;

import java.io.File;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cobalt
 * Date: 01.12.2011
 * Time: 12:27
 * To change this template use File | Settings | File Templates.
 */
public class MotifCatcher {
    public static void main (String[] args) {

        long start = System.currentTimeMillis();
        String seqFile = "seqfile.fasta";
        double clusterThreshold = 20;
        double FPFrequencyThreshold = 0.1;

        DirParam directories = new DirParam();
        boolean directoryRoot = new File(directories.ROOT).mkdirs();
        boolean directoryRSeqs = new File(directories.RSeqs).mkdirs();
        boolean directoryRMeme = new File(directories.RMEME).mkdirs();
        boolean directoryMotifTree = new File(directories.MotifTree).mkdirs();
        boolean directoryFP = new File(directories.FP).mkdirs();

        RParam R_parameters = new RParam();

        MEMEParam MEME_parameters = new MEMEParam();

        ProgParam programLocations = new ProgParam();

        DataSetProfile dsp = new DataSetProfile();
        dsp.motifTreeParam = new MotifTreeParam();
        dsp.progParam = programLocations;
        dsp.dirParam = directories;
        dsp.memeParam = MEME_parameters;

        MakeR mr = new MakeR();
        MakeR.MakeRResult mr_res = mr.run(
              seqFile,
              directories,
              R_parameters,
              MEME_parameters,
              programLocations
        );

        MakeMotifTree mmt = new MakeMotifTree();
        MakeMotifTree.MakeMotifTreeRes mmt_res = mmt.run(dsp);

        TreeWithClusters twc = new TreeWithClusters(dsp, seqFile, mmt_res);
        twc.ComputeClusters(clusterThreshold);
        twc.ComputeFPs(directories.FP,FPFrequencyThreshold);

        MakeMotifMap mmm = new MakeMotifMap();
        String s = mmm.MC_MakeMotifMap(directories, seqFile, MEME_parameters.RevComp);

        int cols = l.get(0).size();
        for (List line : l) {
            for (int i = 0; i < cols; i++) {
                if (line.size()>i && line.get(i) != null)
                    System.out.printf("%20s",line.get(i));
                else
                    System.out.printf("%20s","");
            }
            System.out.println();
        }

        long end = System.currentTimeMillis();
        System.out.println(end-start);

        //TODO
        //CompareOccurenceAndLocalization col = new CompareOccurenceAndLocalization();
        //col.MC_CompareOccurrenceAndLocalization();
    }
}
