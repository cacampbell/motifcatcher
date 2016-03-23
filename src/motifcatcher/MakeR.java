package motifcatcher;

import com.sun.java.swing.plaf.motif.resources.motif;
import com.sun.tools.jdi.GenericListeningConnector;
import motifcatcher.utils.FastaInfo;
import motifcatcher.utils.FastaTools;
import motifcatcher.utils.Sequence;
import phylosort.Path;
import sun.awt.image.ImageWatched;
import sun.tools.asm.CatchData;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: cobalt
 * Date: 28.11.2011
 * Time: 15:23
 * To change this template use File | Settings | File Templates.
 */
public class MakeR {

    Process p;
    boolean prun = true;

    public class MakeRResult {
        public List<List<String>> subsetList;
        public MEMEParam MEME_parameters;
    }

    /**
     * This function generates Ri according to one of the three R-determination
     * protocols.
    %--------------------------------------------------------------------------
    %Inputs:
    %
    %   SeqFile:
    %       input data set of sequences
    %   Directories:
    %       directories for output
    %   R_parameters:
    %       details for creating an 'R'
    %   MEME_parameters:
    %       details for motif search
    %   ProgramLocations:
    %       location of MEME and MAST on system

    %Outputs:
    %   SubsetList:
    %       A single Ri.  All other outputs created implicitly
    */

    public MakeRResult run(String seqFile, DirParam directories, RParam R_parameters, MEMEParam MEME_parameters, ProgParam programLocations) {

        MakeRResult res = new MakeRResult();
        res.MEME_parameters = MEME_parameters;

        // import sequences
        //List<Sequence> sequence = fastaread(seqFile);
        FastaInfo fastaInfo = FastaTools.fastaGetInfo(seqFile);

        List<List<String>> subsetList = null;

        // if a background file is not provided, create one from the input data set.
        if (MEME_parameters.bfile == null) {

            // determine total number of characters in the whole dataset.  This value is
            // used to decide what order of background to use.
            int totalCharacters = fastaInfo.numCharacters;


            // determine order of background file based on the total size of all input
            // sequences.
            String order = null;
            if (totalCharacters >= 2560) {
                order = " 3 ";
            } else if (totalCharacters > 1000) {
                order = " 2 ";
            } else if (totalCharacters > 500) {
                order = " 1 ";
            } else {
                order = " 0 ";
            }

            String type = " ";
            if (MEME_parameters.Alphabet.equals("protein")) {
                type = " -p ";
            }

            //Create background file for meme analysis
            String bfile = directories.RMEME + "/bfile.model";
            String cmd = programLocations.FastaGetMarkov + type + " -m " + order + " < " + seqFile +  " > " + bfile;
            System.out.println(cmd);
            res.MEME_parameters.bfile = bfile;

            try {
                p = Runtime.getRuntime().exec(new String[] {"/bin/bash", "-c", cmd});
                Thread x = new Thread() {
                    @Override
                    public void run() {
                        while (prun) {
                            try {
                                while (p.getErrorStream().available() > 0) {
                                    System.out.print((char)p.getErrorStream().read());
                                }
                            } catch (IOException e) {
                                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                            }
                        }
                    }
                };
                x.start();
                p.waitFor();
                prun = false;

            } catch (IOException e) {
                e.printStackTrace();
            } catch (InterruptedException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

        }

        // TODO: parallelize here !!

        Thread[] mt = new Thread[(int)R_parameters.Seeds];
        for (int i = 1; i <= R_parameters.Seeds; i++) {

            Runnable r = new MemeThread(fastaInfo,R_parameters,directories,programLocations,MEME_parameters,i,seqFile);
            mt[i-1] = new Thread(r);
            mt[i-1].start();
        }

        for (int i = 1; i <= R_parameters.Seeds; i++) {
            try {
                mt[i-1].join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }


        /*
        for (int i = 1; i <= R_parameters.Seeds; i++) {

            Runnable r = new MemeThread(sequence,R_parameters,directories,programLocations,MEME_parameters,i,seqFile);
            mt[i-1] = new Thread(r);
            mt[i-1].start();

            // -------- Generate a random seed --------------------------------------- %

            //generate a random permutation
            //from these randomly organized numbers, select as many as needed to build a
            //subset.
            Vector<Integer> subsetNums = randperm(sequence.size(),R_parameters.SeedSize);

            //build an initial list of genes.
            //each cell contains the genes used for that particular run.
            //the subset List is rebuilt for every run.
            subsetList = new LinkedList<List<String>>();

            subsetList.add(0, new LinkedList<String>());
            for (Integer j : subsetNums) {
                 subsetList.get(0).add(sequence.get(j).header);
            }

            // -------- initializations ---------------------------------------------- %
            int iterations = 0;
            boolean clusterCompleted = false;

            // -------- Iteratively re-cluster --------------------------------------- %
            while(!clusterCompleted) {

                // -------- Extract sequences, write to files ---------------------------- %

                //check to see that the .fasta files contains enough sequences for
                //MEME to run.  If so, build the file and move on.
                if (subsetList.get(iterations).size() >= 2) {

                    //make output file.
                    // TODO check this !!!!
                    String filetitle = directories.RSeqs + "/Seed" + i + "v" + iterations + ".fasta";
                    // TODO: Do that header lookup more efficient either by just using the seqNum or by using a map !!!!!
                    // TODO: !!!! Extremely inefficient !!!!!!
                    for (int j = 0; j < subsetList.get(iterations).size(); j++) {
                        //find sequence matching this header.
                        for ( int q = 0; q < sequence.size(); q++) {
                            if (sequence.get(q).header.equals(subsetList.get(iterations).get(j))) {
                                //optional print statements
                                //disp(['j=' num2str(j) ': ' SubsetList{1,iterations+1}{1,j}])
                                //disp(['q=' num2str(q) ': ' header{1,q}]) {
                                // TODO correct according to original code !!!
                                fastawrite(filetitle,subsetList.get(iterations).get(j),sequence.get(q).content);
                            }
                        }
                    }

                    //MEME analysis
                    //direct the meme program to run according to input parameters.
                    String memeout = directories.RMEME + "/Seed" + i + "v" + iterations;
                    String cmd = programLocations.MEME + " " + filetitle + " -minw " + MEME_parameters.MinWidth + " -maxw " + MEME_parameters.MaxWidth + " -" + MEME_parameters.Alphabet + " -bfile " + MEME_parameters.bfile + " -mod zoops -oc " + memeout + " -nostatus -maxsize 1000000";

                    //option: search the reverse complement, or not.
                    if (MEME_parameters.RevComp) {
                       cmd += " -revcomp";
                    }

                    try {
                        p = Runtime.getRuntime().exec(new String[] {"/bin/bash", "-c", cmd});
                        System.out.println(cmd);
                        prun = true;
                        Thread x = new Thread() {
                            @Override
                            public void run() {
                                while (prun) {
                                    try {
                                        while (p.getErrorStream().available() > 0) {
                                            System.out.print((char)p.getErrorStream().read());
                                        }
                                    } catch (IOException e) {
                                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                                    }
                                }
                            }
                        };
                        x.start();
                        p.waitFor();
                        prun = false;

                    } catch (IOException e) {
                        e.printStackTrace();
                    } catch (InterruptedException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }

                    // a simple ZOOPs search completes a run for SearchStyle 1.
                    // ZOOPS -> MAST -> ZOOPS again is complete for SearchStyle 2.
                    // the iterative process (until convergence) is the default,
                    // which continues
                    if (R_parameters.SearchStyle == 1) {
                        break;
                    } else if (R_parameters.SearchStyle == 2 && iterations == 1) {
                        break;
                    }

                    //MAST analysis
                    String mastinputfile = memeout + "/meme.txt";

                    cmd = programLocations.MAST + ' ' + mastinputfile + " " + seqFile + " -oc " + memeout + " -nostatus -nohtml -bfile " + MEME_parameters.bfile + " -ev " + MEME_parameters.MastThreshold;

                    //option: search for reverse compliment, or not.
                    if (MEME_parameters.RevComp) {
                        cmd += " -norc ";
                    }

                    try {
                       p = Runtime.getRuntime().exec(new String[] {"/bin/bash", "-c", cmd});
                       System.out.println(cmd);
                       prun = true;
                       Thread x = new Thread() {
                           @Override
                           public void run() {
                               while (prun) {
                                   try {
                                       while (p.getErrorStream().available() > 0) {
                                           System.out.print((char)p.getErrorStream().read());
                                       }
                                   } catch (IOException e) {
                                       e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                                   }
                               }
                           }
                       };
                       x.start();
                       p.waitFor();
                       prun = false;
                    } catch (IOException e) {
                        e.printStackTrace();
                    } catch (InterruptedException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }


                    //depending on the variety of search, the user may call for
                    //different outputs.
                    List<String> sl = MC_GetSubsetList(memeout, sequence);

                    //check for missing MAST file
                    if (sl != null) {
                        subsetList.add(iterations+1,sl);
                    }

                    //check for empty MAST-searches; in this case just revert to
                    //previous.
                    if (subsetList.size() <= iterations+1) {
                        subsetList.add(iterations+1, subsetList.get(iterations));
                        clusterCompleted = true;
                        break;
                    } else if (subsetList.get(iterations+1).size() == 1) {
                        subsetList.set(iterations+1, subsetList.get(iterations));
                        clusterCompleted = true;
                        break;
                    }

                    //see if newest gene list is a repeat
                    for (int j = 0; j < subsetList.size(); j++) {
                        for (int k = 0; k < j-1; k++) {
                            // TODO sort items!!!!
                            if (subsetList.get(j).equals(subsetList.get(k))) {
                                clusterCompleted = true;
                                break;
                            }
                        }
                    }

                    //increment the iteration number
                    iterations = iterations + 1;

                } else { //exit the loop
                   clusterCompleted = true;
                }
            }
            //optional print statement
            System.out.println("Related subset " + i + " determined.");
        }
        */
        // TODO ????? subsetlist is overriden in for ???!???!
        res.subsetList = subsetList;
        return res;
    }

    public class MemeThread implements Runnable {

        RParam R_parameters;
        private FastaInfo fastaInfo;
        private DirParam directories;
        private ProgParam programLocations;
        private MEMEParam MEME_parameters;
        private int i;
        private String seqFile;

        private LinkedList<List<String>> subsetList;

        public MemeThread(FastaInfo fastaInfo, RParam r_parameters, DirParam directories, ProgParam programLocations, MEMEParam MEME_parameters, int i, String seqFile) {
            this.fastaInfo = fastaInfo;
            R_parameters = r_parameters;
            this.directories = directories;
            this.programLocations = programLocations;
            this.MEME_parameters = MEME_parameters;
            this.i = i;
            this.seqFile = seqFile;
        }

        public void run() {
            // -------- Generate a random seed --------------------------------------- %

            //generate a random permutation
            //from these randomly organized numbers, select as many as needed to build a
            //subset.
            HashSet<Integer> subsetNums = randperm(fastaInfo.numSequences,R_parameters.SeedSize);

            //build an initial list of genes.
            //each cell contains the genes used for that particular run.
            //the subset List is rebuilt for every run.
            List<List<String>> subsetList = new LinkedList<List<String>>();
            //for (Integer j : subsetNums) {
            //     subsetList.get(0).add(sequence.get(j).header);
            //}
            Map<String, motifcatcher.utils.Sequence> subsetMap = null;
            try {
                subsetMap = FastaTools.getSequences(new FileInputStream(seqFile), subsetNums);
                subsetList.add(0, new LinkedList<String>(subsetMap.keySet()));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }

            // -------- initializations ---------------------------------------------- %
            int iterations = 0;
            boolean clusterCompleted = false;

            // -------- Iteratively re-cluster --------------------------------------- %
            while(!clusterCompleted) {

                // -------- Extract sequences, write to files ---------------------------- %

                //check to see that the .fasta files contains enough sequences for
                //MEME to run.  If so, build the file and move on.
                if (subsetList.get(iterations).size() >= 2) {

                    //make output file.
                    // TODO check this !!!!
                    List<String> lookup = new LinkedList<String>();
                    String filetitle = directories.RSeqs + "/Seed" + i + "v" + iterations + ".fasta";
                    for (String sh : subsetList.get(iterations)) {
                        //find sequence matching this header.
                        //System.out.println("trying to find:" + sh);
                        Sequence s = subsetMap.get(sh);
                        if (s != null) {
                            FastaTools.fastaAppend(filetitle,s.header,s.content);
                        } else {
                            //System.out.println("Sequence not found in map !!!!");
                            lookup.add(sh);
                        }
                    }
                    List<Sequence> list = null;
                    try {
                        list = FastaTools.getSequences(new FileInputStream(seqFile), lookup);
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }
                    for (Sequence s : list) {
                        subsetMap.put(s.header,s);
                        FastaTools.fastaAppend(filetitle,s.header,s.content);
                    }

                    //MEME analysis
                    //direct the meme program to run according to input parameters.
                    String memeout = directories.RMEME + "/Seed" + i + "v" + iterations;
                    String cmd = programLocations.MEME + " " + filetitle + " -minw " + MEME_parameters.MinWidth + " -maxw " + MEME_parameters.MaxWidth + " -" + MEME_parameters.Alphabet + " -bfile " + MEME_parameters.bfile + " -mod zoops -oc " + memeout + " -nostatus -maxsize 1000000";

                    //option: search the reverse complement, or not.
                    if (MEME_parameters.RevComp) {
                       cmd += " -revcomp";
                    }

                    try {
                        p = Runtime.getRuntime().exec(new String[] {"/bin/bash", "-c", cmd});
                        System.out.println(cmd);
                        prun = true;
                        Thread x = new Thread() {
                            @Override
                            public void run() {
                                while (prun) {
                                    try {
                                        while (p.getErrorStream().available() > 0) {
                                            System.out.print((char)p.getErrorStream().read());
                                        }
                                    } catch (IOException e) {
                                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                                    }
                                }
                            }
                        };
                        x.start();
                        p.waitFor();
                        prun = false;

                    } catch (IOException e) {
                        e.printStackTrace();
                    } catch (InterruptedException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }

                    // a simple ZOOPs search completes a run for SearchStyle 1.
                    // ZOOPS -> MAST -> ZOOPS again is complete for SearchStyle 2.
                    // the iterative process (until convergence) is the default,
                    // which continues
                    if (R_parameters.SearchStyle == 1) {
                        break;
                    } else if (R_parameters.SearchStyle == 2 && iterations == 1) {
                        break;
                    }

                    //MAST analysis
                    String mastinputfile = memeout + "/meme.txt";

                    cmd = programLocations.MAST + ' ' + mastinputfile + " " + seqFile + " -oc " + memeout + " -nostatus -nohtml -bfile " + MEME_parameters.bfile + " -ev " + MEME_parameters.MastThreshold;

                    //option: search for reverse compliment, or not.
                    if (MEME_parameters.RevComp) {
                        cmd += " -norc ";
                    }

                    try {
                       p = Runtime.getRuntime().exec(new String[] {"/bin/bash", "-c", cmd});
                       System.out.println(cmd);
                       prun = true;
                       Thread x = new Thread() {
                           @Override
                           public void run() {
                               while (prun) {
                                   try {
                                       while (p.getErrorStream().available() > 0) {
                                           System.out.print((char)p.getErrorStream().read());
                                       }
                                   } catch (IOException e) {
                                       e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                                   }
                               }
                           }
                       };
                       x.start();
                       p.waitFor();
                       prun = false;
                    } catch (IOException e) {
                        e.printStackTrace();
                    } catch (InterruptedException e) {
                        e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
                    }


                    //depending on the variety of search, the user may call for
                    //different outputs.
                    List<String> sl = null;
                    try {
                        sl = MC_GetSubsetList(memeout, new FileInputStream(seqFile));
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }

                    //check for missing MAST file
                    if (sl != null) {
                        if (sl.size() >= (int)(0.9 * fastaInfo.numSequences))
                        {
                            clusterCompleted = true;
                            System.out.println("Too many sequences in MAST file!");
                        } else {
                            subsetList.add(iterations+1,sl);
                        }
                    }

                    //check for empty MAST-searches; in this case just revert to
                    //previous.
                    if (subsetList.size() <= iterations+1) {
                        subsetList.add(iterations+1, subsetList.get(iterations));
                        clusterCompleted = true;
                        System.out.println("MAST file is empty.");
                        break;
                    } else if (subsetList.get(iterations+1).size() == 1) {
                        subsetList.set(iterations+1, subsetList.get(iterations));
                        clusterCompleted = true;
                        System.out.println("Only one sequence left.");
                        break;
                    }

                    //see if newest gene list is a repeat
                    for (int j = 0; j < subsetList.size(); j++) {
                        for (int k = 0; k < j-1; k++) {
                            // TODO sort items!!!!
                            if (subsetList.get(j).equals(subsetList.get(k))) {
                                clusterCompleted = true;
                                System.out.println("Subset is a repeat.");
                                break;
                            }
                        }
                    }

                    //increment the iteration number
                    //TODO: make a parameter!!!!
                    if (iterations < 10) {
                        iterations = iterations + 1;
                    } else {
                        clusterCompleted = true;
                        System.out.println("Maximum number of iterations reached.");
                    }


                } else { //exit the loop
                    clusterCompleted = true;
                    System.out.println("Less than two sequences left.");
                }
            }
            //optional print statement
            System.out.println("Related subset " + i + " determined.");
        }

        public List<String> MC_GetSubsetList(String mastdirectory, InputStream sequence) {
            // Extract list of sequence entries from a mast .txt output file.  The format of the
            // file depends on whether or not reverse-complement strands are allowed.

            HashSet<String> geneList = new HashSet<String>();

            String file = mastdirectory + "/mast.txt";

            try {
                BufferedReader br = new BufferedReader(new FileReader(file));
                List<String> extraction = new LinkedList<String>();

                //check for a valid file identifier
                if (br != null) {

                    String line;
                    boolean startIndex = false;
                    boolean stopIndex = false;
                    int startCounter = 2;
                    while ((line = br.readLine()) != null) {

                        if (!startIndex && line.startsWith("SEQUENCE") &&
                                line.contains("NAME") &&
                                line.contains("DESCRIPTION") &&
                                line.contains("E-VALUE") &&
                                line.contains("LENGTH"))
                        {
                            startIndex = true;
                            startCounter = 2;
                        }

                        if (startIndex && line.equals("********************************************************************************")) {
                            stopIndex = true;
                        }

                        if (startIndex && !stopIndex) {
                            if (startCounter > 0) {
                                startCounter--;
                            } else {
                                extraction.add(line);
                            }
                        }
                    }

                    br.close();
                }


                // -1 to skip last empty line
                for (int j = 0; j < extraction.size() - 1; j++) {
                    //System.out.println("MAST Adding to Subset:" + extraction.get(i).split(" ")[0]);
                    geneList.add(extraction.get(j).split(" ")[0]);
                }

                System.out.println("Headers in " + i + " before: "+ geneList.size());
                List<String> res = FastaTools.completeHeaders(sequence, geneList);
                System.out.println("Headers in " + i + " after: "+ res.size());

                //sort the output list in faux-alphabetical order
                Collections.sort(res);

                return res;

            } catch (FileNotFoundException e) {
                //e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }
            return null;
        }
    }

    private HashSet<Integer> randperm(int maxSize, int seedSize) {
        //TODO optimize !!!
        HashSet<Integer> rands = new HashSet<Integer>(seedSize);
        Random r = new Random();
        while (rands.size() < seedSize) {
            rands.add(r.nextInt(maxSize));
        }
        return rands;
    }


    public static void main (String[] args) {
        String seqFile = "/Users/cobalt/Downloads/MotifCatcher_v1.01/MotifCatcher_v1.01/Synthetic_Data.fasta";
        DirParam directories = new DirParam();

        new File(directories.ROOT).mkdirs();
        new File(directories.RSeqs).mkdirs();
        new File(directories.RMEME).mkdirs();
        //new File(directories.MotifTree).mkdirs();
        //new File(directories.FP).mkdirs();

        RParam R_parameters = new RParam();
        MEMEParam MEME_parameters = new MEMEParam();
        ProgParam programLocations = new ProgParam();
        new MakeR().run(seqFile,directories,R_parameters,MEME_parameters,programLocations);
    }
}