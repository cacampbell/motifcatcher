package motifcatcher;

import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: cobalt
 * Date: 30.11.2011
 * Time: 11:44
 * To change this template use File | Settings | File Templates.
 */
public class MakeMotifTree {

    boolean prun;
    Process p;

    public class MC_GetTransfacPWMsRes {
        public boolean EmptyFlag;
        public List<String> DirectoryList;
    }

    public class MakeMotifTreeRes {
        public PhyTree tree;
        public List<String> DirectoryList;
    }

    public class EListEntry {
        public String dir;
        public Double eval;

        public EListEntry(String dir, Double eval) {
            this.dir = dir;
            this.eval = eval;
        }
    }

    public MakeMotifTreeRes run(DataSetProfile dsp) {
        //Description: this function creates a motif tree from the motifs associated
        //with a set of related subsets.
        List<String> DirectoryList;
        PhyTree tree;

        //Retrieve R, and write out a list of directories.
        DirectoryList = MC_GetR(dsp.dirParam.RMEME);

        //output file for significant PWMs.
        String outputfile = dsp.dirParam.MotifTree + "/SignificantR_PWMs";

        //update directory list - cutoff for tree building.
        MC_GetTransfacPWMsRes res = MC_GetTransfacPWMs(DirectoryList,outputfile,dsp.dirParam,dsp.motifTreeParam);

        if (!res.EmptyFlag) {

            //(2) Export TransFac matrices to STAMP, Build tree.
            // ----------------------------------------------------------------------- %

            //explanations of cmd:
            //-tf = TransFac format, input motifs
            //-sd = score distribution file
            //-cc = Pearson Correlation Coefficient (PCC)
            //-align SWU = Smith-Waterman ungapped (SWU)
            //-match = needs to be there for the STAMP program to run properly

            String cmd = dsp.progParam.STAMP + " -tf " + outputfile + " -sd " + dsp.progParam.STAMPScores + " -cc " + dsp.motifTreeParam.ColComparison + "  -align " + dsp.motifTreeParam.AlignmentScheme + " -out " + outputfile + " -match " + dsp.progParam.STAMPDatabase;

            try {
                p = Runtime.getRuntime().exec(new String[] {"/bin/bash", "-c", cmd});
                System.out.println(cmd);
                prun = true;
                Thread x = new Thread() {
                    @Override
                    public void run() {
                        while (prun) {
                           try {
                               while (p.getInputStream().available() > 0) {
                                   System.out.print((char)p.getInputStream().read());
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

            //(3) Cladogram tree of ES, and Motif Families
            // ----------------------------------------------------------------------- %
            //import output tree
            tree = new PhyTree(outputfile + ".tree");
            //PhyTree.ClusterRes r = tree.cluster(10);
            //PhyTree.PhyTreeInfo i = tree.getPhyTreeInfo();
            //System.out.println(i.LeafNames);
        } else {
            tree = new PhyTree();
            System.out.println("No significant related subset motifs to build a tree out of.");
        }
        MakeMotifTreeRes r = new MakeMotifTreeRes();
        r.tree = tree;
        r.DirectoryList = DirectoryList;

        return r;
    }

    //Description:
    //   This function retrieves the Ri-associated motifs (PWMs), and formats
    //   them in Transfac format, so they may be processed by the STAMP
    //   platform.
    private MC_GetTransfacPWMsRes MC_GetTransfacPWMs(List<String> directoryList, String outputfile, DirParam dirParam, MotifTreeParam motifTreeParam) {

        MC_GetTransfacPWMsRes res = new MC_GetTransfacPWMsRes();
        //preset: EmptyFlag is set to true.
        res.EmptyFlag = true;

        //Initialize EList
        //List<EListEntry> EList = new LinkedList<EListEntry>();
        PriorityQueue<EListEntry> eQueue = new PriorityQueue<EListEntry>(100, new Comparator<EListEntry>() {
            public int compare(EListEntry eListEntry, EListEntry eListEntry1) {
                if (eListEntry.eval < eListEntry1.eval) {
                    return -1;
                } else if (eListEntry.eval.equals(eListEntry1.eval)) {
                    return 0;
                } else {
                    return 1;
                }
            }
        });

        //The initial DirectoryList is all directories; this is trimmed down to be
        //only the directories used for motif -construction.
        //int NumAccepted = 0;
        //int EListCandidate = 0;

        for (String dir : directoryList) {

            //retrieve meme output file.
            String file = dirParam.RMEME + "/" + dir + "/meme.txt";

            //extract information from meme output file.
            List<String> extraction = new LinkedList<String>();
            try {
                BufferedReader br = new BufferedReader(new FileReader(file));
                //check for a valid file identifier
                if (br != null) {
                    String line;
                    while ((line = br.readLine()) != null) {
                        extraction.add(line);
                    }
                    br.close();
                }
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }

            //create a dummy Eval, in the case that no PWM exists for this motif.
            Double Eval = null;

            //retrieve PWM coordinates from extraction.
            for (String line : extraction) {
                String[] cols = line.split(" ");
                if ( cols.length >= 9 &&
                        cols[0].equals("letter-probability") &&
                        cols[1].equals("matrix:") &&
                        cols[2].equals("alength=") &&
                        cols[4].equals("w=") &&
                        cols[6].equals("nsites=") &&
                        cols[8].equals("E=")) {

                    //determine the E-value for the motif associated with this R.
                    Eval = Double.valueOf(cols[9]);
                }
            }

            //if the determined Eval is not the dummy Eval, note this value.
            if (Eval != null) {
                eQueue.add(new EListEntry(dir,Eval));
                //Note the name and E-value of this motif.
                //EList.add(EListCandidate++, new EListEntry(dir,Eval));

                //incrememnt the number of accepted (by E-value criterion) R.
                //if (Eval <= motifTreeParam.EValThreshold) {
                //    NumAccepted++;
                //}
            }
        }

        //retrieve relevant E-value Ri.
        //if (EList.size() > 0) {
        if (!eQueue.isEmpty()) {

            //valid R were discovered!
            res.EmptyFlag = false;

            //sort the E-values into ascending order (best matches near the top)
            //Collections.sort(EList, new Comparator<EListEntry>() {
            //    public int compare(EListEntry eListEntry, EListEntry eListEntry1) {
            //        if (eListEntry.eval < eListEntry1.eval) {
            //            return -1;
            //        } else if (eListEntry.eval == eListEntry1.eval) {
            //            return 0;
            //        } else {
            //            return 1;
            //        }
            //    }
            //});

            //revise the Directory List to include only the appropriate R PWMs.
            //int num = Math.min(NumAccepted,motifTreeParam.MaxR);

            res.DirectoryList = new LinkedList<String>();
            for (int i = 0; i < motifTreeParam.MaxR && !eQueue.isEmpty(); i++) {
                res.DirectoryList.add(eQueue.poll().dir);
            }
            //for (int i = 0; i < num; i++) {
            //    res.DirectoryList.add(EList.get(i).dir);
            //}

            //using the corrected DirectoryList, extract the appropriate files.
            for (String dir : res.DirectoryList) {

                //retrieve meme output file.
                String file = dirParam.RMEME + "/" + dir + "/meme.txt";

                //extract information from meme output file.
                List<String> extraction = new LinkedList<String>();
                try {
                    BufferedReader br = new BufferedReader(new FileReader(file));
                    //check for a valid file identifier
                    if (br != null) {
                        String line;
                        while ((line = br.readLine()) != null) {
                            extraction.add(line);
                        }
                        br.close();
                    }
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                } catch (IOException e) {
                    e.printStackTrace();
                }

               //retrieve PWM coordinates from extraction.
                int MatStart = 0;
                float Sites = 0;
                int MatStop = 0;

               for (int j = 0; j < extraction.size(); j++) {

                   String[] cols = extraction.get(j).split(" ");

                    if (cols.length > 8 &&
                            cols[0].equals("letter-probability") &&
                            cols[1].equals("matrix:") &&
                            cols[2].equals("alength=") &&
                            cols[4].equals("w=") &&
                            cols[6].equals("nsites=") &&
                            cols[8].equals("E=")) {

                        MatStart = j + 2;
                        Sites = Float.valueOf(cols[7]);
                    }

                    if ( cols.length > 3 &&
                            cols[0].equals("\tMotif") &&
                            cols[1].equals("1") &&
                            cols[2].equals("regular") &&
                            cols[3].equals("expression") ) {
                        MatStop = j - 3;
                   }
               }

               //write the PWM to the output file.
               try {
                   BufferedWriter bs = new BufferedWriter(new FileWriter(outputfile, true));

                   //header
                   bs.write("DE\t"+dir+"\n");

                   //matrix
                   int PosCounter = -1;
                   for (int j = MatStart; j < MatStop; j++) {
                       //record the position counter - start counting at 0
                       PosCounter++;

                       String[] cols = extraction.get(j).split(" ");
                       //columns are frequencies, assembled from the number of sites.
                       int col1 = Math.round(Float.valueOf(cols[1])*Sites);
                       int col2 = Math.round(Float.valueOf(cols[3])*Sites);
                       int col3 = Math.round(Float.valueOf(cols[5])*Sites);
                       int col4 = Math.round(Float.valueOf(cols[7])*Sites);

                       bs.write(PosCounter + "\t");
                       bs.write(col1 + "\t");
                       bs.write(col2 + "\t");
                       bs.write(col3 + "\t");
                       bs.write(col4 + "\n");
                   }

                   //footer
                   bs.write("XX\n");

                   //close file.
                   bs.close();
               } catch (FileNotFoundException e) {
                   e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
               } catch (IOException e) {
                   e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
               }
            }
        }
        return res;
    }

    public List<String> MC_GetR(String RDirectory) {
        //Description: This function determines the appropriate related subsets
        //for future analysis
        //
        //Inputs:
        //   RDirectory:
        //       Directory containing all related subsets
        //
        //Outputs:
        //   LastIterations:
        //       MEME output dirctories of last iterations (related subsets)
        //   FileList:
        //       list of files/directories in EndState folder, directly after
        //       End State generation step.

        //create a list of files in the end state directory.
        File dir = new File(RDirectory);
        File[] list = dir.listFiles();

        //determine total number of clusters to extract
        int numcluster = 0;
        for (int i = 1; i < list.length; i++) {
            String filename = list[i].getName();
            if (filename.contains("Seed")) {

                String [] A = filename.split("v");
                String [] B = A[0].split("Seed");

                int value = Integer.valueOf(B[1]);

                if (value >= numcluster) {
                   numcluster = value;
                }
            }
        }

        //initialize LastIterations output
        List<String> LastIterations = new LinkedList<String>();

        //for each cluster, find the top iteration.
        for (int i = 1; i <= numcluster; i++) {

            String lookfor = "Seed" + i + "v";
            int LastIteration = 0;

            //find top iteration
            for (File f : list) {
               if (f.getName().contains(lookfor)) {
                   String[] A = f.getName().split("v");

                   int value = Integer.valueOf(A[1]);

                   if (value >= LastIteration) {
                      LastIteration = value;
                   }

               }
            }

            //write final iteration information
            LastIterations.add("Seed"+ i +"v"+ LastIteration);
        }
        return LastIterations;
    }
}
