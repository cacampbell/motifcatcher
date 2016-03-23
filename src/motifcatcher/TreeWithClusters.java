package motifcatcher;

import motifcatcher.utils.FastaTools;
import motifcatcher.utils.Sequence;
import java.io.*;
import java.util.*;

public class TreeWithClusters {

    //define all global variables
    MakeMotifTree.MakeMotifTreeRes mMTRes;
    MakeClustersRes mkr;

    boolean prun;
    Process p;
    private DataSetProfile dataSetProfile;
    private String seqFile;


    public TreeWithClusters(DataSetProfile dataSetProfile, String seqFile, MakeMotifTree.MakeMotifTreeRes makeMotifTreeRes) {
        this.dataSetProfile = dataSetProfile;
        this.seqFile = seqFile;
        this.mMTRes = makeMotifTreeRes;
    }

    public class MakeClustersRes {
        List<Family> Families;
        AllMotifsAndLocations AllMotifsAndLocations;
    }

    public class SeedEntry {
        private Integer seedNum;
        private String seedPath;

        public SeedEntry(Integer seedNum, String seedPath) {
            this.seedNum = seedNum;
            this.seedPath = seedPath;
        }
    }

    public class Family {
        private List<Integer> a;
        private List<String> inputSet;
        private List<SeedEntry> seedList;

        public List<SeedEntry> getSeedList() {
            return seedList;
        }

        public void setSeedList(List<SeedEntry> seedList) {
            this.seedList = seedList;
        }

        public Family() {
            a = new LinkedList<Integer>();
        }

        public void add(int k) {
            a.add(k);
        }

        public int size() {
            return a.size();
        }

        public List<String> getInputSet() {
            return inputSet;
        }

        public void setInputSet(List<String> inputSet) {
            this.inputSet = inputSet;
        }

        public Integer get(int j) {
            return a.get(j);
        }

    }

    public class AllMotifsAndLocations {
        LinkedList<Family> fs;

        public void add(Family f) {
            fs.add(f);
        }
    }

    public void ComputeClusters(double clusteringThreshold) {

        System.out.println("Processing request...");

        //evaluate clusters
        PhyTree.ClusterRes cr = mMTRes.tree.cluster(clusteringThreshold);
        System.out.println("found " + Collections.max(cr.LeafClusters) + " Clusters.");

        //TODO visualize
        //cmap = colormap(lines);
        //for k = 1:max(LeafClusters)
        //    set(h.BranchLines(NodeClusters == k),"Color",cmap(k,:))
        //}

        System.out.println("Computing variables necessary for R-associated motif display...");

        //Create clusters
        mkr = MC_MakeClusters(dataSetProfile, cr, seqFile, clusteringThreshold);
    }

    private MakeClustersRes MC_MakeClusters(DataSetProfile dataSetProfile, PhyTree.ClusterRes clusterRes, String seqFile, double clusteringThreshold) {

        // (1) Compute clusters and retrieve Motif Tree
        // ----------------------------------------------------------------------- //

        //retrieve information from tree
        PhyTree.PhyTreeInfo phyTreeInfo = mMTRes.tree.getPhyTreeInfo();

        // (2) Create "Families" output.
        // ----------------------------------------------------------------------- //
        
        //initialize "Families" structure
        List<Family> Families = new LinkedList<Family>();

        //transfer information from Motif Tree output to AllCoreMotifs form,
        //and note family grouping.
        for (int j = 1; j <= Collections.max(clusterRes.LeafClusters); j++) {
            Family f = new Family();
            for (int i = 0; i < phyTreeInfo.LeafNames.size(); i++) {
               if (clusterRes.LeafClusters.get(i) == j) {
                   //find "AllCoreMotifs" number, from the "Leaf Nodes" name
                   //argument.
                    for (int k = 0; k < mMTRes.DirectoryList.size(); k++) {
                        if (mMTRes.DirectoryList.get(k).contains(phyTreeInfo.LeafNames.get(i))) {
                            f.add(k);
                        }
                    }
               }
           }
           Families.add(f);
        }
            
        //Sort families by length (bubble sort).
        Collections.sort(Families,new Comparator<Family>() {
            public int compare(Family f1, Family f2) {
                return (f1.size() - f2.size());
            }
        });
            
        // (3) Create "AllMotifsAndLocations" output.
        // ----------------------------------------------------------------------- //

        List AllMotifsAndLocations = new ArrayList(Families.size());
        TreeWithClusters.AllMotifsAndLocations aml = new AllMotifsAndLocations();
        // Each cell entry: a Final Cluster Family.  Within each cell entry, each
        // row is a particular ES, and each column a particular site.  In other
        // words - each cell is a mini "Membership Matrix" for each family.
        for (Family f : Families) {
            List<String> InputSet = new LinkedList<String>();
            List<SeedEntry> seed = new LinkedList<SeedEntry>();

            for (int j = 0; j < f.size(); j++) {
                InputSet.add(dataSetProfile.dirParam.RMEME + "/" + mMTRes.DirectoryList.get(f.get(j)));

                //build header for future output structure.
                seed.add(new SeedEntry(f.get(j), mMTRes.DirectoryList.get(f.get(j))));
            }
            f.setSeedList(seed);

            //retrieve each familial membership matrix.
            MC_MakeMembershipMatrix(f, InputSet, seqFile, dataSetProfile.memeParam.RevComp);

            aml.add(f);
        }

        MakeClustersRes res = new MakeClustersRes();
        res.Families = Families;
        res.AllMotifsAndLocations = aml;
        return res;
    }

    private List MC_MakeMembershipMatrix(Family seed, List<String> inputSet, String seqFile, boolean revComp) {

        //This function produces a membership matrix for cases where top motifs are
        //determined by direct motif comparison (and clustering).
        //
        //Inputs:
        //   InputsSet:
        //       cell array of MEME output directories (whole path)
        //   SeqFileHeader:
        //       .fasta file of all sequences in the superset.
        //
        //Ouputs:
        //   MM:
        //       merged membership matrix
        //   MotifLengths:
        //       length of motifs, ordered as same as InputSet.
        //
        //Algorithm:
        //   (1) Extract all information from all directories in input set, and
        //       create a miniature membership matrix for each set. 
        //   (2) Merge all sites from each directory into a single membership
        //       matrix.
        //   (3) sort membership matrix by multiplicity of sites included.
        
        //(1) Extract all information each file in the InputSet
        // ----------------------------------------------------------------------- //
        
        //initalize extraction cell array, and MotifLengths.
        List<List<List>> mats = new ArrayList<List<List>>(inputSet.size());
        List<Integer> MotifLengths = new ArrayList<Integer>(inputSet.size());

        for (int i = 0; i < inputSet.size(); i++) {

            //retrieve meme output file.
            String file = inputSet.get(i) + "/meme.txt";
            mats.add(i,new LinkedList<List>());

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
            Integer startlook = null;
            Integer stoplook = null;
            Integer MotifLength = 0;


            //retrieve PWM coordinates from extraction.
            for (int j = 0; j < extraction.size(); j++) {
                String line = extraction.get(j);
                String[] cols = line.split("\\s+");

                //revcomp version
                if (revComp) {
                    if ( cols.length >= 5 &&
                            cols[0].equals("Sequence") &&
                            cols[1].equals("name") &&
                            cols[2].equals("Strand") &&
                            cols[3].equals("Start") &&
                            cols[4].equals("P-value") &&
                            cols[5].equals("Site")) {

                        startlook = j + 1;
                    }
                } else {
                    if ( cols.length >= 4 &&
                            cols[0].equals("Sequence") &&
                            cols[1].equals("name") &&
                            cols[2].equals("Start") &&
                            cols[3].equals("P-value") &&
                            cols[4].equals("Site")) {

                        startlook = j + 2;
                    }
                }

                if ( cols.length > 4 &&
                    cols[1].equals("Motif") &&
                    cols[2].equals("1") &&
                    cols[3].equals("block") &&
                    cols[4].equals("diagrams")) {

                        stoplook = j - 3;
                }

                if ( cols.length >= 4 &&
                    cols[0].equals("MOTIF") &&
                    cols[1].equals("1") &&
                    cols[2].equals("width") &&
                    cols[3].equals("=")) {

                    MotifLength = Integer.valueOf(cols[4]);
                    MotifLengths.add(i,MotifLength);
                }

            }

            //Find all sites.  Look for key character ":" in the appropriate
            //range, or key character "[" in the appropriate range.
            List<List<String>> Sites = new LinkedList<List<String>>();

            for (int j = startlook; j < stoplook; j++) {
                String line = extraction.get(j);
                String[] cols = line.split("\\s+");

                if (  cols[0].contains("[") ||
                      cols[0].contains(":") ||
                      cols[0].contains("|") ||
                      cols[0].contains("iY") ) {
                    //the "iY" search key was added for a yeast analysis, "i" for intergenic, "Y" for yeast

                    List<String> entry = new LinkedList<String>();
                    entry.add(cols[0]); //Site name
                    if (revComp) {
                        entry.add(cols[2]); //Start
                    } else {
                        entry.add(cols[1]); //Start
                    }
                    Sites.add(entry);
                }
            }
                   
            /*
            Find the full name of each Site by searching for it in the
            SeqFile.

            Note that here, (+) and (-)-strandedness are not considered;
            this is because the interesting form of the motif might
            actually be in reverse-compliment orientation.

            Instead, each motif length is evaluated by forward and reverse
            compliment, at later stages of analysis.
            TODO counter = Sites.length
            */

            for (int k = 0; k < Sites.size(); k++) {
                mats.get(i).add(k, new LinkedList<String>());

                String seq = FastaTools.findNameFor(seqFile, Sites.get(k).get(0));

                //record site
                mats.get(i).get(k).add(seq);

                //record span of motif
                LinkedList entry = new LinkedList();
                entry.add(Sites.get(k).get(1));
                entry.add(Integer.valueOf(Sites.get(k).get(1)) + (MotifLength-1));
                mats.get(i).get(k).add(entry);
            }
        }
            
        //(2) Merge all extracted motifs into a single membership matrix
        // ----------------------------------------------------------------------- //

        // TODO: How does MM and mats structure in MATLAB look like ????
        //initialize: MM takes the first input set.
        List MM = new LinkedList(mats.get(0));

        int Motifs = ((List)MM.get(3)).size();
        int Rows = MM.size();
            
        //incrementally add additional motifs to the MM structure.
        for (int i = 1; i < mats.size(); i++) {

            Motifs++; //each "mat" group involves a new motif.

            for (int j = 0; j < mats.get(i).size(); j++) {
                boolean MakeNewRow = true;
                for (int k = 0; k < Rows; k++) {
                    if ((mats.get(i).get(j).get(0).equals(((List)MM.get(k)).get(0)))) {

                        //instead of making a new row, write the information to
                        //the appropriate place.
                        MakeNewRow = false;
                        while (((List)MM.get(k)).size() < Motifs) {
                            ((List)MM.get(k)).add(null);
                        }
                        ((List) MM.get(k)).add(Motifs-1,mats.get(i).get(j).get(1));
                    }
                }

                if (MakeNewRow) {
                    //add a row
                    Rows++;

                    //write information to new row.
                    LinkedList entry = new LinkedList();
                    entry.add(mats.get(i).get(j).get(0));
                    while (entry.size() < Motifs) {
                        entry.add(null);
                    }
                    entry.add(Motifs-1,mats.get(i).get(j).get(1));
                    MM.add(entry);
                }
            }
        }
            
        //(3) Sort sites by number of motifs contained
        // ----------------------------------------------------------------------- //
        // TODO check if sorting is correct as in original MATLAB version

        Collections.sort(MM,new Comparator<List>() {
            public int compare(List list, List list1) {
                return ((List)list.get(1)).size() - ((List)list1.get(1)).size();
            }
        });

        return MM;
    }


    //Create FP directory, compute clusters, families, and FPs, save struct
    public void ComputeFPs(String fpPath, Double frequencyThreshold) {
        //start message
        System.out.println("Processing request...");

        //create FP directory
        new File(fpPath).mkdirs();

        System.out.println("Clusters successfully retrieved.");
        System.out.println("Building FP profiles...");

        //obtain specifications for FP construction

        //build FPs
        MC_MakeFPs(seqFile,mkr.AllMotifsAndLocations, fpPath, frequencyThreshold, dataSetProfile);

        System.out.println("Job successfully completed!");
    }

    private void MC_MakeFPs(String seqFile, AllMotifsAndLocations allMotifsAndLocations, String fPdir, double sharedSitesThreshold, DataSetProfile dataSetProfile) {

        //Description: This function creates FPs from the GUI tree viewer.

        //start off by making a directory to contain all the FBP sequence.
        String FPseq_dir = fPdir + "/Sequences";
        new File(FPseq_dir).mkdirs();

        //obtain sites
        List<List<String>> Sites = MC_GetSharedSites(allMotifsAndLocations,sharedSitesThreshold);

        //create a .fasta file with key members from a given family.
        for (int i = 0; i < Sites.size(); i++) {

           //initialize a family file for each family
           String familyfile = FPseq_dir + "/Family" + i + ".fasta";

           //build each family file with the appropriate sequences
            for (int j = 0; j < Sites.get(i).size(); j++) {
                List<Sequence> seqs = null;
                try {
                    seqs = FastaTools.getSequences(new FileInputStream(seqFile), Sites.get(i));
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                }
                for (Sequence s : seqs) {
                    FastaTools.fastaAppend(familyfile,s.header,s.content);
                }
            }

            //find the best motif possible for each familyfile.
            //make an output directory
            String outputdir = fPdir + "/Family" + i;
            new File(outputdir).mkdirs();

            //map to correct bfile
            String bfile = dataSetProfile.memeParam.bfile;

            //input to meme program.
            //note that here we use 'oops' instead of 'zoops', this requires that
            //all input sites produce a motif.
            String cmd = dataSetProfile.progParam.MEME + " " + familyfile + " -minw " + dataSetProfile.memeParam.MinWidth + " -maxw " + dataSetProfile.memeParam.MaxWidth + " -" + dataSetProfile.memeParam.Alphabet + " -bfile " + bfile + " -mod oops -oc " + outputdir + " -nostatus -maxsize 1000000";

            //option: search the reverse complement, or not.
            if (dataSetProfile.memeParam.RevComp) {
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
        }
    }

    private List<List<String>> MC_GetSharedSites(AllMotifsAndLocations allMotifsAndLocations, double sharedSitesThreshold) {
        //This function creates a list of sequence entries from a MM-like structure
        //(AllMotifsAndLocations) that appear at least (FP frequency threshold) 
        //percent of the time among all runs that are included in a given motif group.
        
        List<List<String>> Sites = new LinkedList<List<String>>();

        //TODO!!!!!!
        /*
        for (int i = 0; i < allMotifsAndLocations.size(); i++) {
            //each family from a given data set has a membership matrix.
            int NumSites = 0;
            Sites.add(i, new LinkedList<String>());
            
            for (int j = 1; j < allMotifsAndLocations.get(i).size(); j++) {
                //each row
                int NonEmptySites = 0;
                for (int k = 1; k < allMotifsAndLocations.get(i).get(j).size(); k++) {
                    if (allMotifsAndLocations.get(i).get(j).get(k) != null) {
                        NonEmptySites++;
                    }
                }
                
                //if enough sites pass, note the regions
                //TODO: check what happens in original code
                if (NonEmptySites / (allMotifsAndLocations.get(i).get(0).size()-1) >= sharedSitesThreshold) {
                    Sites.get(i).add(allMotifsAndLocations.get(i).get(j).get(0));
                }
            }
        }
        */
        return Sites;
    }

    /*
    public List MakeMotifMap(String ethc, DirParam dirParam, String SeqFile, MEMEParam memeParam) {
        System.out.println("Processing request...");

        //retrieve FPdir

        // TODO: Why here and not in next step ???
        MakeMotifMap mmm = new MakeMotifMap();
        List mm = mmm.MC_MakeMotifMap(dirParam, SeqFile, memeParam.RevComp);

        System.out.println("Job successfully completed!");
        return mm;
    }
    */
}
