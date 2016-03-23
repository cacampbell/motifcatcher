package motifcatcher;
import motifcatcher.utils.FastaTools;
import org.apache.commons.lang3.exception.ExceptionUtils;
import java.io.*;
import java.util.*;

//TODO: Fix and Rewrite MC_MakeMotifMap
public class MakeMotifMap {

    private static void DebugException(Exception e) {
        Throwable cause = e.getCause();
        String[] exceptionTrace = ExceptionUtils.getRootCauseStackTrace(cause);
        String exceptionTraceString = Arrays.toString(exceptionTrace);
        System.out.print(exceptionTraceString);
    } //debugException

    public List MC_MakeMotifMap(DirParam dirParam, String seqFile, boolean revComp) {
        /*
        After clusters, families, etc have been determined on a given data set,
        this command generates the final, overall MotifMap.  a MotifMap is an
        object that describes significant motifs in the entire data set.

        Inputs:
        Directory
        Folder containing all information from previous analyses, but most
        importantly, meme output directories with FBP stuff.
        FBP (R, x) taken from each Family output file.
        SeqFile
        .fasta file containing all sequences in the superset.
        RevComp
        motifs may be found on (+) or (-) strand, or just (+) strand

        Outputs:
        MotifMap
        Output structure in the same form of "AllMotifsAndLocations",
        except instead of runs, labeled by Family.
        As per design in the "directory" folder, Families are automatically
        numbered by their rank of convergence (Family 1 describes related
        ES most often converged, Family 2 second-most-often, etc).
        row 1: Family, with Family Number
        row 2: regular expression
        row 3-}: sites (col 1), localization (other columns)
        build a list of all MEME output directories.
        */

        File dir = new File(dirParam.FP);
        File[] fileArray = dir.listFiles();

        //determine total number of clusters to extract
        int MaxFam = -1;
        assert fileArray != null;
        for (File file : fileArray) {
            //find the total number of families.
            if (file.getName().contains("Family")) {
                String[] num = file.getName().split("Family");
                if (num[1] != null) {
                    int maxN = Integer.valueOf(num[1]);
                    if (maxN > MaxFam) {
                        MaxFam = maxN;
                    }
                }
            }
        }
        MaxFam++;
        List motifMap = null;

        if (MaxFam > 0) {
            //Tried to generify, and changed names from "list1-list5" to canonical names of their expected data
            LinkedList<LinkedList<String>> Headers = new LinkedList<LinkedList<String>>();
            LinkedList<String> InputSet = new LinkedList<String>();
            LinkedList<String> familyNumber = new LinkedList<String>();
            LinkedList<String> motifRegularExpression = new LinkedList<String>();
            LinkedList<String> motif = new LinkedList<String>();
            familyNumber.add(null);
            motifRegularExpression.add(null);
            motif.add(null);

            for (int i = 0; i < MaxFam; i++) {
                InputSet.add(i, dirParam.FP + File.pathSeparator + "Family" + i);
            }

            for (int i = 0; i < MaxFam; i++) {
                String file = InputSet.get(i) + File.pathSeparator + "meme.txt";
                familyNumber.add(i, "Family" + i);
                motifRegularExpression.add(i, MC_GetRegularExpression(file));
                motif.add(i, MC_GetMotif(file)); // "TODO: Not a string .... !!!!"
                //Thanks Sven, that's super helpful
                //TODO: What the fuck is happening here?
            }

            Headers.add(familyNumber);
            Headers.add(motifRegularExpression);
            Headers.add(motif);
            //Recover the Motif Map
            motifMap = MC_MakeMembershipMatrixwHeader(InputSet, seqFile, Headers, revComp);
            //Add non-motif containing sequences to the 'structure' end
            FastaTools.addAllMissing(seqFile, motifMap, Headers.size());

        } else {
            System.err.println("No Motif Families were found in this directory");
        }

        return motifMap;
    } //MC_MakeMotifMap

    //TODO: Fix and Rewrite MC_MakeMembershipMatrixHeader
    //TODO: Rename / Refactor this method
    private static List MC_MakeMembershipMatrixwHeader(List<String> inputSet, String sequences, LinkedList<LinkedList<String>> headers, boolean revComp) {
        /*
        function [MM MotifLengths] =
        This function produces a membership matrix for cases where top motifs are
        determined by direct motif comparison (and clustering).

        Inputs:
        InputsSet:
        cell array of MEME output directories (whole path)
        SeqFileHeader:
        .fasta file of all sequences in the superset.
        Headers:
        Labels for each data set.

        Ouputs:
        MM:
        merged membership matrix
        MotifLengths:
        length of motifs, ordered as same as InputSet.

        Algorithm:
        (1) Extract all information from all directories in input set, and
        create a miniature membership matrix for each set.
        (2) Merge all sites from each directory into a single membership
        matrix.
        (3) sort membership matrix by multiplicity of sites included.
        (1) Extract all information each file in the InputSet
        ----------------------------------------------------------------------- //
        initalize extraction cell array, and MotifLengths.
        */

        /*
        * TODO: Fix Wonky List strucutres. Cell Arrays are not real.
        * TODO: Take a stab at improving the macro performance, this must be O(n^3) in at least two places
        * These horrible methods could be why the program is so slow with large data sets, even with threading.
        * I am going to try to make some sensible infrastructures for handling the data
        * It appears that Sven meant to literally replicate the MATLAB Cell
        * Array structure that was originally used to transport the data.
        * If necessary, this can be changed.
        *
        * CHRIS: Why is this happening to me???
        * SVEN: Trust me, I'm a doctor.
        */
        ArrayList<LinkedList<LinkedList<String>>> mats = new ArrayList<LinkedList<LinkedList<String>>>(inputSet.size());
        ArrayList<Integer> MotifLengths = new ArrayList<Integer>(inputSet.size());

        for (int i = 0; i < inputSet.size(); i++) {
            String file = inputSet.get(i) + File.pathSeparator + "meme.txt";
            mats.add(i, new LinkedList<LinkedList<String>>());
            LinkedList<String> extraction = new LinkedList<String>();
            BufferedReader br = null;

            /* I am not entirely sure whether or not it matters
             * that the file read operation fails for any one particular MEME
             * output file. I would guess that it doesn't, since the idea is to
             * combine the relatively useless subsets of MEME calls into a
             * more powerful meta motif. I think it is appropriate to continue
             * after exception then.
             */
            try {
                br = new BufferedReader(new FileReader(file));
                if (br.ready()) {
                    String line;
                    while ((line = br.readLine()) != null) {
                        extraction.add(line);
                    } //read in each line of MEME output

                    br.close();
                }
            } catch (FileNotFoundException e) { //Could not find MEME output
                DebugException(e);
            } catch (IOException e) { //IO error in MEME output operations
                DebugException(e);
            } finally {
                try {
                    if (br != null) { br.close(); }
                } catch (IOException e) {
                    DebugException(e);
                } // Catch IO Exception from attempting to close BufferedReader
            } //Finally, close the BufferedReader if open

            /*
            'Create a dummy Eval, in the case that no PWM exists for this motif.'
            Here, Sven says "In the case" without any surrounding logic
            Could this be improved, or was this just meant to be "in the eventuality"
            */
            Integer start = null;
            Integer stop = null;
            Integer MotifLength = 0;

            //'retrieve PWM coordinates from extraction.'
            /*
             * Surely, there is a better way to do this than what Sven has laid
             * out below.
             * TODO: Rewrite this section to inspect the extraction, then format PWM
             * Like, seriously, wtf
             * Maybe explain the format of MEME output in different cases here
             * TODO: Iterate over the LinkedList 'extraction', rather than looping
             * TODO: Move this to separate method, and reformat iteration
             * this is a war crime
             */
            for (int j = 0; j < extraction.size(); j++) {
                String line = extraction.get(j);
                String[] cols = line.split("\\s+"); // Splitting each line by whitespace (1 or more chars)

                //revcomp version
                if (revComp) { //Can move this check elsewhere
                    if ( cols.length >= 5 &&
                            cols[0].equals("Sequence") &&
                            cols[1].equals("name") &&
                            cols[2].equals("Strand") &&
                            cols[3].equals("Start") &&
                            cols[4].equals("P-value") &&
                            cols[5].equals("Site")) {

                        start = j + 2;
                    }
                } else { //not revcomp
                    if ( cols.length >= 4 &&
                            cols[0].equals("Sequence") &&
                            cols[1].equals("name") &&
                            cols[2].equals("Start") &&
                            cols[3].equals("P-value") &&
                            cols[4].equals("Site")) {

                        start = j + 2;
                    } // if
                } //else

                if ( cols.length >= 3 &&
                    cols[1].equals("Motif") &&
                    cols[2].equals("1") &&
                    cols[3].equals("block") &&
                    cols[4].equals("diagrams")) {

                        stop = j - 3;
                } //if

                if ( cols.length >= 4 &&
                    cols[0].equals("MOTIF") &&
                    cols[1].equals("1") &&
                    cols[2].equals("width") &&
                    cols[3].equals("=")) {

                    MotifLength = Integer.valueOf(cols[4]);
                    MotifLengths.add(i,MotifLength);
                } //if

            } // for j

           //Find all sites.  Look for key character ":" in the appropriate
           //range.
           // TODO: Add this to the parse PWM method -- or whatever
           LinkedList<LinkedList<String>> Sites = new LinkedList<LinkedList<String>>();

            for (int j = start; j < stop; j++) {
                String line = extraction.get(j);
                String[] cols = line.split("\\s+");

                if (  cols[0].contains("[") ||
                      cols[0].contains(":") ||
                      cols[0].contains("|") ||
                      cols[0].contains("iY") ) {
                    //the "iY" search key was added for a yeast analysis, "i" for intergenic, "Y" for yeast
                    //seriously, wtf, hard coding the species tags? Are you for real Sven?

                    LinkedList<String> entry = new LinkedList<String>();
                    entry.add(cols[0]); //Site name

                    if (revComp) {
                        entry.add(cols[2]); //Start
                    } else {
                        entry.add(cols[1]); //Start
                    }

                    Sites.add(entry);
                } //if
            } //for j

           //add the header to the top of the file.
            //TODO: Headers is a linkedlist, so ... change to iteration
           for (int j = 0; j < headers.size(); j++) {
               // TODO check
               //I'm unsure about what we are checking, but seems like a good idea
               LinkedList<String> entry = new LinkedList<String>();
               entry.add(null);
               entry.add(headers.get(j).get(i + 1));  //what
               mats.get(i).add(j, entry);
           } //for j

           /*
           Find the full name of each Site by searching for it in the
           SeqFile.

           Note that here, (+) and (-)-strandedness are not considered;
           this is because the interesting form of the motif might
           actually be in reverse-compliment orientation.

           Instead, each motif length is evaluated by forward and reverse
           compliment, at later stages of analysis.
           */
            int offset = headers.size();
            for (int k = 0; k < Sites.size(); k++) {
                mats.get(i).add(k + offset, new LinkedList<String>());
                String s = FastaTools.findNameFor(sequences, Sites.get(k).get(0));
                if (s != null) {
                    //record site
                    mats.get(i).get(k + offset).add(s);
                    //record span of motif
                    LinkedList<String> entry = new LinkedList<String>();
                    entry.add(Sites.get(k).get(1));
                    String string = String.valueOf(Integer.valueOf(Sites.get(k).get(1)) + (MotifLength - 1));
                    entry.add(string);
                    mats.get(i).add(k + offset, entry); //This may be a source of bugs in the future
                    //Sven originally had a TO DO here that read 'optimize correct .. ?'
                    //TODO clean up and check this section for logical errors
                } else {
                    System.err.println("Sequence name for 'Sites'[" + k + "] is null");
                } //else
            } //for k
        } //for i

        //(2) Merge all extracted motifs into a single membership matrix
        // ----------------------------------------------------------------------- //

        //initialize: MM takes the first input set.
        LinkedList<LinkedList<String>> MM = new LinkedList<LinkedList<String>>(mats.get(0));
        /*
         * I think that it may be a good idea to move all of the list massacres to one unholy plague class
         * that handles the data for just this method. At least I would have a better idea what I am doing
         * while fixing the logic.
         */
        int Motifs = MM.get(3).size(); //TODO Fix this magic, related to fixing weird list structure
        int Rows = MM.size();

        //incrementally add additional motifs to the MM structure.
        for (int i = 1; i < mats.size(); i++) { //already initialized (so it starts with 1 here)
            /*
            TODO: Again, what the fucking fuck? Using for and get with LinkedLists is much slower than iterating
            http://stackoverflow.com/questions/1879255/performance-of-traditional-for-loop-vs-iterator-foreach-in-java
            Change this section and all other for loops to Iterate!!!
            */

            Motifs++; //each "mat" group involves a new motif.

            //and, a new header.
            for (int j = 0; j < headers.size(); j++) {
                MM.get(j).add(mats.get(i).get(j).get(1)); //
            }


            for (int j = 3; j < mats.get(i).size(); j++) { //TODO: change to iteration
                boolean MakeNewRow = true;

                for (int k = 3; k < Rows; k++) { //TODO: change to iteration
                    String nextRowName = mats.get(i).get(j).get(0);
                    String existingRowName = MM.get(k).get(0);

                    if (nextRowName.equals(existingRowName)) { //no new row, that row exists
                        MakeNewRow = false;
                        LinkedList<String> entry = MM.get(k);

                        while (entry.size() < Motifs - 1) {
                            entry.add(null);
                        } //TODO: why are you initializing this linked list???

                        String motif = mats.get(i).get(j).get(1);
                        entry.add(Motifs - 1, motif); //TODO: Holy shit. Convert to Push, same below
                        break; //Breaks out of (for k) (??)
                    }
                }

                if (MakeNewRow) { //Row needs to be added
                    Rows++;
                    LinkedList<String> entry = new LinkedList<String>();
                    String name = mats.get(i).get(j).get(0);
                    String motif = mats.get(i).get(j).get(1);

                    while (entry.size() < Motifs - 1) {
                        entry.add(null);
                    } //TODO: Refector this code to iteration -- move 'make new row' logic upward in scope

                    entry.add(name);
                    entry.add(Motifs - 1, motif);
                    MM.add(entry);
                }
            }
        }

        //(3) Sort sites by number of motifs contained
        // ----------------------------------------------------------------------- //
        // TODO check if sorting is correct
        // This note is from Sven

        Collections.sort(MM.subList(3, MM.size() - 1), new Comparator<List>() {
            public int compare(List list1, List list2) {
                int cL1 = 0;
                int cL2 = 0;
                for (Object i : list1) {
                    if (i != null) cL1++;
                }
                for (Object i : list2) {
                    if (i != null) cL2++;
                }

                return cL2 - cL1;
            }
        });

        return MM;
}

    private static float[][] MC_GetMotif(String file) {
        /*
        function Motif = MC_GetMotif(file)
        This function extracts a motif from a given meme.txt file.

        Inputs:
        file:
        meme.txt outputfile

        Outputs:
        motif:
        double array of motif, nx4, n = length of motif.
        perform file extraction
        TODO: Should this method be returning a String? Should it be returning a float[][]?
        I am unsure what is *supposed* to be happening here, so it looks like its
        time to rewrite. I think that it is supposed to be a float, based on what it does, but it
        is use as if what it should be returning is a string. Should I maybe be converting internally?
        Is the wrong method being called?
        */

        //TODO Move read file boiler plate to its own method from the various other methods
        LinkedList<String> extraction = new LinkedList<String>();
        BufferedReader br = null;

        try {
            br = new BufferedReader(new FileReader(file));
            if (br.ready()) {
                String line;

                while ((line = br.readLine()) != null) {
                    extraction.add(line);
                } //while

                br.close();
            } //if br is ready
        } catch (FileNotFoundException e) {
            DebugException(e);
        } catch (IOException e) {
            DebugException(e);
        } finally {
            try {
                if (br != null) br.close();
            } catch (IOException e) {
                DebugException(e);
            } //catch
        } //finally

        int start = -1, stop = -1, width = -1; //Initialize with program breaking value
        //retrieve PWM coordinates from extraction.
        for (int j = 0; j < extraction.size(); j++) { //TODO: Convert to Iteration
            String line = extraction.get(j); //part of converting to iteration
            String[] cols = line.split("\\s+"); //Splitting on whitespace

            //revcomp version
            if (cols.length >= 5 &&
                    cols[1].equals("Motif") &&
                    cols[2].equals("1") &&
                    cols[3].equals("position-specific") &&
                    cols[4].equals("probability") &&
                    cols[5].equals("matrix")) {

                String cols2[] = extraction.get(j + 2).split(" ");
                width = Integer.valueOf(cols2[5]); //width of matrix
                start = j + 3;
                stop = start + width;
            }  //TODO: Use parse matrix method (to be written) -- do away with this copypasta
        }

        //initialize the motif matrix.
        float motif[][] = new float[width][4];

        for (int i = start; i < stop; i++) {
            String cols3[] = extraction.get(i).split("\\s+");
            //Write to output matrix
            for (int j = 0; j < 4; j++) {
                motif[i - start][j] = Float.valueOf(cols3[j + 1]);
            } //for j
        } //for i
        return motif;
    } //MC_GetMotif

    private static String MC_GetRegularExpression(String file) {

        //function Seq = MC_GetRegularExpression(file)
        //This function returns a regular expression string from a meme output file.

        //TODO Move read file boiler plate to its own method from the various other methods
        LinkedList<String> extraction = new LinkedList<String>();
        BufferedReader br = null;

        try {
            br = new BufferedReader(new FileReader(file));
            if (br.ready()) {
                String line;

                while ((line = br.readLine()) != null) {
                    extraction.add(line);
                } //while

                br.close();
            } //if br is ready
        } catch (FileNotFoundException e) {
            DebugException(e);
        } catch (IOException e) {
            DebugException(e);
        } finally {
            try {
                if (br != null) br.close();
            } catch (IOException e) {
                DebugException(e);
            } //catch
        } //finally

        String expression = "";

        //retrieve PWM coordinates from extraction.
        for (int j = 0; j < extraction.size(); j++) { //TODO: Convert to iteration
            String line = extraction.get(j); //part of conversion to iteration
            String[] cols = line.split("\\s+");

            //revcomp version
            if (cols.length >= 4 &&
                    cols[1].equals("Motif") &&
                    cols[2].equals("1") &&
                    cols[3].equals("regular") &&
                    cols[4].equals("expression")) {

                expression = extraction.get(j + 2);
            }
        }

        //"translate to common language" -sven
        String Seq = "";
        int j = 0;

        // TO DO replace the following program with regular expression substitution !!!
        //Original TO DO is listed above
        //TODO Replace the following with a regular expression builder
        //TODO make sense of just what the hell sven was thinking when he wrote this
        //TODO: find / write [ACGT]* to nucleic acid notation library

        while (j < expression.length()) {
           //Holy shit sven are you for real. THIS IS NOT HOW YOU BUILD A REGEX
           //check for regular nucleotide expressions
           if (expression.charAt(j) == 'A') {
               Seq += "A";
               j++;
           } else if (expression.charAt(j) == 'C') {
               Seq += "C";
               j++;
           } else if (expression.charAt(j) == 'G') {
               Seq += "G";
               j++;
           } else if (expression.charAt(j) == 'T') {
               Seq += "T";
               j++;
            } else {
               //evaluate nested expression
               int stop = expression.indexOf(']', j);

               if (stop - j - 1 == 2) { //2-ples
                   String frag = expression.substring(j + 1, stop);
                   if (frag.contains("A") && frag.contains("C")) {
                       Seq += "M";
                   } else if (frag.contains("A") && frag.contains("G")) {
                       Seq += "R";
                   } else if (frag.contains("A") && frag.contains("T")) {
                       Seq += "W";
                   } else if (frag.contains("C") && frag.contains("G")) {
                       Seq += "S";
                   } else if (frag.contains("C") && frag.contains("T")) {
                       Seq += "Y";
                   } else if (frag.contains("G") && frag.contains("T")) {
                       Seq += "K";
                   } //if-else translating motif fragments into DNA chars
                   j += 4;

                } else if (stop - j - 1 == 3) { //3-ples

                   String frag = expression.substring(j+1,stop);
                   if (frag.contains("C") && frag.contains("G") && frag.contains("T")) {
                       Seq += "B";
                   } else if (frag.contains("A") && frag.contains("G") && frag.contains("T")) {
                       Seq += "D";
                   } else if (frag.contains("A") && frag.contains("C") && frag.contains("T")) {
                       Seq += "H";
                   } else if (frag.contains("A") && frag.contains("C") && frag.contains("G")) {
                       Seq += "V";
                   }
                   j += 5;

                } else { //4-ples
                  Seq += "N";
                  j += 6;
               }
           }

       }

        //return value "Seq" is in the appropriate format
        return Seq;
    }
}
