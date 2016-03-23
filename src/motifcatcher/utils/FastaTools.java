package motifcatcher.utils;

import java.io.BufferedWriter;
import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: cobalt
 * Date: 04.01.2012
 * Time: 11:44
 * To change this template use File | Settings | File Templates.
 */
public class FastaTools {

    public static void fastaAppend(String filetitle, String s, String content) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filetitle, true));
            bw.write('>');
            bw.write(s);
            bw.newLine();
            bw.write(content);
            bw.newLine();
            bw.newLine();
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public static FastaInfo fastaGetInfo(String seqFile) {

        FastaInfo fi = new FastaInfo();
        try {
            BufferedReader br = new BufferedReader(new FileReader(seqFile));
            if (br != null) {
                String line;
                while ((line = br.readLine()) != null) {
                    if (line.startsWith(">")) {
                        fi.numSequences++;
                        line = br.readLine();
                        if (line != null) {
                            fi.numCharacters += line.trim().length();
                        }
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return fi;
    }

    public static List<Sequence> fastaRead(String seqFile) {

        List<Sequence> res = new LinkedList<Sequence>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(seqFile));
            if (br != null) {
                String line;
                String oldH = null;
                StringBuffer seq = new StringBuffer();
                while ((line = br.readLine()) != null) {
                    if (line.startsWith(">")) {
                        //close previous
                        if (oldH != null && seq.length() != 0) {
                            Sequence s = new Sequence();
                            s.header = oldH;
                            s.content = seq.toString();
                            res.add(s);
                        }
                        // start new
                        oldH = line.substring(1);
                        seq = new StringBuffer();
                    } else {
                        seq.append(line.trim());
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return res;
    }

    public static Map<String,Sequence> getSequences(InputStream dis, HashSet<Integer> j) {
        HashMap<String,Sequence> res = new HashMap<String, Sequence>();
        BufferedReader br = new BufferedReader(new InputStreamReader(dis));
        int counter = 0;
        try {
            String line;
            String oldH = null;
            StringBuffer seq = new StringBuffer();
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    counter++;
                    //close previous
                    if (oldH != null && seq.length() != 0) {
                        if (j.contains(counter - 1)) {
                            Sequence s = new Sequence();
                            s.header = oldH;
                            s.content = seq.toString();
                            res.put(s.header,s);
                            if (res.size() == j.size()) {
                                break;
                            }
                        }
                    }
                    // start new
                    oldH = line.substring(1);
                    seq = new StringBuffer();
                } else {
                    seq.append(line.trim());
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return res;
    }

    public static List<String> completeHeaders(InputStream dis, HashSet<String> j) {
        List<String> res = new LinkedList<String>();
        BufferedReader br = new BufferedReader(new InputStreamReader(dis));
        try {
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    String sh = line.substring(1);
                    for (String x : j) {
                        //TODO better match strategy
                        if (x.equals(sh)) {
                            res.add(sh);
                            break;
                        }
                    }
                    if (res.size() == j.size()) {
                        break;
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return res;
    }

    public static List<Sequence> getSequences(InputStream stream, List<String> lookup) {
        List<Sequence> res = new LinkedList<Sequence>();
        BufferedReader br = new BufferedReader(new InputStreamReader(stream));
        try {
            String line;
            String oldH = null;
            StringBuffer seq = new StringBuffer();
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    //close previous
                    if (oldH != null && seq.length() != 0) {
                        for (String part : lookup) {
                            //TODO better match strategy
                            if (oldH.equals(part)) {
                                Sequence s = new Sequence();
                                s.header = oldH;
                                s.content = seq.toString();
                                res.add(s);
                            }
                        }
                    }
                    // start new
                    oldH = line.substring(1);
                    seq = new StringBuffer();
                } else {
                    seq.append(line.trim());
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return res;
    }

    public static String findNameFor(String seqFile, String s) {
        BufferedReader br = null;
        String res = null;
        try {
            br = new BufferedReader(new FileReader(seqFile));
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    String sh = line.substring(1);
                    //TODO better match strategy
                    if (sh.equals(s)) {
                        res = sh;
                        break;
                    }
                }
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return res;
    }

    public static void addAllMissing(String seqfile, List motifMap, int size) {

        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(seqfile));
            String line;
            while ((line = br.readLine()) != null) {
                if (line.startsWith(">")) {
                    String sh = line.substring(1);
                    //default: add all sequences
                    boolean DoNotAdd = false;
                    //if the sequence is already written in the structure, however, no need
                    //to add it.
                    for (int j = size; j < motifMap.size(); j++) {
                        if (((String)((List)motifMap.get(j)).get(0)).contains(sh)) {
                            DoNotAdd = true;
                            break;
                        }
                    }
                    //write to structure.
                    if (!DoNotAdd) {
                        List entry = new LinkedList();
                        entry.add(sh);
                        motifMap.add(entry);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
