package motifcatcher;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

//TODO: Phylosort is outdated and unmaintained; find different solution
//import phylosort.NewickParser;
//import phylosort.TreeNode;

//TODO: Oh boy, rewrite PhyTree.
//Dendrosort is a recent package that purports to parse phylogenetic trees
//Write my own Newick Parser?
public class PhyTree {

    public String filename = null;
    public TreeNode root;
    private ClusterRes res;
    private Double cacheCT;

    public PhyTree(String filename) {
        this.filename = filename;

        //TODO: Read this file in a more robust way
        String ns = "";
        try {
            BufferedReader b = new BufferedReader(new FileReader(filename));
            while (b.ready()) {
                ns += b.readLine();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        root = NewickParser.parse(ns);
    }

    //TODO: Okay, not sure what happened here. Necessary for the API?
    public PhyTree() {
    }

    //TODO: This too will need to be rewritten with whichever new package used
    public ClusterRes cluster(double clusteringThreshold) {

        if (cacheCT != null && cacheCT == clusteringThreshold) {
            return res;
        } else {
            res = new ClusterRes();
            res.LeafClusters = realCluster(root,clusteringThreshold);
            cacheCT = clusteringThreshold;
        }
        return res;
    }

    private List<Integer> realCluster(TreeNode tr, Double v) {
        /*
        function [clus,nclus,steps] = cluster(tr,v,varargin)
        %CLUSTER Construct clusters from a phylogenetic treeFile.
        %
        % I = CLUSTER(T,V) returns a column vector with cluster indices for each
        % species (leaf) of the phylogenetic treeFile object T. CLUSTER tests for every
        % k (number of possible clusters) all feasible treeFile partitions and chooses
        % the one that optimizes a predefined criterion. Then, CLUSTER selects the k-th
        % partition with the optimal value of the criterion. CLUSTER starts with
        % two clusters (k=2), and continues splitting the treeFile until the criterion
        % satisfies the scalar threshold value V, or the maximum number of possible
        % clusters is reached.
        %
        % CLUSTER(...,'CRITERION',C) sets the criterion used for determining the
        % number of clusters as a function of the species pairwise distances. The
        % available options are:
        %
        %   'maximum'    - Maximum within cluster pairwise distance (Wmax). Smaller
        %   (default)      is better, cluster splitting stops when Wmax<=V.
        %
        %   'median'     - Median within cluster pairwise distance (Wmed). Smaller
        %                  is better, cluster splitting stops when Wmed<=V.
        %
        %   'average'    - Average within cluster pairwise distance (Wavg). Smaller
        %                  is better, cluster splitting stops when Wavg<=V.
        %
        %   'ratio'      - Between/within cluster pairwise distance ratio defined
        %                  as:
        %                        BWrat = (trace(B)/(k-1)) / (trace(W)/(n-k))
        %
        %                  where B and W are the between/within-scatter matrices, k
        %                  is number of clusters, and n is number of species in the
        %                  treeFile. Larger is better, cluster splitting stops when
        %                  BWrat>=V.
        %
        %   'gain'       - Within cluster pairwise distance gain defined as:
        %
        %                        Wgain = (trace(Wold)/trace(W)-1) * (n-k-1)
        %
        %                  where W and Wold are the within-scatter matrices for k
        %                  and k-1, k is number of clusters, and n is number of
        %                  species in the treeFile. Larger is better, cluster splitting
        %                  stops when Wgain<=V.
        %
        %   'silhouette' - Average silhouette width (SWavg). SWavg ranges from -1
        %                  to +1, larger is better, cluster splitting stops when
        %                  SWavg>=V.
        %
        % CLUSTER(...,'MAXCLUST',N) sets the maximum number of possible clusters
        % for the tested partitions. N defaults to number of leaves in the treeFile.
        %
        % Notes:
        % 1. When using the 'ratio', 'gain', or 'silhouette' criterions, it is hard
        % to know a good threshold value V in advance. Set V to [] (empty) to find
        % the optimal number of clusters below N. Set also N to a small value to
        % avoid expensive computation by testing all possible number of clusters.
        % 2. When using the 'maximum', 'median', or 'average' criterions, the
        % optimal number of clusters is equal to N because such metric monotonically
        % decrease as k increases. Set V to [] (empty) to force CLUSTER to return N
        % clusters.
        %
        % CLUSTER(...,'DISTANCES',D) substitutes the patristic (treeFile) distances
        % with a user provided pairwise distance matrix D. For example, D can be
        % the real sample pairwise distances.
        %
        % [I,J] = CLUSTER(...) returns cluster indices for all the nodes in the
        % treeFile, leaf nodes and branch nodes.
        %
        % [I,J,S] = CLUSTER(...) returns the branch being considered for the cut
        % and the value of criterion at each step of the algorithm. Set V to []
        % (empty) and N to its default to obtain the whole curve of the criterion
        % versus the number of clusters. Observe that the computation of some
        % criterions may be computationally intensive.

        % Copyright 2009 The MathWorks, Inc.
        % $Revision: 1.1.6.4 $Author: batserve $ $Date: 2009/05/07 18:15:49 $

        % References:
        %
        % [1] Dudoit S, Fridlyan J, A prediction-based resampling method for
        %     estimating the number of clusters in a dataset. Genome Biology,
        %     3(7):research0036.10036.21, 2002.
        %
        % [2] Theodoridis, S. and Koutroumbas, K. Pattern Recognition, Academic
        %     Press, pp. 434-435, 1999.
        %
        % [3] Kaufman L, Rousseeuw PJ, Finding Groups in Data: An Introduction to
        %     Cluster Analysis. New York, Wiley, 1990.
        %
        % [4] Calinski R, Harabasz J, A dendrite method for cluster analysis.
        %     Commun Statistics, 3:1-27, 1974.
        %
        % [5] Hartigan JA, Statistical theory in clustering, J Classification,
        %     2:63-76, 1985.
        */


        int numLeaves = tr.getLeaves().size();
        int numBranches = numLeaves - 1;
        int numLabels = numBranches + numLeaves;


        //translate treeFile into matlab style
        PhyTreeTranslator ptt = new PhyTreeTranslator(tr);
        int[][] tree = ptt.getTree();
        double[] treeDist = ptt.getTreeDist();
        char criteria = 'm'; // TODO default: maximum

        //build P (pairwise distance matrix)
        //double P [][] = new double[numLeaves][numLeaves];
        double P [][] = ptt.getPDist();
        int n = Integer.MAX_VALUE;
        // TODO [criteria,P,n] = parse_inputs(tr,numLeaves,varargin{:});

        // Check v input and initialize properly in case v=[]
        if (v == null) {
            if (criteria == 's' || criteria == 'r') {
                v = Double.POSITIVE_INFINITY;
            } else {
                v = Double.NEGATIVE_INFINITY;
            }
        }
        if (criteria == 's'  || criteria == 'r') {
            v = -v;
        }

        // Find all possible binary clusterizations based on the input treeFile
        boolean bc[][] = new boolean[numLabels][numLabels];
        for (int i = 0; i < numLeaves; i++) {
            bc[i][i] = true;
        }
        for (int i = 0; i < numBranches; i++) {
            //TODO: any
            int col1 = tree[i][0]-1;
            int col2 = tree[i][1]-1;
            for (int r = 0; r < numLabels; r++) {
                bc[r][i+numLeaves] = bc[r][col1] | bc[r][col2];
            }
            //TODO: end
            bc[i+numLeaves][i+numLeaves] = true;
        }
        printM(bc);

        int selLeaves[] = new int[numLabels];
        for (int i = 0; i < numLeaves; i++) {
            selLeaves[i] = 1;
        }

        // Initialize all species and edges to the first cluster
        int eclus[] = new int[numLabels];
        Arrays.fill(eclus, 1);

        // Initialize output "steps"
        double steps[][] = new double[numLeaves][2];
        switch (criteria) { // for k = 1
            case 'm':
                steps[0][0] = Double.NaN;
                steps[0][1] = Collections.max(squareform(P));
                break;
            /*
            case 'd':
                steps(1,:) = [NaN median(squareform(P))];
                break;
            case 'a':
                steps(1,:) = [NaN mean(squareform(P))];
                break;
            case 's':
                steps(1,:) = [NaN 1];
                break;
            case 'r':
                steps(1,:) = [NaN 0];
                P2 = Math.P.^2;
                TP2 = sum(P2(:));
                break;
            case 'g':
                P2 = P.^2;
                Wold = sum(squareform(P2))./numLeaves;
                steps(1,:) = [NaN inf];
                break;
            */
        }

        // Initialize bcsel; used for marking edges that had been removed for
        // creating partitions. bcsel(end) is the root and is invalidated from the
        // beginning, no partition is valid at the root.
        boolean bcsel[] = new boolean[numLabels];
        bcsel[bcsel.length-1] = true;

        List<Double> cache = new ArrayList<Double>();
        cache.add(0.0);
        int k = 1;
        // Main loop: repeat for every k>1 until exit conditions apply
        while ((k < numLeaves) && (k<n) && (steps[k-1][1] > v)) {

            k++;
            List<Integer> pospart = findNot(bcsel);
            double mas[] = new double[pospart.size()];
            double nma[][] = new double[pospart.size()][2];


            double[][] oma = matMul(matNE(newMat(1, k - 1),subMat(eclus,pospart)),cache);


            switch (criteria) {
                case 'm':
                    for (int i = 0;i < pospart.size(); i++) {
                        int[] t1 = subMatRow(eclus,selLeaves);
                        int t2 = subMatRow(eclus, pospart.get(i));
                        int[] thisclu = subMatEqual(t1, t2);
                        boolean[] b1 = subMat(bc, selLeaves, pospart.get(i));
                        int[] idx = matAnd(thisclu, b1);
                        if (sum(idx)>1) {
                            double[][] p = subMat(P, idx, idx);
                            List<Double> squareform = squareform(p);
                            nma[i][0] = Collections.max(squareform);
                        }
                        b1 = subMat(bc, selLeaves, pospart.get(i));
                        idx = matAnd(thisclu, matNot(b1));
                        if (sum(idx)>1) {
                            double[][] p = subMat(P, idx, idx);
                            List<Double> squareform = squareform(p);
                            nma[i][1] = Collections.max(squareform);
                        }
                    }
                    printM(nma);
                    mas = maxInCol(oma,nma);
                    /*
                case "d"
                    for i = 1:numel(pospart)
                        thisclu = eclus(selLeaves)==eclus(pospart(i));
                        idx = thisclu & bc(selLeaves,pospart(i));
                        if sum(idx)>1
                            nma(i,1) = median(squareform(P(idx,idx)));
                        }
                        idx = thisclu & ~bc(selLeaves,pospart(i));
                        if sum(idx)>1
                            nma(i,2) = median(squareform(P(idx,idx)));
                        }

                    }
                    mas = max([oma nma],[],2);
                case "a"
                    for i = 1:numel(pospart)
                        thisclu = eclus(selLeaves)==eclus(pospart(i));
                        idx = thisclu & bc(selLeaves,pospart(i));
                        if sum(idx)>1
                            nma(i,1) = mean(squareform(P(idx,idx)));
                        }
                        idx = thisclu & ~bc(selLeaves,pospart(i));
                        if sum(idx)>1
                            nma(i,2) = mean(squareform(P(idx,idx)));
                        }

                    }
                    mas = max([oma nma],[],2);
                case "s" // Reference [3] (we use the negative of the criterion since
                           // this implementation searches for the cut that
                           // minimizes the criterion while the reference indicates
                           // that larger values of the silhouette  are better)
                    for i = 1:numel(pospart)
                        tclus = eclus(selLeaves) + k.*bc(selLeaves,pospart(i));
                        sili = silhouette((1:numLeaves)",tclus,@(x,y) P(x,y));
                        mas(i) = - mean(sili);
                    }
                case "r" // Reference [4] (we use the negative of the criterion since
                           // this implementation searches for the cut that
                           // minimizes the criterion while the reference indicates
                           // that larger values of the ratio are better)
                    for i = 1:numel(pospart)
                        thisclu = eclus(selLeaves)==eclus(pospart(i));
                        idx = thisclu & bc(selLeaves,pospart(i));
                        if sum(idx)>1
                            nma(i,1) = sum(squareform(P2(idx,idx)))./sum(idx);
                        }
                        idx = thisclu & ~bc(selLeaves,pospart(i));
                        if sum(idx)>1
                            nma(i,2) = sum(squareform(P2(idx,idx)))./sum(idx);
                        }
                    }
                    mas = (1-TP2./sum([oma nma],2)) .* (numLeaves-k) ./ (k-1);
                    mas(isnan(mas)) = -inf;
                case "g" // Reference [5]
                    for i = 1:numel(pospart)
                        thisclu = eclus(selLeaves)==eclus(pospart(i));
                        idx = thisclu & bc(selLeaves,pospart(i));
                        if sum(idx)>1
                            nma(i,1) = sum(squareform(P2(idx,idx)))./sum(idx);
                        }
                        idx = thisclu & ~bc(selLeaves,pospart(i));
                        if sum(idx)>1
                            nma(i,2) = sum(squareform(P2(idx,idx)))./sum(idx);
                        }
                    }
                    mas = (numLeaves-k-1).*(Wold./sum([oma nma],2)-1);
                    */
            }

            // Find the cut that leads to the minimization of the criterion
            double bcri;
            //if (criteria == 'g') {
            //    bcri = max(mas);
            //} else {
            bcri = min(mas);
            //}
            int mins[] = subMatEqlMat(pospart, mas, bcri);
            int bidx;
            if (mins.length > 1) {
                // From all the possible cuts that give the minimum criterion select
                // the longest edge for the cut:
                bidx = maxPos(treeDist, mins);
            } else {
                bidx = mins[0];
            }
            // Find cluster being cut:
            int kk = eclus[bidx];
            // Create new cluster:
            for (int x = 0; x < eclus.length; x++) {
                if (eclus[x] == kk && bc[x][bidx]) {
                    eclus[x] = k;
                }
            }
            // Update pre-calculated criteria on current clusters
            int index = -1;
            for (int x=0; x<pospart.size(); x++) {
                if (pospart.get(x) == bidx) {
                    index = x;
                    break;
                }
            }
            int maxindex = Math.max(k,kk) - 1;
            while (maxindex >= cache.size()) {
                cache.add(0.0);
            }
            cache.set(k-1, nma[index][0]);
            cache.set(kk-1,nma[index][1]);

            // In case of criteria is "gain" recalculate the best within-scatter
            // matrix:
            //if (criteria == 'g') {
            //    Wold = sum([oma(pospart==bidx,:) nma(pospart==bidx,:)]);
            //}

            // Mark edge being cut
            bcsel[bidx] = true;
            // Remove edges in branches where the other two adjacent edges had
            // already been cut:
            boolean[] tbcsel = matNot(bcsel);
            while (matNE(tbcsel,bcsel)) {
                tbcsel = bcsel;

                int[][] ints1 = matSelElm(tree, bcsel);
                int[] ints2 = subMat(bcsel, numLeaves, bcsel.length);
                int[][] cats = matCatHor(ints1, ints2);
                int[] ints = matSumRow(cats);
                List<Integer> q = matEql(ints,2);
                Set<Integer> nodes = find(tree,q);
                for (Integer x : nodes)
                    bcsel[x-1] = true;
                for (int x : q)
                    bcsel[numLeaves+x] = true;
            }
            // Update output "steps"
            steps[k-1][0] = bidx;
            steps[k-1][1] = bcri;
        }

        steps = Arrays.copyOfRange(steps,0,k);
        //if criteria == "g" {
        //    [~,optk] = max(steps(2:end,2));
        //    optk = optk +1;
        //} else {
        int optk = minPos2(steps);
        //}
        optk = Math.max(2,optk);

        // calculate cluster for each node with propagation towards leaves
        int[] nclus = new int [numLabels];
        Arrays.fill(nclus,1);
        double[][] cutstemp = Arrays.copyOfRange(steps,1,optk); //TODO ???????
        Double[] cuts = new Double[cutstemp.length];
        for (int x = 0; x < cuts.length; x++) {
            try {
                cuts[x] = cutstemp[x][0];
            } catch (Exception e) {
                e.printStackTrace();
                cuts[x] = 0.0;
            }
        }
        Arrays.sort(cuts,new Comparator<Double>() {
            public int compare(Double doubles, Double doubles1) {
                if (doubles.equals(doubles1)) {
                    return 0;
                } else if (doubles < doubles1) {
                    return 1;
                } else if (doubles > doubles1) {
                    return -1;
                }
                return 0;
            }
        });

        for (int i = 0; i < cuts.length; i++) {
            k = i+2;
            nclus[(int)Math.round(cuts[i])] = k;
            for (int j = (int)Math.round(cuts[i]) - numLeaves; j >= 0; j--) {
                if (nclus[j+numLeaves] == k) {
                    nclus[tree[j][0]-1] = k;
                    nclus[tree[j][1]-1] = k;
                }
            }
        }

        // first output argument is only cluster index for each leaf
        int[] clus = subMatRow2(nclus,newMat(0,numLeaves-1));

        //if (criteria == 's' || criteria == 'r') {
        //    for (int x = 0; x < steps.length; x++) {
        //        steps[x][2] = -steps[x][2];
        //    }
        //}


        List<Integer> res = new ArrayList<Integer>(clus.length);
        for (int i : clus) {
            res.add(i);
        }
        return res;

    }

    public class PhyTreeTranslator {

        TreeNode tr;
        private int[][] tree;
        private double[] treeDist;
        private double[] rootDist;
        private double D[][];
        private boolean isDcomputed = false;
        private HashMap<TreeNode,Integer> map = new HashMap<TreeNode,Integer>();
        int cN;
        int cL;
        int cB;
        int numLeaves;
        int numLabels;

        public PhyTreeTranslator(TreeNode tr) {

            this.tr = tr;

            numLeaves = tr.getLeaves().size();
            int numBranches = numLeaves - 1;
            numLabels = numBranches + numLeaves;

            cL = 0;
            cB = numLeaves;
            cN = 1;

            tree = new int[numBranches][2];
            treeDist = new double[numLabels];
            rootDist = new double[numLabels];

            D = new double[numLabels][numLabels];
            for (int y = 0; y < numLabels; y++) {
                for (int x = 0; x < numLabels; x++) {
                    if (x == y)
                        D[y][x] = 0;
                    else
                        D[y][x] = Double.POSITIVE_INFINITY;
                }
            }
            translate();
        }

        public int[][] getTree() {
            return tree;
        }

        public double[] getTreeDist() {
            return treeDist;
        }

        private void translate() {
            translateTree(tr,0);
        }

        private int translateTree(TreeNode tr, double dist) {
            if (!tr.isLeaf()) {
                int count = 0;
                int[] bx = new int[tr.getChildrenCount()];
                double[] tx = new double[tr.getChildrenCount()];
                for (TreeNode t : tr.getChildren()) {
                    int b = translateTree(t,dist+tr.getLength());
                    bx[count] = b;
                    tx[count++] = t.getLength();
                }
                //initialize D
                for (int i = 0; i < bx.length; i++) {
                    D[cB][bx[i]-1] = tx[i];
                    D[bx[i]-1][cB] = tx[i];
                }
                tree[cB-numLeaves] = bx;
                treeDist[cB] = tr.getLength();
                rootDist[cB] = tr.getLength() + dist;
                map.put(tr,cB);
                cB++;
                return cB;
            } else {
                treeDist[cL] = tr.getLength();
                rootDist[cL] = tr.getLength() + dist;
                map.put(tr,cL);
                cL++;
                return cL;
            }
        }

        public double[][] getPDist() {
            if (!isDcomputed) {
                for (int k = 0; k < numLabels; k++) {
                    for (int i = 0; i < numLabels; i++) {
                        for (int j = 0; j < numLabels; j++) {
                            if (D[i][k] + D[k][j] < D[i][j]) {
                                D[i][j] = D[i][k] + D[k][j];
                            }
                        }
                    }
                }
                isDcomputed = true;
            }
            double[][] res = new double[numLeaves][numLeaves];
            for (int i = 0 ; i < numLeaves; i++) {
                res[i] = Arrays.copyOf(D[i],numLeaves);
            }
            printM(res);
            return res;
        }
    }


    private Set<Integer> find(int[][] tree, List<Integer> q) {
        HashSet<Integer> res = new HashSet<Integer>();
        for (Integer y : q) {
            for (int i = 0; i < tree[y].length; i++) {
                res.add(tree[y][i]);
            }
        }
        return res;
    }

    private List<Integer> matEql(int[] ints, int i) {
        List<Integer> res = new ArrayList<Integer>(ints.length);
        for (int y = 0; y < ints.length; y++) {
            if (ints[y] == i) {
                res.add(y);
            }
        }
        return res;
    }

    private int[] matSumRow(int[][] ints) {
        int[] res = new int[ints.length];
        for (int y = 0; y < ints.length; y++) {
            for (int x = 0; x < ints[0].length; x++) {
                res[y] += ints[y][x];
            }
        }
        return res;
    }

    private int[][] matCatHor(int[][] ints1, int[] ints2) {
        int res[][] = new int[ints1.length][ints1[0].length+1];
        for (int y = 0; y < ints1.length; y++) {
            System.arraycopy(ints1[y], 0, res[y], 0, ints1[0].length);
            res[y][ints1[0].length] = ints2[y];
        }
        return res;
    }

    private int[][] matSelElm(int[][] tree, boolean[] bcsel) {
        int[][] res = new int[tree.length][tree[0].length];
        for (int y = 0; y < res.length; y++) {
            for (int x = 0; x < res[0].length; x++) {
                if (tree[y][x] < bcsel.length && bcsel[tree[y][x]-1]) {
                    res[y][x] = 1;
                }
            }
        }
        return res;
    }

    private int[] subMat(boolean[] bcsel, int numLeaves, int i) {
        int[] res = new int[i-numLeaves];
        for (int x = numLeaves; x < i; x++) {
            if (bcsel[x])
                res[x-numLeaves] = 1;
            else
                res[x-numLeaves] = 0;
        }
        return res;
    }

    private int minPos2(double[][] steps) {

        //double res = Double.NaN;
        int pos = 2;
        double min = Double.MAX_VALUE;
        for (int i = 0; i < steps.length; i++) {
            if (steps[i][1] < min) {
                min = steps[i][1];
                //res = steps[i][0];
                pos = i+1;
            }
        }
        //return (int)Math.round(res);
        return pos;
    }

    private boolean matNE(boolean[] tbcsel, boolean[] bcsel) {
        for (int i = 0; i < tbcsel.length; i++) {
            if (tbcsel[i] != bcsel[i]) {
                return true;
            }
        }
        return false;
    }

    private int maxPos(double[] treeDist, int[] mins) {

        int maxPos = -1;
        double max = Double.MIN_VALUE;

        for (int x : mins) {
            if (treeDist[x] > max) {
                maxPos = x;
                max = treeDist[x];
            }
        }

        return maxPos;
    }

    private double min(double[] mas) {
        double res = Double.MAX_VALUE;
        for (double v : mas) {
            if (v < res) res = v;
        }
        return res;
    }

    private int[] subMatEqlMat(List<Integer> pospart, double[] mas, double bcri) {
        List<Integer> r = new ArrayList<Integer>();
        for (int i = 0; i < pospart.size(); i++) {
            if (mas[i] == bcri) r.add(pospart.get(i));
        }
        int[] res = new int[r.size()];
        for (int i = 0; i < r.size(); i++) {
            res[i] = r.get(i);
        }
        return res;
    }


    private double[] maxInCol(double[][] oma, double[][] nma) {
        double [] res = new double[oma.length];
        for (int y = 0; y < oma.length; y++) {
            res[y] = Double.MIN_VALUE;
            for (int x = 0; x < oma[y].length; x++) {
                if (res[y] < oma[y][x]) res[y] = oma[y][x];
            }
        }
        for (int y = 0; y < nma.length; y++) {
            for (int x = 0; x < nma[y].length; x++) {
                if (res[y] < nma[y][x]) res[y] = nma[y][x];
            }
        }
        return res;
    }

    private double[][] matMul(int[][] ints, List<Double> cache) {
        double [][] res = new double[ints.length][cache.size()];
        for (int y = 0; y < ints.length; y++) {
            for (int x = 0; x < cache.size(); x++) {
                res[y][x] = ints[y][x] * cache.get(x);
            }
        }
        return res;
    }

    private int[][] matNE(int[] ints, int[] ints1) {
        int[][] res = new int[ints1.length][ints.length];
        for (int x = 0; x < ints.length; x++) {
            for (int y = 0; y < ints1.length; y++) {
                if (ints1[y] != ints[x]) {
                    res[y][x] = 1;
                } else {
                    res[y][x] = 0;
                }
            }
        }
        return res;

    }

    private int[] newMat(int i, int i1) {
        int[] res = new int[i1-i+1];
        int c = 0;
        for (int x = i; x <= i1; x++) {
            res[c++] = x;
        }
        return res;
    }

    private int[] subMat(int[] eclus, List<Integer> pospart) {
        int[] res = new int[pospart.size()];
        int c = 0;
        for (Integer i : pospart) {
            res[c++] = eclus[i];
        }
        return res;
    }

    private boolean[] matNot(boolean[] booleans) {
        boolean [] res = new boolean[booleans.length];
        for (int i = 0; i < res.length; i++) {
            res[i] = !booleans[i];
        }
        return res;
    }

    private int[] matAnd(int[] thisclu, boolean[] booleans) {
        int[] r = new int[thisclu.length];
        for (int x = 0; x < thisclu.length; x++) {
            if (booleans[x] && thisclu[x] != 0) {
                r[x] = thisclu[x];
            }
        }
        return r;
    }


    private boolean[] subMat(boolean[][] p, int[] idx, int idx1) {

        List<Boolean> r = new ArrayList<Boolean>(idx.length);
        for (int y = 0; y < idx.length; y++) {
            if (idx[y] == 1) {
                r.add(p[y][idx1]);
            }
        }
        boolean res[] = new boolean[r.size()];
        int cx = 0;
        for (Boolean aR : r) {
            res[cx++] = aR;
        }
        return res;
    }


    private double[][] subMat(double[][] p, int[] idx, int[] idx1) {
        List<List<Double>> ry = new ArrayList<List<Double>>(idx.length);
        for (int y = 0; y < p.length && y < idx.length; y++) {
            if (idx[y] == 0) continue;
            List<Double> rx = new ArrayList<Double>(idx1.length);
            for (int x = 0; x < p[y].length && x < idx1.length; x++) {
                if (idx1[x] == 0) continue;
                rx.add(p[y][x]);
            }
            ry.add(rx);
        }
        double res[][] = new double[ry.size()][ry.get(0).size()];
        int cx = 0;
        for (List<Double> y : ry) {
            int cy = 0;
            for (Double x : y) {
                res[cx][cy++] = x;
            }
            cx++;
        }
        return res;
    }

    private int sum(int[] idx) {
        int res = 0;
        for (int i : idx) {
            res += i;
        }
        return res;
    }

    private int[] subMatEqual(int[] ints, int i) {
        List<Integer> res = new ArrayList<Integer>();
        for (int x : ints) {
            if (x == i)
                res.add(x);
        }
        int[] r = new int[res.size()];
        int c = 0;
        for (Integer x : res) {
            r[c++] = x;
        }
        return r;
    }

    private int[] subMatRow(int[] eclus, int[] selLeaves) {
        List<Integer> r = new ArrayList<Integer>(selLeaves.length);
        for (int i = 0; i < selLeaves.length; i++) {
            if (selLeaves[i] == 1) {
                r.add(eclus[i]);
            }
        }
        int[] res = new int[r.size()];
        int count = 0;
        for (Integer i : r) {
            res[count++] = i;
        }
        return res;
    }

    private int[] subMatRow2(int[] eclus, int[] selrows) {
        List<Integer> r = new ArrayList<Integer>(selrows.length);
        for (int i = 0; i < selrows.length; i++) {
            r.add(eclus[selrows[i]]);
        }
        int[] res = new int[r.size()];
        int count = 0;
        for (Integer i : r) {
            res[count++] = i;
        }
        return res;
    }

    private int subMatRow(int[] eclus, int row) {
        return eclus[row];
    }

    private List<Integer> findNot(boolean[] bcsel) {
        List<Integer> res = new ArrayList<Integer>(bcsel.length);

        for (int i = 0; i < bcsel.length; i++) {
            if (!bcsel[i]) res.add(i);
        }

        return res;
    }

    private List<Double> squareform(double[][] p) {
        int size = p.length * (p.length-1) / 2;
        List<Double> res = new ArrayList<Double>(size);

        for (int y = 0; y < p.length-1; y++) {
            for (int x = y+1; x < p.length; x++) {
                res.add(p[x][y]);
            }
        }
        return res;
    }

    /*

        P = squareform(pdist(tr)); // Pairwise distances, stored in square form
                                   // since it will be indexed by row and columns
    */

    public PhyTreeInfo getPhyTreeInfo() {

        PhyTreeInfo i = new PhyTreeInfo();
        i.LeafNames = new ArrayList<String>();
        for (TreeNode n : root.getLeaves()) {
            i.LeafNames.add(n.getLabel());
        }
        return i;
    }


    public class ClusterRes {
        public List<Integer> LeafClusters;
        List<Integer> NodeClusters;
    }

    public class PhyTreeInfo {
        int numLeaves;
        int numBranches;
        int numNodes;
        List<List<Integer>> Pointers; // N x 2 matrix
        List<Double> Distances;
        public List<String> LeafNames;
        public List<String> BranchNames;
        public List<String> NodeNames;
    }

    public void printM(double[][] l) {
        int cols = l[0].length;
        for (double[] line : l) {
            for (int i = 0; i < cols; i++) {
                if (line.length>i)
                    System.out.printf("%10.5f",line[i]);
                else
                    System.out.printf("%10s","");
            }
            System.out.println();
        }
    }

    public void printM(boolean[][] l) {
        int cols = l[0].length;
        for (boolean[] line : l) {
            for (int i = 0; i < cols; i++) {
                if (line.length>i)
                    System.out.printf("%10.5B",line[i]);
                else
                    System.out.printf("%10s","");
            }
            System.out.println();
        }
    }

    public static void main(String[] args) {
        PhyTree t = new PhyTree(args[0]);
        t.cluster(Double.valueOf(args[1]));
    }
}

