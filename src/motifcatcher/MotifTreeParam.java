package motifcatcher;

//TODO: Update these parameters
//TODO: Add options for different underlying motif search programs
//TODO: Generify motif search program, make / extend API for different motif finders
public class MotifTreeParam {
    public String ColComparison = "ALLR"; // PCC, ALLR, ALLRLL, CS, KL, SSD
    public String AlignmentScheme = "SWA"; // NW, SW, SWA, SWU
    public Double EValThreshold = 0.001;
    public int MaxR = 100;
    public Double gapOpenPen = 1.00;
    public Double gapExtPen = 0.50;
}
