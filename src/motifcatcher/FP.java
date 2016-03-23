package motifcatcher;

/**
* Created with IntelliJ IDEA.
* User: cobalt
* Date: 12.06.2012
* Time: 14:36
* To change this template use File | Settings | File Templates.
*/
public class FP {
    String name;
    String dir;
    String file;
    Float matrix[][];

    public FP(String s) {
        this.name = s;
    }

    public FP(String s, String dir) {
        this.name = s;
        this.dir = dir;
    }

    public String getDir() {
        return dir;
    }

    public void setDir(String dir) {
        this.dir = dir;
    }

    public String getFile() {
        return file;
    }

    public void setFile(String file) {
        this.file = file;
    }

    public Float[][] getMatrix() {
        return matrix;
    }

    public void setMatrix(Float[][] matrix) {
        this.matrix = matrix;
    }
}
