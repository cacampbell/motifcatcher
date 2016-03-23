package motifcatcher;

public class SiteEntry {
    private String name;
    private String fullName;

    TreeWithClusters.SeedEntry seedEntry;

    private Integer start;
    private Integer end;

    public SiteEntry(String name, TreeWithClusters.SeedEntry se) {
        this.name = name;
        this.seedEntry = se;
    }

    public Integer getStart() {
        return start;
    }

    public void setStart(Integer start) {
        this.start = start;
        this.end = start + seedEntry.getMotifLength();
    }

    public String getName() {
        return name;
    }

    public String getFullName() {
        return fullName;
    }

    public void setFullName(String fullName) {
        this.fullName = fullName;
    }
}
