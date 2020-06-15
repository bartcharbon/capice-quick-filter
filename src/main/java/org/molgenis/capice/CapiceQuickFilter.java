package org.molgenis.capice;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;

import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFSampleHeaderLine;
import java.io.File;
import java.util.*;

/**
 * CapiceQuickFilter tool. With only CAPICE and GnomAD annotations, perform
 * basic but effective filtering and reporting of potentially clinically
 * interesting variants.
 */
public class CapiceQuickFilter {

    /*
     * Class variables
     */
    private File input;
    private File output;
    private double capiceThreshold;
    private double gnomadThreshold;
    private String caseSampleID;
    private List<String> controlSampleIDs;

    /*
     * Static variables
     */
    private static final String version = "v0.0.1";
    private static final String DE_NOVO = "Potential de novo/uncontrolled hetzygote: ";
    private static final String HOM_ALT = "Potential homozygous                    : ";
    private static final String NON_AUT = "Potential non-autosomal                 : ";
    private static final String COMPHET = "Potential compound heterozygote         : ";

    /*
     * Constructor
     */
    public CapiceQuickFilter(File input, File output, double capiceThreshold, double gnomadThreshold, String caseSampleID, List<String> controlSampleIDs) {
        this.input = input;
        this.output = output;
        this.capiceThreshold = capiceThreshold;
        this.gnomadThreshold = gnomadThreshold;
        this.caseSampleID = caseSampleID;
        this.controlSampleIDs = controlSampleIDs;
    }

    /**
     * Run the CapiceQuickFilter after constructing.
     */
    void run() throws Exception
    {
        /*
         * Counters for reporting
         */
        int totalVariantCount = 0;
        int droppedByGnomAD = 0;
        int droppedByCAPICE = 0;
        int droppedByNullOrRefCaseGeno = 0;
        int droppedByHomZygAltControlGeno = 0;
        int droppedByHetZygAltNoHetComp = 0;
        int variantWithoutGnomAD = 0;
        int variantWithoutCAPICE = 0;

        /*
         * Initialize the VCF reader
         */
        VCFFileReader r = new VCFFileReader(input, false);
        /*
         * Prepare objects to store output
         */
        HashMap<String, List<VariantContext>> geneToHetZyg = new HashMap<>();
        HashMap<String, List<VariantContext>> reportedVariants = new HashMap<>();
        reportedVariants.put(DE_NOVO, new ArrayList<>());
        reportedVariants.put(HOM_ALT, new ArrayList<>());
        reportedVariants.put(NON_AUT, new ArrayList<>());
        reportedVariants.put(COMPHET, new ArrayList<>());

        /*
         * Get the sample names from the VCF meta-data
         * TODO: verify that order is guaranteed
         */
        List<String> sampleNames = new ArrayList<>();
        for(String sample: r.getFileHeader().getSampleNamesInOrder()){
            sampleNames.add(sample);
        }

        /*
         * Sanity checks: are the sample and control IDs present in the VCF?
         * Also, store indices of case and control samples for later use and
         * a list of all indices for convenience
         */
        if(!sampleNames.contains(caseSampleID))
        {
            throw new Exception("index sample id not found: " + caseSampleID);
        }

        /*
         * Start iterating over the input VCF file
         */
        List<String> trio = new ArrayList<>();
        trio.add(caseSampleID);
        trio.addAll(controlSampleIDs);
        for(VariantContext vc : r)
        {
            totalVariantCount++;

            /*
             * Retrieve GnomAD and CAPICE values if present
             */
            Double lowestGnomAD =
                    Helper.getLowestGnomAD(vc);
            Double highestCapice = Helper.getHighestCapice(vc);

            /*
             * Keep track of missing GnomAD and CAPICE values
             */
            if(highestCapice == null)
            {
                variantWithoutCAPICE++;
            }
            if(lowestGnomAD == null)
            {
                variantWithoutGnomAD++;
            }

            /*
             * If not missing, we have reasons to drop variants
             */
            if(highestCapice != null && highestCapice < capiceThreshold)
            {
                droppedByCAPICE++;
                continue;
            }
            if(lowestGnomAD != null && lowestGnomAD > gnomadThreshold)
            {
                droppedByGnomAD++;
                continue;
            }

            /*
             * Now we need to investigate genotypes
             * Store case and control genotype in objects
             */
            Genotype caseGeno = null;
            HashMap<String, Genotype> controlGeno = new HashMap<>();

            /*
             * Iterate over samples and extract the case and control genotypes
             */
            for(String sample : vc.getSampleNamesOrderedByName())
            {
                if(sample.equals(caseSampleID))
                {
                    caseGeno = vc.getGenotype(caseSampleID);
                }
                else if(controlSampleIDs.contains(sample))
                {
                    controlGeno.put(sample, vc.getGenotype(sample));
                }
            }

            /*
             * Drop variant if case genotype consists of only reference
             * alleles and/or missing alleles
             */
            int caseAltCount = 0;
            for(Allele a: caseGeno.getAlleles())
            {
                if(a.isNonRefAllele()) {
                    caseAltCount++;
                }
            }
            if(caseAltCount == 0)
            {
                droppedByNullOrRefCaseGeno++;
                continue;
            }

            /*
             * Analyse further. If one control sample is homozygous, drop the
             * variant. If not, track if 1+ control(s) are heterozygous.
             */
            boolean atLeastOneCtrlWithOneAlt = false;
            for(String key : controlGeno.keySet())
            {
                Genotype controlAlleles = controlGeno.get(key);
                int controlAltCount = 0;
                for(Allele a: controlAlleles.getAlleles())
                {
                    if(a.isNonRefAllele()) {
                        controlAltCount++;
                    }
                }
                if(controlAltCount == 2)
                {
                    droppedByHomZygAltControlGeno++;
                    continue;
                }
                else if(controlAltCount == 1)
                {
                    atLeastOneCtrlWithOneAlt = true;
                }
            }


            /*
             * There are no homozygous controls. So if the case is homozygous
             * alternative, report it and continue.
             */
            if(caseAltCount == 2)
            {
                reportedVariants.get((Helper.isAutosomal(vc) ? HOM_ALT : NON_AUT)).add(Helper.retainIndices(vc, trio));
                continue;
            }

            /*
             * If case has 1 alt, and there are no controls with alt alleles,
             * it is either de novo (controls present), or an 'uncontrolled'
             * heterozygote (controls not present). Report and continue.
             */
            if(caseAltCount == 1 && !atLeastOneCtrlWithOneAlt)
            {
                reportedVariants.get((Helper.isAutosomal(vc) ? DE_NOVO : NON_AUT)).add(Helper.retainIndices(vc, trio));
                continue;
            }

            /*
             * If case has 1 alt but there are also controls with 1 alt, it
             * could still be compound heterozygous. But we can only tell
             * after we have seen all variants from this gene. Save for later.
             * Exception is variants on allosomes, always report these.
             */
            Set<String> genes = Helper.getGenes(vc);
            if(caseAltCount == 1)
            {
                if(!Helper.isAutosomal(vc))
                {
                    reportedVariants.get(NON_AUT).add(Helper.retainIndices(vc, trio));
                    continue;
                }
                else
                {
                    for(String gene : genes)
                    {
                        if(!geneToHetZyg.containsKey(gene))
                        {
                            geneToHetZyg.put(gene, new ArrayList<>());
                        }
                        geneToHetZyg.get(gene).add(vc);
                    }
                      continue;
                }
            }

            /*
             * We should have covered all states when looping over all
             * variants in the input VCF. If not, crash the program.
             */
            throw new Exception("Bad state: all possibilities should be covered by now. Offending variant: " + vc.toStringDecodeGenotypes());
        }


        /*
         * Iterate over the heterozygous variants that may become compound
         * heterozygous variants if there are two in one gene. Keep track
         * which are reported to prevent duplicates.
         */
        Set<VariantContext> hasBeenReported = new HashSet<>();
        Set<VariantContext> hasBeenDropped = new HashSet<>();
        for(String gene : geneToHetZyg.keySet()) {
            if (geneToHetZyg.get(gene).size() > 1) {
                for (VariantContext variantContext : geneToHetZyg.get(gene)) {
                    if (!hasBeenReported.contains(variantContext)) {
                        reportedVariants.get(COMPHET).add(Helper.retainIndices(variantContext, trio));
                        hasBeenReported.add(variantContext);
                    }
                }
            }
        }

        /*
         * Also iterate over the leftovers to make sure all numbers add up.
         * For this, we must also considered those already reported for a
         * different gene.
         */
        for(String gene : geneToHetZyg.keySet()) {
            if (geneToHetZyg.get(gene).size() == 1)
            {
                VariantContext variantContext = geneToHetZyg.get(gene).get(0);
                if(!hasBeenDropped.contains(variantContext) && !hasBeenReported.contains(variantContext)) {
                    droppedByHetZygAltNoHetComp++;
                    hasBeenDropped.add(variantContext);
                }
            }
        }

        /*
         * Count total reported and total dropped
         */
        int totalRep =
                reportedVariants.get(HOM_ALT).size() + reportedVariants.get(DE_NOVO).size() + reportedVariants.get(COMPHET).size() + reportedVariants.get(NON_AUT).size();
        int totalDrop =
                droppedByGnomAD + droppedByCAPICE + droppedByNullOrRefCaseGeno + droppedByHomZygAltControlGeno + droppedByHetZygAltNoHetComp;

        /*
         * Print the header with information in the output VCF file.
         * TODO: retain original header for proper meta-data
         * TODO: assign VCF version
         * TODO: perhaps sort output by chrom/pos instead of category
         */
        VCFHeader header = r.getFileHeader();
        header.addMetaDataLine(new VCFHeaderLine("1","Output of CapiceQuickFilter " + version));
        header.addMetaDataLine(new VCFHeaderLine("2","Settings:"));
        header.addMetaDataLine(new VCFHeaderLine("3","- Input file: " + input.getAbsolutePath()));
        header.addMetaDataLine(new VCFHeaderLine("4","- Output file: " + output.getAbsolutePath()));
        header.addMetaDataLine(new VCFHeaderLine("5","- CAPICE threshold: " + capiceThreshold));
        header.addMetaDataLine(new VCFHeaderLine("6","- GnomAD threshold: " + gnomadThreshold));
        header.addMetaDataLine(new VCFHeaderLine("7","- Case sample ID: " + caseSampleID));
        header.addMetaDataLine(new VCFHeaderLine("8","- Control sample IDs: " + controlSampleIDs));
        header.addMetaDataLine(new VCFHeaderLine("9","Total number of variants processed: " + totalVariantCount));
        header.addMetaDataLine(new VCFHeaderLine("10","Total number of potential candidates found: " + totalRep));
        header.addMetaDataLine(new VCFHeaderLine("11","Breakdown of potential candidates by type:"));
        header.addMetaDataLine(new VCFHeaderLine("12","- " + HOM_ALT + reportedVariants.get(HOM_ALT).size()));
        header.addMetaDataLine(new VCFHeaderLine("13","- " + DE_NOVO + reportedVariants.get(DE_NOVO).size()));
        header.addMetaDataLine(new VCFHeaderLine("14","- " + COMPHET + reportedVariants.get(COMPHET).size()));
        header.addMetaDataLine(new VCFHeaderLine("15","- " + NON_AUT + reportedVariants.get(NON_AUT).size()));
        header.addMetaDataLine(new VCFHeaderLine("16","Total number of variants dropped: " + totalDrop));
        header.addMetaDataLine(new VCFHeaderLine("17","Breakdown of dropped variants by reason:"));
        header.addMetaDataLine(new VCFHeaderLine("18","- CAPICE score below threshold = " + droppedByCAPICE));
        header.addMetaDataLine(new VCFHeaderLine("19","- GnomAD allele frequency over threshold = " + droppedByGnomAD));
        header.addMetaDataLine(new VCFHeaderLine("20","- Case genotype null or reference = " + droppedByNullOrRefCaseGeno));
        header.addMetaDataLine(new VCFHeaderLine("21","- Homozygous control was present = " + droppedByHomZygAltControlGeno));
        header.addMetaDataLine(new VCFHeaderLine("22","- Flagged for compound but no second hit: " + droppedByHetZygAltNoHetComp));
        header.addMetaDataLine(new VCFHeaderLine("23","Additional information:"));
        header.addMetaDataLine(new VCFHeaderLine("24","- Variants without GnomAD annotation: " + variantWithoutGnomAD));
        header.addMetaDataLine(new VCFHeaderLine("25","- Variants without CAPICE annotation: " + variantWithoutCAPICE));
        header.addMetaDataLine(new VCFHeaderLine("26","Potential candidates categorized by type (full info below, can be copy-pasted side by side):"));
        for(String key : reportedVariants.keySet())
        {
            int i = 26;
            for(VariantContext variantContext : reportedVariants.get(key))
            {
                i++;
                String variant = variantContext.toString();
                header.addMetaDataLine(new VCFHeaderLine(""+i,key + (variant.length() > 70 ?
                        variant.substring(0, 70).replace("\t", " ") :
                        variant.replace("\t", " "))));
            }
        }

        Set<VCFHeaderLine> lines = new HashSet<>();
        for(VCFHeaderLine line : header.getMetaDataInInputOrder()){
            if(!(line instanceof VCFSampleHeaderLine)){
                lines.add(line);
            }
        }
        List<String> sortedSampleNames = new LinkedList<>();
        for(String sample: header.getSampleNamesInOrder()){
            if(trio.contains(sample)){
                sortedSampleNames.add(sample);
            }
        }
        VCFHeader vcfHeader = new VCFHeader(lines, sortedSampleNames);
        VariantContextWriter writer = createVCFWriter(output, vcfHeader);
        /*
         * Print the VCF columns with sample names and then all variant data.
         * By sorting the indices first, we maintain the same order used in
         * retainIndices() to print the genotypes.
         */
        Collections.sort(trio);
        for(String key : reportedVariants.keySet())
        {
            for(VariantContext variantContext : reportedVariants.get(key)) {
                writer.add(variantContext);
            }
        }
        writer.close();
    }

    static VariantContextWriter createVCFWriter(final File outFile, VCFHeader header) {
        VariantContextWriterBuilder vcWriterBuilder =
            new VariantContextWriterBuilder().clearOptions().setOutputFile(outFile);
        VariantContextWriter writer = vcWriterBuilder.build();

        writer.writeHeader(header);
        return writer;
    }
}
