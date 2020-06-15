package org.molgenis.capice;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Main class for running CapiceQuickFilter from command-line.
 */
public class Main
{
    public static void main(String args[]) throws Exception
    {
        /*
         * Print help if incorrect number of arguments are given
         */
        if(args.length != 5 && args.length != 6)
        {
            System.out.println("Please supply 5 or 6 arguments:");
            System.out.println("- File location of your input .VCF.GZ file.");
            System.out.println("- Output file location. May not exist yet.");
            System.out.println("- CAPICE score threshold. Lower scoring " +
                    "variants are dropped. Suggesting 0.2 for 90% sensitivity.");
            System.out.println("- GnomAD allele frequency threshold. Higher " +
                    "frequency variants are dropped. Suggesting 0.05 to be safe.");
            System.out.println("- Case sample ID (ie. proband,index).");
            System.out.println("- [optional] Control sample ID(s), " +
                    "comma-separated if multiple.");
            System.exit(0);
        }

        /*
         * Input .VCF.GZ file
         */
        File input = new File(args[0]);
        if(!(input.getName().endsWith(".vcf.gz") || input.getName().endsWith(".vcf")))
        {
            System.out.println("Input VCF file name '" + input.getName() + "' does not end in .vcf' or '.vcf.gz'. Are you sure this is a valid input?");
            System.exit(0);
        }
        if(!input.exists())
        {
            System.out.println("Input VCF file not found at " + input.getAbsolutePath()+".");
            System.exit(0);
        }

        /*
         * Output file
         */
        File output = new File(args[1]);
        if(output.exists())
        {
            System.out.println("Output VCF file already exists at " + output.getAbsolutePath()+". Please delete it first, or supply a different output file name.");
            System.exit(0);
        }

        /*
         * CAPICE threshold
         */
        String capiceThresholdStr = args[2];
        try {
            Double.parseDouble(capiceThresholdStr);
        } catch(NumberFormatException e){
            System.out.println("CAPICE threshold is not a decimal number: " + capiceThresholdStr);
            System.exit(0);
        }
        double capiceThreshold =  Double.parseDouble(capiceThresholdStr);
        if(capiceThreshold < 0.0 || capiceThreshold > 1.0)
        {
            System.out.println("CAPICE threshold must be between 0.0 and 1.0 " +
                    "instead of " + capiceThreshold);
            System.exit(0);
        }

        /*
         * GnomAD threshold
         */
        String gnomadThresholdStr = args[3];
        try {
            Double.parseDouble(gnomadThresholdStr);
        } catch(NumberFormatException e){
            System.out.println("GnomAD threshold is not a number: " + gnomadThresholdStr);
            System.exit(0);
        }
        double gnomadThreshold =  Double.parseDouble(gnomadThresholdStr);
        if(gnomadThreshold < 0.0 || gnomadThreshold > 1.0)
        {
            System.out.println("GnomAD threshold must be between 0.0 and 1.0 " +
                    "instead of " + gnomadThreshold);
            System.exit(0);
        }

        /*
         * Case sample ID within input VCF
         */
        String caseSampleID = args[4];
        if(caseSampleID.isEmpty())
        {
            System.out.println("Case sample ID may not be empty.");
            System.exit(0);
        }

        /*
         * Control sample ID(s) within input VCF
         */
        List<String> controlSampleIDs;
        if(args.length == 6)
        {
            String controlSampleIDsStr = args[5];
            String[] controlSampleIDsArr = controlSampleIDsStr.split(",", -1);
            for(String controlID : controlSampleIDsArr)
            {
                if(controlID.isEmpty())
                {
                    System.out.println("Control sample ID may not be empty.");
                    System.exit(0);
                }
            }
            controlSampleIDs = Arrays.asList(controlSampleIDsArr);
        }
        else
        {
            controlSampleIDs = new ArrayList<>();
        }

        /*
         * Argument parsing done, start CapiceQuickFilter
         */
        System.out.println("Arguments OK. Starting...");
        long start = System.nanoTime();
        CapiceQuickFilter cqf = new CapiceQuickFilter(input, output,
                capiceThreshold, gnomadThreshold, caseSampleID, controlSampleIDs);
        cqf.run();
        System.out.println("...completed in " + ((System.nanoTime()-start)/1000000)+"ms.");
    }
}
