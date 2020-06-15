package org.molgenis.capice;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import htsjdk.variant.variantcontext.VariantContextBuilder;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Helper class with some self-contained functions used in CapiceQuickFilter.
 */
public class Helper {

  /**
   * Retrieve genes symbols are annoted by VEP.
   */
  static Set<String> getGenes(VariantContext vc) {
    Set<String> genes = new HashSet<>();

    List<String> val = getInfoValAsList(vc, "CSQ");
    for (String forEachAltAndTranscript : val) {
      String[] csqField =
          forEachAltAndTranscript.split("\\|", -1);
      String gene = csqField[3];
      if (!gene.isEmpty()) {
        genes.add(gene);
      }
    }
    return genes;
  }

  /**
   * Get highest CAPICE score, or NULL if CAPICE is not present. Note that we are not matching exact
   * alt allele here. If one variant has a high enough score, it passes for further interpretation.
   */
  static Double getHighestCapice(VariantContext vc) {
    Double highestCapice = null;
    List<String> val = getInfoValAsList(vc, "CAPICE");
    for (String CS : val) {
      if(CS == null){
        System.err.println("null capice score for line: " + vc.toStringDecodeGenotypes());
        return Double.valueOf(-1);
      }
      double csDouble = Double.parseDouble(CS);
      if (highestCapice == null || csDouble > highestCapice) {
        highestCapice = csDouble;
      }
    }

    return highestCapice;
  }

  private static List<String> getInfoValAsList(VariantContext vc, String capice) {
    Object raw = vc.getAttribute(capice);
    if(raw instanceof List){
      return (List<String>) raw;
    }
    return Collections.singletonList((String) raw);
  }

  /**
   * Get lowest GnomAD allele frequency, or NULL if GnomAD is not present. Note that we are not
   * matching exact alt allele here. If one variant is rare enough, it passes for further
   * interpretation.
   */
  static Double getLowestGnomAD(VariantContext vc) {
    Double lowestGnomAD = null;
    List<String> val = getInfoValAsList(vc, "CSQ");
    for (String forEachAltAndTranscript : val) {
      String[] csqField =
          forEachAltAndTranscript.split("\\|", -1);
      String gnomadAFStr = csqField[26];
      if (!gnomadAFStr.isEmpty()) {
        double gnomadAF = Double.parseDouble(gnomadAFStr);
        if (lowestGnomAD == null || gnomadAF < lowestGnomAD) {
          lowestGnomAD = gnomadAF;
        }
      }
    }

    return lowestGnomAD;
  }

  /**
   * Return a VCF line but only genotypes for selected indices.
   */
  static VariantContext retainIndices(VariantContext variant,
      List<String> trio) {
      VariantContextBuilder variantContextBuilder = new VariantContextBuilder(variant);
      List<Genotype> genotypes = new ArrayList<>();
      for(String sample : trio){
          genotypes.add(variant.getGenotype(sample));
      }
      variantContextBuilder.genotypes(genotypes);
    return variantContextBuilder.make();
  }

  /**
   * Check if a chromosome is autosomal or not.
   */
  static boolean isAutosomal(VariantContext vc) {
    return vc.getContig().matches("\\d+(\\.\\d+)?");
  }
}
