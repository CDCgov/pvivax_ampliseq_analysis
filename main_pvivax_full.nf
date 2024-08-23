nextflow.enable.dsl=2

include { Trim_reads } from './modules/mars_modules/trim_read'
include { PreFastqC; pre_multiqc; PostFastqC; multiqc} from './modules/mars_modules/fastqc'
include {BWA_index} from './modules/mars_modules/index'
include {BWA_align} from './modules/mars_modules/align'
include {Sam_sort; Picard_add_read; VCF_call; Get_Bed} from './modules/mars_modules/vcf_call'
include {buildsnpeff_db; annotation; vartype  } from './modules/mars_modules/annotation'
include {vcf_to_DF; csv_merge} from './modules/mars_modules/csv_merge'
include {getcoverage; WT_cov ; Trim_Stats; Reads_merge} from './modules/mars_modules/coverage'
include {Snpfilter; Summary_merge; Summary; Dataviz_Reportable_snps; DataViz_Novel_snps; Introns_merge } from './modules/mars_modules/final_snp'


include { trainBALK; mapPhix; fastqc_beforeTrim; trimReads; humanMap; noHumanReads; indexRefs; mapToRefs; mapToRefs_bowtie2; cleanSAM_bowtie2; gatkHapCaller_bowtie2; indexBAM_bowtie2; cleanSAM; markDupes; gatkHapCaller; bcftoolsMpileup; bcftoolsCall; bcftoolsFilter_Index; gatherGVCFs; genotypeGATK; variantsToTable_gatk; gatkBarcode; predictGeo_gatk;  mergeBCFtools; variantsToTable_bcftools; bcftoolsBarcode; predictGeo_bcftools; mergePredictions; indexBAM } from './modules/geoPrediction_modules/geoAmpliseq_qc_reads'

include {bbduk; fastqc_beforeMerge; fastqc_afterMerge } from './modules/treeClustering_modules/workflow1_readInFastqs'

include { makeRefs; makeRefs_v2; finalNewRefs; finalNewRefs_v2; mapping; mapping2fastq; bowtie2; bedExtraction; bedExtraction_plusNewHaps; makeInitialReads; makeFinalReads; makeFinalReads_v2; newHaps1; nameNewHaps; makeRefs_assignHaps2; makeRefs_assignHaps2_v2; finalNewRefs_assignHaps2; blastDB; firstClustering; blastMulti; hapsBlast; firstClustering_v2; blastMulti_v2; newHaps1_v2; bedExtraction_v2; nameNewHaps_sameLength } from './modules/treeClustering_modules/workflow1_SNPmarkers'

include { combineGenotypes; copyGenotypes; haplotypeSheet; haplotypeSheet_specific } from './modules/treeClustering_modules/workflow1_writeHaplotypes'

include { individualMats; individualChroms; scaleMatrices } from './modules/treeClustering_modules/workflow2'
include { runModule3_Pvivax; runModule3_parnas } from './modules/treeClustering_modules/workflow3'

include { makeTree_full; updateGeoPredictionDB; makeTree_states } from './modules/reporting_modules/pvivax_tree'
include { marsVOI_allStates; marsVOI_singleState; runReportGen_allStates; runReportGen_singleState } from './modules/reporting_modules/pvivax_reports'


tmpDir = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2")
tmpDir.deleteDir()
tmpDir.mkdir()

processReads = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/PROCESSED_READS")
if (! processReads.exists()){
  processReads.mkdir()
}

qcFolder = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/QClogs")



renamedReads = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/RENAMED_READS")
renamedReads.mkdir()

refHapsDir = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/REF_HAPS")
refHapsDir.mkdir()
tmpGenotypes = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/TMP_GENOTYPES")
tmpGenotypes.mkdir()

nonJunc_mapping = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/NON_JUNC_MAPPING")
nonJunc_mapping.mkdir()
if (! nonJunc_mapping.exists()){
  nonJunc_mapping.mkdir()
}

nonJunc_NEWHAPS = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/NON_JUNC_MAPPING/nonJunc_NEW_HAPS")
nonJunc_NEWHAPS.mkdir()
nonJunc_clippedReads = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/NON_JUNC_MAPPING/CLIPPED_READS")
nonJunc_clippedReads.mkdir()

NEW_HAPS_DIR = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/REF_SEQS/BLASTING/NEW_HAPS")
if (! NEW_HAPS_DIR.exists()){
  NEW_HAPS_DIR.mkdir()
}

// Junc_mapping = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/JUNC_MAPPING")
// Junc_mapping.mkdir()
// if (! Junc_mapping.exists()){
//   Junc_mapping.mkdir()
// }

// NEW_HAPS_DIR_Junc = file("$baseDir/HAPLOTYPE_CALLER_PVIVAX/TMP2/JUNC_New_haps")
// if (! NEW_HAPS_DIR_Junc.exists()){
//   NEW_HAPS_DIR_Junc.mkdir()
// }


workflow P_vivax_MaRS_GeoPrediction {
   // cerate a channel for ref and reads
    myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
    ref_ch =  Channel.fromPath(params.ref, checkIfExists: true)
    // read_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).filter{ it.size()>0 }
    adapter_ch = Channel.fromPath(params.adapters, checkIfExists: true)
    snpeff_ch = Channel.fromPath(params.snpeff_config, checkIfExists: true)
    gff_ch = Channel.fromPath( params.gff, checkIfExists: true )
    voi_ch = Channel.fromPath( params.voi, checkIfExists: true )
  //  bed_ch = Channel.fromPath( params.bed, checkIfExists: true )

    PreFastqC(myReads)
    Trim_reads(adapter_ch.combine(myReads))
    humanMap(Trim_reads.out.Trimmed_fastq)
    noHumanReads(humanMap.out.sortedNoHumanBAM)
    PostFastqC(noHumanReads.out.humanRemoveFASTQ)
    multiqc(PostFastqC.out.collect())

  

    buildsnpeff_db()
    //PreFastqC(read_ch)
    // Trim_reads(adapter_ch.combine(myReads))
    // PostFastqC(Trim_reads.out.Trimmed_fastq)
    // multiqc(PostFastqC.out.collect())


    BWA_index(ref_ch)
    BWA_align(BWA_index.out.collect().combine(noHumanReads.out.humanRemoveFASTQ))
    Sam_sort(BWA_align.out.Bwa_samfile)
    Picard_add_read(BWA_align.out.Bwa_samfile)
    Get_Bed (gff_ch)
    VCF_call(BWA_index.out.combine(Get_Bed.out).combine(Picard_add_read.out.Picard_out_bam))
    annotation(snpeff_ch.combine(VCF_call.out.variants), buildsnpeff_db.out.buildDB)
    vartype(annotation.out.var_annotation)
    getcoverage(Picard_add_read.out.Picard_out_bam)
    Trim_Stats(Trim_reads.out.Trimmed_stats.join(getcoverage.out.Samtools_coverage))

    WT_cov(ref_ch.combine(gff_ch).combine(voi_ch).combine(getcoverage.out.samtools_depth))
    Reads_merge(Trim_Stats.out.Reads_cov.collect())
    vcf_to_DF(vartype.out.vartype_annotation)
    csv_merge(vcf_to_DF.out.csv_annotate)
    Snpfilter(voi_ch.combine(csv_merge.out.CSV_merge.join(WT_cov.out.WT_coverage)))
    Summary_merge(Snpfilter.out.snp_report.collect())
    Summary(Summary_merge.out)
    Introns_merge(csv_merge.out.Introns_file.collect())
    Dataviz_Reportable_snps(Summary.out.Reportable_snps)
    DataViz_Novel_snps(Summary.out.Novel_snps)
    
    mapToRefs_bowtie2(noHumanReads.out.humanRemoveFASTQ)
    cleanSAM_bowtie2(mapToRefs_bowtie2.out.bowtieMapped)
    // markDupes(cleanSAM.out.cleanedSAM)
    indexBAM_bowtie2(cleanSAM_bowtie2.out.cleanedSAM)
    gatkHapCaller_bowtie2(indexBAM_bowtie2.out.indexOut)

    gatk_allFiles = gatkHapCaller_bowtie2.out.gatkOut.collect()
    gatherGVCFs(gatk_allFiles.collect())
    genotypeGATK(gatherGVCFs.out.combinedGVCFs)
    variantsToTable_gatk(genotypeGATK.out.gatkVCF)
    gatkBarcode(variantsToTable_gatk.out.gatkTable, gatk_allFiles.collect())
    predictGeo_gatk(gatkBarcode.out.barcodeOut)


    bcftoolsMpileup(indexBAM_bowtie2.out.indexOut)
    bcftoolsCall(bcftoolsMpileup.out.pileupOut)
    bcftoolsFilter_Index(bcftoolsCall.out.callOut)
    bcftools_allFiles = bcftoolsFilter_Index.out.vcfOnly.collect()
    mergeBCFtools(bcftools_allFiles, bcftoolsFilter_Index.out.vcfIndex.collect())
    variantsToTable_bcftools(mergeBCFtools.out.bcftoolsMerged)
    bcftoolsBarcode(variantsToTable_bcftools.out.bcftoolsTable, bcftools_allFiles)
    predictGeo_bcftools(bcftoolsBarcode.out.barcodeOut)
    mergePredictions(predictGeo_gatk.out.gatkCountry, predictGeo_gatk.out.gatkRegion, predictGeo_bcftools.out.bcfCountry, predictGeo_bcftools.out.bcfRegion, gatkBarcode.out.barcodeOut, bcftoolsBarcode.out.barcodeOut)

    marsVOI_allStates(mergePredictions.out.regionPredictionMerge, Summary.out.Summary_Snp)

  

}

workflow P_vivax_MaRS {

    // cerate a channel for ref and reads
     myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
    ref_ch =  Channel.fromPath(params.ref, checkIfExists: true)
    // read_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).filter{ it.size()>0 }
    adapter_ch = Channel.fromPath(params.adapters, checkIfExists: true)
    snpeff_ch = Channel.fromPath(params.snpeff_config, checkIfExists: true)
    gff_ch = Channel.fromPath( params.gff, checkIfExists: true )
    voi_ch = Channel.fromPath( params.voi, checkIfExists: true )
  //  bed_ch = Channel.fromPath( params.bed, checkIfExists: true )

    buildsnpeff_db()
    //PreFastqC(read_ch)
    Trim_reads(adapter_ch.combine(myReads))
    PostFastqC(Trim_reads.out.Trimmed_fastq)
    multiqc(PostFastqC.out.collect())
    BWA_index(ref_ch)
    BWA_align(BWA_index.out.collect().combine(Trim_reads.out.Trimmed_fastq))
    Sam_sort(BWA_align.out.Bwa_samfile)
    Picard_add_read(BWA_align.out.Bwa_samfile)
    Get_Bed (gff_ch)
    VCF_call(BWA_index.out.combine(Get_Bed.out).combine(Picard_add_read.out.Picard_out_bam))
    annotation(snpeff_ch.combine(VCF_call.out.variants), buildsnpeff_db.out.buildDB)
    vartype(annotation.out.var_annotation)
    getcoverage(Picard_add_read.out.Picard_out_bam)
    Trim_Stats(Trim_reads.out.Trimmed_stats.join(getcoverage.out.Samtools_coverage))

    WT_cov(ref_ch.combine(gff_ch).combine(voi_ch).combine(getcoverage.out.samtools_depth))
    Reads_merge(Trim_Stats.out.Reads_cov.collect())
    vcf_to_DF(vartype.out.vartype_annotation)
    csv_merge(vcf_to_DF.out.csv_annotate)
    Snpfilter(voi_ch.combine(csv_merge.out.CSV_merge.join(WT_cov.out.WT_coverage)))
    Summary_merge(Snpfilter.out.snp_report.collect())
    Summary(Summary_merge.out)
    Introns_merge(csv_merge.out.Introns_file.collect())
    Dataviz_Reportable_snps(Summary.out.Reportable_snps)
    DataViz_Novel_snps(Summary.out.Novel_snps)



  }


workflow P_vivax_GeoClassifier_complete_pipeline {


  myReads =Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
  

  trimReads(myReads)

  humanMap(trimReads.out.trimmedFASTQ)

  noHumanReads(humanMap.out.sortedNoHumanBAM)

  mapToRefs_bowtie2(noHumanReads.out.humanRemoveFASTQ)
  cleanSAM_bowtie2(mapToRefs_bowtie2.out.bowtieMapped)
  // markDupes(cleanSAM.out.cleanedSAM)
  indexBAM_bowtie2(cleanSAM_bowtie2.out.cleanedSAM)
  gatkHapCaller_bowtie2(indexBAM_bowtie2.out.indexOut)

  gatk_allFiles = gatkHapCaller_bowtie2.out.gatkOut.collect()
  gatherGVCFs(gatk_allFiles.collect())
  genotypeGATK(gatherGVCFs.out.combinedGVCFs)
  variantsToTable_gatk(genotypeGATK.out.gatkVCF)
  gatkBarcode(variantsToTable_gatk.out.gatkTable, gatk_allFiles.collect())
  predictGeo_gatk(gatkBarcode.out.barcodeOut)


  bcftoolsMpileup(indexBAM_bowtie2.out.indexOut)
  bcftoolsCall(bcftoolsMpileup.out.pileupOut)
  bcftoolsFilter_Index(bcftoolsCall.out.callOut)
  bcftools_allFiles = bcftoolsFilter_Index.out.vcfOnly.collect()
  mergeBCFtools(bcftools_allFiles, bcftoolsFilter_Index.out.vcfIndex.collect())
  variantsToTable_bcftools(mergeBCFtools.out.bcftoolsMerged)
  bcftoolsBarcode(variantsToTable_bcftools.out.bcftoolsTable, bcftools_allFiles)
  predictGeo_bcftools(bcftoolsBarcode.out.barcodeOut)
  mergePredictions(predictGeo_gatk.out.gatkCountry, predictGeo_gatk.out.gatkRegion, predictGeo_bcftools.out.bcfCountry, predictGeo_bcftools.out.bcfRegion, gatkBarcode.out.barcodeOut, bcftoolsBarcode.out.barcodeOut)

  

}


workflow P_vivax_geoCombined_BALK_useExistingVCF {

/*
  myReads =Channel.fromFilePairs("${params.fastqDir}/*_{1,2}.fastq.gz")
 
 //refs are already indexed, but can include that as part of script. Maybe check if already indexed or something, probably best

  trimReads(myReads)

  humanMap(trimReads.out.trimmedFASTQ)

  noHumanReads(humanMap.out.sortedNoHumanBAM)

  

  mapToRefs(noHumanReads.out.humanRemoveFASTQ)

  cleanSAM(mapToRefs.out.mappedSAM)

  markDupes(cleanSAM.out.cleanedSAM)

  indexBAM(markDupes.out.singleFile)

  gatkHapCaller(indexBAM.out.indexOut)

  indexBAM.out.indexOut.view()
 bcftoolsMpileup(indexBAM.out.indexOut)
 bcftoolsCall(bcftoolsMpileup.out.pileupOut)

  gatk_allFiles = gatkHapCaller.out.gatkOut.collect()
*/

  gatk_allFiles = Channel.fromPath("${params.variantFolder}/*ampliseq_gatk_bowtie2.g.vcf")

  gatherGVCFs(gatk_allFiles.collect())

  genotypeGATK(gatherGVCFs.out.combinedGVCFs)
  variantsToTable_gatk(genotypeGATK.out.gatkVCF)
  gatkBarcode(variantsToTable_gatk.out.gatkTable, gatk_allFiles.collect())

  predictGeo_gatk(gatkBarcode.out.barcodeOut)

  bcfFiles =  Channel.fromPath("${params.variantFolder}/*bcf.vcf")

  bcftoolsFilter_Index(bcfFiles)

  bcftools_allFiles = bcftoolsFilter_Index.out.vcfOnly.collect()
 
  mergeBCFtools(bcftools_allFiles, bcftoolsFilter_Index.out.vcfIndex.collect())

  variantsToTable_bcftools(mergeBCFtools.out.bcftoolsMerged)

  bcftoolsBarcode(variantsToTable_bcftools.out.bcftoolsTable, bcftools_allFiles)

  predictGeo_bcftools(bcftoolsBarcode.out.barcodeOut)

  mergePredictions(predictGeo_gatk.out.gatkCountry, predictGeo_gatk.out.gatkRegion, predictGeo_bcftools.out.bcfCountry, predictGeo_bcftools.out.bcfRegion, gatkBarcode.out.barcodeOut, bcftoolsBarcode.out.barcodeOut)

}


//Haplotype calling only - will create a haplotype sheet at the end
workflow Workflow1_v2_vivax{
  myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_*.fastq.gz")
  newHaps = Channel.fromPath("${params.origNewHaps}/*.fasta")

  Channel
    .fromPath("${params.loci}")
    .splitText()
    .set{ loci_ch }

  bbduk(myReads)
  // fastqc_beforeMerge(myReads)

  // fastqc_afterMerge(bbduk.out.readyReads)

  makeRefs_v2()
  finalNewRefs_v2(makeRefs_v2.out.completeFile.collectFile())
  mapping(finalNewRefs_v2.out.outFinalRefs,bbduk.out.readyReads)

  mapping2fastq(mapping.out.mappedOnlySam)
  bowtie2(mapping2fastq.out.mappedFastq)

  bedExtraction_plusNewHaps(loci_ch, bowtie2.out.finishSorting, finalNewRefs_v2.out.outShortRefs)
  makeFinalReads_v2(bedExtraction_plusNewHaps.out.tryBedFastq.groupTuple())
  
  nameNewHaps(bedExtraction_plusNewHaps.out.allNewHaps.collectFile(), finalNewRefs_v2.out.outShortRefs)

  makeRefs_assignHaps2_v2()

  finalNewRefs_assignHaps2(makeRefs_v2.out.completeFile.collectFile(), nameNewHaps.out.newHaps_singleFile.collectFile().ifEmpty("$baseDir/modules/files_forNextflow/empty.fasta"))
  blastDB(finalNewRefs_assignHaps2.out.outFinalRefs)
  
  firstClustering_v2(makeFinalReads_v2.out.finalClipped)
  blastMulti_v2(firstClustering_v2.out.multiTuple, firstClustering_v2.out.depthTuple, blastDB.out.outFinalRefs_db)
  hapsBlast(blastMulti_v2.out.blastFasta, blastDB.out.outFinalRefs_db)

  // place the genotypes in the correct folder with the correct name
  copyGenotypes(hapsBlast.out.genotypeNoJunction.ifEmpty("$baseDir/modules/files_forNextflow/emptySNP_genotype").collectFile())

  // create a haplotype sheet and copy the newly identified haplotypes to the correct folder
  haplotypeSheet(copyGenotypes.out.endFile.collectFile())

}

//Use specimen genotypes folders to create a haplotype sheet
workflow hapSheetOnly{

  haplotypeSheet_specific()
}


//make an ensemeble matrix, determine parnas clusters, and create a tree for all specimens in the dataset
workflow Workflow2and3_Pvivax{
  chromosomes = Channel.fromList(['Nu_CHR01','Nu_CHR02','Nu_CHR03','Nu_CHR04','Nu_CHR05','Nu_CHR06','Nu_CHR07','Nu_CHR08','Nu_CHR09','Nu_CHR10','Nu_CHR11','Nu_CHR12','Nu_CHR13','Nu_CHR14'])
  updatedDB = Channel.fromPath("$baseDir/treeFolder/*geoPrediction_simple.txt").toSortedList().flatten().last()
  hapSheets = Channel.fromPath("$baseDir/haplotype_sheets/*haplotype_sheet_largeFiles.txt").toSortedList().flatten().last()
  paramInput_wf2 = Channel.fromPath(params.wf2HapSheet)

  if( params.wf2HapSheet == "$baseDir/modules/files_forNextflow/empty.fasta" ){
    individualChroms(hapSheets, chromosomes)
    individualMats(individualChroms.out.chromHapSheet)
    scaleMatrices(individualMats.out.matrixOut.collect(), individualChroms.out.chromHapSheet.collect())
    runModule3_Pvivax(scaleMatrices.out.finalMatrix)
    runModule3_parnas(runModule3_Pvivax.out.tree, runModule3_Pvivax.out.threshold)
    updateGeoPredictionDB(updatedDB)
    makeTree_full(scaleMatrices.out.finalMatrix, runModule3_parnas.out.clustersOutput, updateGeoPredictionDB.out.geoPredictionUpdate,  updateGeoPredictionDB.out.travelUpdate)
    
  }
  else{
    individualChroms(paramInput_wf2, chromosomes)
    individualMats(individualChroms.out.chromHapSheet)
    scaleMatrices(individualMats.out.matrixOut.collect(), individualChroms.out.chromHapSheet.collect())
    runModule3_Pvivax(scaleMatrices.out.finalMatrix)
    runModule3_parnas(runModule3_Pvivax.out.tree, runModule3_Pvivax.out.threshold)
    updateGeoPredictionDB(updatedDB)
    makeTree_full(scaleMatrices.out.finalMatrix, runModule3_parnas.out.clustersOutput, updateGeoPredictionDB.out.geoPredictionUpdate, updateGeoPredictionDB.out.travelUpdate)
    
  }


}


workflow Workflow1through3_vivax {
  //workflow 1
  myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_*.fastq.gz")
  newHaps = Channel.fromPath("${params.origNewHaps}/*.fasta")

  Channel
    .fromPath("${params.loci}")
    .splitText()
    .set{ loci_ch }

  bbduk(myReads)
  // fastqc_beforeMerge(myReads)

  // fastqc_afterMerge(bbduk.out.readyReads)

  makeRefs_v2()
  finalNewRefs_v2(makeRefs_v2.out.completeFile.collectFile())
  mapping(finalNewRefs_v2.out.outFinalRefs,bbduk.out.readyReads)

  mapping2fastq(mapping.out.mappedOnlySam)
  bowtie2(mapping2fastq.out.mappedFastq)

  bedExtraction_plusNewHaps(loci_ch, bowtie2.out.finishSorting, finalNewRefs_v2.out.outShortRefs)
  makeFinalReads_v2(bedExtraction_plusNewHaps.out.tryBedFastq.groupTuple())
  
  nameNewHaps(bedExtraction_plusNewHaps.out.allNewHaps.collectFile(), finalNewRefs_v2.out.outShortRefs)

  makeRefs_assignHaps2_v2()

  finalNewRefs_assignHaps2(makeRefs_v2.out.completeFile.collectFile(), nameNewHaps.out.newHaps_singleFile.collectFile().ifEmpty("$baseDir/modules/files_forNextflow/empty.fasta"))
  blastDB(finalNewRefs_assignHaps2.out.outFinalRefs)
  
  firstClustering_v2(makeFinalReads_v2.out.finalClipped)
  blastMulti_v2(firstClustering_v2.out.multiTuple, firstClustering_v2.out.depthTuple, blastDB.out.outFinalRefs_db)
  hapsBlast(blastMulti_v2.out.blastFasta, blastDB.out.outFinalRefs_db)

  // place the genotypes in the correct folder with the correct name
  copyGenotypes(hapsBlast.out.genotypeNoJunction.ifEmpty("$baseDir/modules/files_forNextflow/emptySNP_genotype").collectFile())

  // create a haplotype sheet and copy the newly identified haplotypes to the correct folder
  haplotypeSheet(copyGenotypes.out.endFile.collectFile())

  //workflow 2 and 3
  chromosomes = Channel.fromList(['Nu_CHR01','Nu_CHR02','Nu_CHR03','Nu_CHR04','Nu_CHR05','Nu_CHR06','Nu_CHR07','Nu_CHR08','Nu_CHR09','Nu_CHR10','Nu_CHR11','Nu_CHR12','Nu_CHR13','Nu_CHR14'])
  updatedDB = Channel.fromPath("$baseDir/treeFolder/*geoPrediction_simple.txt").toSortedList().flatten().last()
  //hapSheets = Channel.fromPath("$baseDir/haplotype_sheets/*.txt").toSortedList().flatten().last()
  //paramInput_wf2 = Channel.fromPath(params.wf2HapSheet)

  individualChroms(haplotypeSheet.out.finalHapSheetLarge, chromosomes)
  individualMats(individualChroms.out.chromHapSheet)
  scaleMatrices(individualMats.out.matrixOut.collect(), individualChroms.out.chromHapSheet.collect())
  runModule3_Pvivax(scaleMatrices.out.finalMatrix)
  runModule3_parnas(runModule3_Pvivax.out.tree, runModule3_Pvivax.out.threshold)
  updateGeoPredictionDB(updatedDB)
  makeTree_full(scaleMatrices.out.finalMatrix, runModule3_parnas.out.clustersOutput, updateGeoPredictionDB.out.geoPredictionUpdate,  updateGeoPredictionDB.out.travelUpdate)
  
}

workflow everything {
  //MaRS and geoPrediction
  // cerate a channel for ref and reads
    myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_001.fastq.gz")
    ref_ch =  Channel.fromPath(params.ref, checkIfExists: true)
    // read_ch = Channel.fromFilePairs(params.reads, checkIfExists: true).filter{ it.size()>0 }
    adapter_ch = Channel.fromPath(params.adapters, checkIfExists: true)
    snpeff_ch = Channel.fromPath(params.snpeff_config, checkIfExists: true)
    gff_ch = Channel.fromPath( params.gff, checkIfExists: true )
    voi_ch = Channel.fromPath( params.voi, checkIfExists: true )
  //  bed_ch = Channel.fromPath( params.bed, checkIfExists: true )

    PreFastqC(myReads)
    Trim_reads(adapter_ch.combine(myReads))
    humanMap(Trim_reads.out.Trimmed_fastq)
    noHumanReads(humanMap.out.sortedNoHumanBAM)
    PostFastqC(noHumanReads.out.humanRemoveFASTQ)
    multiqc(PostFastqC.out.collect())

  

    buildsnpeff_db()
    //PreFastqC(read_ch)
    // Trim_reads(adapter_ch.combine(myReads))
    // PostFastqC(Trim_reads.out.Trimmed_fastq)
    // multiqc(PostFastqC.out.collect())


    BWA_index(ref_ch)
    BWA_align(BWA_index.out.collect().combine(noHumanReads.out.humanRemoveFASTQ))
    Sam_sort(BWA_align.out.Bwa_samfile)
    Picard_add_read(BWA_align.out.Bwa_samfile)
    Get_Bed (gff_ch)
    VCF_call(BWA_index.out.combine(Get_Bed.out).combine(Picard_add_read.out.Picard_out_bam))
    annotation(snpeff_ch.combine(VCF_call.out.variants), buildsnpeff_db.out.buildDB)
    vartype(annotation.out.var_annotation)
    getcoverage(Picard_add_read.out.Picard_out_bam)
    Trim_Stats(Trim_reads.out.Trimmed_stats.join(getcoverage.out.Samtools_coverage))

    WT_cov(ref_ch.combine(gff_ch).combine(voi_ch).combine(getcoverage.out.samtools_depth))
    Reads_merge(Trim_Stats.out.Reads_cov.collect())
    vcf_to_DF(vartype.out.vartype_annotation)
    csv_merge(vcf_to_DF.out.csv_annotate)
    Snpfilter(voi_ch.combine(csv_merge.out.CSV_merge.join(WT_cov.out.WT_coverage)))
    Summary_merge(Snpfilter.out.snp_report.collect())
    Summary(Summary_merge.out)
    Introns_merge(csv_merge.out.Introns_file.collect())
    Dataviz_Reportable_snps(Summary.out.Reportable_snps)
    DataViz_Novel_snps(Summary.out.Novel_snps)
    
    mapToRefs_bowtie2(noHumanReads.out.humanRemoveFASTQ)
    cleanSAM_bowtie2(mapToRefs_bowtie2.out.bowtieMapped)
    // markDupes(cleanSAM.out.cleanedSAM)
    indexBAM_bowtie2(cleanSAM_bowtie2.out.cleanedSAM)
    gatkHapCaller_bowtie2(indexBAM_bowtie2.out.indexOut)

    gatk_allFiles = gatkHapCaller_bowtie2.out.gatkOut.collect()
    gatherGVCFs(gatk_allFiles.collect())
    genotypeGATK(gatherGVCFs.out.combinedGVCFs)
    variantsToTable_gatk(genotypeGATK.out.gatkVCF)
    gatkBarcode(variantsToTable_gatk.out.gatkTable, gatk_allFiles.collect())
    predictGeo_gatk(gatkBarcode.out.barcodeOut)


    bcftoolsMpileup(indexBAM_bowtie2.out.indexOut)
    bcftoolsCall(bcftoolsMpileup.out.pileupOut)
    bcftoolsFilter_Index(bcftoolsCall.out.callOut)
    bcftools_allFiles = bcftoolsFilter_Index.out.vcfOnly.collect()
    mergeBCFtools(bcftools_allFiles, bcftoolsFilter_Index.out.vcfIndex.collect())
    variantsToTable_bcftools(mergeBCFtools.out.bcftoolsMerged)
    bcftoolsBarcode(variantsToTable_bcftools.out.bcftoolsTable, bcftools_allFiles)
    predictGeo_bcftools(bcftoolsBarcode.out.barcodeOut)
    mergePredictions(predictGeo_gatk.out.gatkCountry, predictGeo_gatk.out.gatkRegion, predictGeo_bcftools.out.bcfCountry, predictGeo_bcftools.out.bcfRegion, gatkBarcode.out.barcodeOut, bcftoolsBarcode.out.barcodeOut)

    // marsVOI_report(mergePredictions.out.regionPredictionMerge, Summary.out.Summary_Snp)

//  Clustering 
  //workflow 1 - haplotype calling + haplotype sheet
  // myReads = Channel.fromFilePairs("${params.fastqDir}/*_{R1,R2}_*.fastq.gz")
  newHaps = Channel.fromPath("${params.origNewHaps}/*.fasta")

  Channel
    .fromPath("${params.loci}")
    .splitText()
    .set{ loci_ch }

  // bbduk(myReads)
  // fastqc_beforeMerge(myReads)

  // fastqc_afterMerge(bbduk.out.readyReads)

  makeRefs_v2()
  finalNewRefs_v2(makeRefs_v2.out.completeFile.collectFile())
  //
  // mapping(finalNewRefs_v2.out.outFinalRefs,Trim_reads.out.readyReads)
  mapping(finalNewRefs_v2.out.outFinalRefs, noHumanReads.out.humanRemoveFASTQ)

  mapping2fastq(mapping.out.mappedOnlySam)
  bowtie2(mapping2fastq.out.mappedFastq)

  bedExtraction_plusNewHaps(loci_ch, bowtie2.out.finishSorting, finalNewRefs_v2.out.outShortRefs)
  makeFinalReads_v2(bedExtraction_plusNewHaps.out.tryBedFastq.groupTuple())
  
  nameNewHaps(bedExtraction_plusNewHaps.out.allNewHaps.collectFile(), finalNewRefs_v2.out.outShortRefs)

  makeRefs_assignHaps2_v2()

  finalNewRefs_assignHaps2(makeRefs_v2.out.completeFile.collectFile(), nameNewHaps.out.newHaps_singleFile.collectFile().ifEmpty("$baseDir/modules/files_forNextflow/empty.fasta"))
  blastDB(finalNewRefs_assignHaps2.out.outFinalRefs)
  
  firstClustering_v2(makeFinalReads_v2.out.finalClipped)
  blastMulti_v2(firstClustering_v2.out.multiTuple, firstClustering_v2.out.depthTuple, blastDB.out.outFinalRefs_db)
  hapsBlast(blastMulti_v2.out.blastFasta, blastDB.out.outFinalRefs_db)

  // place the genotypes in the correct folder with the correct name
  copyGenotypes(hapsBlast.out.genotypeNoJunction.ifEmpty("$baseDir/modules/files_forNextflow/emptySNP_genotype").collectFile())

  // create a haplotype sheet and copy the newly identified haplotypes to the correct folder
  haplotypeSheet(copyGenotypes.out.endFile.collectFile())

  //workflow 2 and 3 - matrix + tree
  chromosomes = Channel.fromList(['Nu_CHR01','Nu_CHR02','Nu_CHR03','Nu_CHR04','Nu_CHR05','Nu_CHR06','Nu_CHR07','Nu_CHR08','Nu_CHR09','Nu_CHR10','Nu_CHR11','Nu_CHR12','Nu_CHR13','Nu_CHR14'])
  // updatedDB = Channel.fromPath("$baseDir/treeFolder/*geoPrediction_simple.txt").toSortedList().flatten().last()
  //hapSheets = Channel.fromPath("$baseDir/haplotype_sheets/*.txt").toSortedList().flatten().last()
  //paramInput_wf2 = Channel.fromPath(params.wf2HapSheet)

  individualChroms(haplotypeSheet.out.finalHapSheetLarge, chromosomes)
  individualMats(individualChroms.out.chromHapSheet)
  scaleMatrices(individualMats.out.matrixOut.collect(), individualChroms.out.chromHapSheet.collect())
  runModule3_Pvivax(scaleMatrices.out.finalMatrix)
  runModule3_parnas(runModule3_Pvivax.out.tree, runModule3_Pvivax.out.threshold)


}