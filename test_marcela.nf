#!/usr/bin/env nextflow

/* Executar antes no terminal: 
- Para criar a imagem docker do QIIME2: docker pull quay.io/qiime2/core:2023.9
- Para instalação do QIIME2, seguir as orientações de: https://docs.qiime2.org/2024.5/install/native/#install-qiime-2-within-a-conda-environment
- Para ativar o ambiente conda do QIIME2: conda activate qiime2-amplicon-2024.5 
- Para rodar o script abaixo  use: nextflow run test_marcela.nf -profile docker
- Para visualizar os arquivos qzv gerados, realizar o download dele e usar em https://view.qiime2.org/ */


//Parâmetros
params.reads = "/workspaces/training/nf-training/data/data_Marcela/SRR*_16S_L001_R{1,2}_001.fastq.gz"
params.outdir = "results"

reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)

process IMPORT {
  tag "${sample_id}"
  publishDir "${params.outdir}/import", mode: 'copy'

  input:
  tuple val(sample_id), path(reads)

  output:
  path "demux.qza", emit: demux

  script:
  """
  R1_PATH=\$(readlink -f ${reads[0]})
  R2_PATH=\$(readlink -f ${reads[1]})
  
  mkdir -p manifest
  echo "sample-id,absolute-filepath,direction" > manifest/manifest.csv
  echo "${sample_id},\${R1_PATH},forward" >> manifest/manifest.csv
  echo "${sample_id},\${R2_PATH},reverse" >> manifest/manifest.csv

  qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path manifest/manifest.csv \
    --output-path demux.qza \
    --input-format 'PairedEndFastqManifestPhred33'
  """
}

process DADA2 {
  cpus 4 
  publishDir "${params.outdir}/dada2", mode: 'copy'
  
  input:
  path demux

  output:
  path "table.qza", emit: table
  path "rep-seqs.qza", emit: rep_seqs
  path "stats.qza", emit: stats

  script:
  """
  qiime dada2 denoise-paired \
    --i-demultiplexed-seqs ${demux} \
    --p-trunc-len-f 220 \
    --p-trunc-len-r 200 \
    --p-n-threads ${task.cpus} \
    --o-table table.qza \
    --o-representative-sequences rep-seqs.qza \
    --o-denoising-stats stats.qza
  """
}

process QZV {
  publishDir "${params.outdir}/visualizations", mode: 'copy'
  
  input:
  path demux
  path table
  path rep_seqs
  path stats

  output:
  path "*.qzv", emit: visualizations
  

  script:
  """
  qiime demux summarize \
    --i-data ${demux} \
    --o-visualization demux.qzv

  qiime feature-table summarize \
    --i-table ${table} \
    --o-visualization table.qzv

  qiime feature-table tabulate-seqs \
    --i-data ${rep_seqs} \
    --o-visualization rep-seqs.qzv

  qiime metadata tabulate \
    --m-input-file ${stats} \
    --o-visualization stats.qzv
  """
}

workflow {
  IMPORT(reads_ch)
  DADA2(IMPORT.out.demux)
  QZV(IMPORT.out.demux,DADA2.out.table, DADA2.out.rep_seqs, DADA2.out.stats)
}