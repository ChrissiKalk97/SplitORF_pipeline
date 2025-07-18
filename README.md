# Split-ORF pipeline

Split Open Reading frames (Split-ORFs) exist on transcripts that encode two or more open reading frames.
These multiple open reading frames each encode a part of a full length protein, splitting it into pieces. 

The Split-ORF pipeline predicts possible Split-ORF transcripts in a user defined set of transcript sequences.
Additionally, with DNA and protein unique regions are predicted for the Split-ORFs. These regions do not occur in any
other annotated protein coding transcript. 

## Input Files

- Transcripts sequences for Split-ORF prediction (FASTA)
- Reference protein coding transcript sequences (DNA sequences, FASTA)
- Reference protein sequences (amino acid sequences, FASTA)
- PFAM annotation (with columns: transcript ID, start position annotation, stop position annotation, PFAM ID, TSV)
- Exon coordinates (can be downloaded from Ensembl Biomart Strucutres with the following columns:
  Gene stable ID,	Transcript stable ID,	Exon region start (bp),	Exon region end (bp),	Transcript start (bp),
  	Transcript end (bp),	Strand,	Chromosome/scaffold name, TSV)
- Alignment algorithm (blast or diamond)


### How to get Input files
All of the Input files can be downloaded directly from Ensembl Biomart.

1) Download cDNA and amino acid protein coding sequences sequences:
Filters: Genes -> Transcript Type -> protein_coding
Attributes: Sequences -> select cDNA sequences or peptide
Header information: gene stable ID, transcript stable ID


2) download transcripts for SO prediction, here with nonsense mediated decay transcripts:
Filters: Genes -> Transcript type -> nonsense_mediated_decay
Attributes: Sequences -> select cDNA sequences
Header information: gene stable ID, transcript stable ID

3) download PFAM annotation in tsv format:
Attributes: Features -> Transcripts stable ID, PFAM start, PFAM end, PFAM ID
Protein Domains and Families: Pfam ID

4) Exon coordinates:
Attributes: Structures -> Gene stable ID, Transcript stable ID, Exon region start (bp),	Exon region end (bp),	Transcript start (bp),
  	Transcript end (bp), Strand, Chromosome/scaffold name



## Run pipeline
```bash
bash run_pipeline_gen_coords_with_ORFs_21_03_25.sh \
Reference_protein_sequences .fa \
Transcripts_sequences_for_Split-ORF_prediction.fa \
PFAM_annotation.bed \
Reference_protein_coding_transcript_sequences.fa \
Exon_coordinates.bed \
blast
