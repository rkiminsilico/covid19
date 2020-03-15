# covid19

Tool for multiple sequence alignment of the proteins of COVID19 (Sars-CoV2)

Requirements:
* Python 3.7 or higher
* pip install aiohttp

Usage:
1. Download Sars-CoV2 protein sequences at https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&VirusLineage_ss=SARS-CoV-2,%20taxid:2697049. Click "Download" button and choose "Sequence data (FASTA Format) Protein" at Step 1, "Download All Records" at Step 2, and "Use default: Accession GenBank Title" at Step 3. Then, click "Download" button.
2. The filename of the downloaded file is usually sequences.fasta.
3. Run `python covid19.py -i sequences.fasta` Optionally, aligning specific proteins can be skipped with `--skip` option.
4. Alignment will be done for each protein using the included clustalo and a browser page will be launched to view the multiple sequence alignments of the proteins.
