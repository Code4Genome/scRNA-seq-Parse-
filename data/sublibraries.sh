#Parse pipeline
#!/bin/bash
  
1. First request resources from SLURM 

salloc --partition express --nodes 4 --cpus-per-task 36 -t 03:00:00

2. Then, activate conda enviroment

conda activate spipe-v1.0.0

3. Run the commands: (generate indexed genome)

split-pipe \
    --mode mkref \
    --genome_name hg38 \
    --output_dir /home/c/c/newvolume_test/Genomes/hg38 \
    --genes /home/c/c/newvolume_test/Genomes/Homo_sapiens.GRCh38.109.gtf.gz \
    --fasta /home/c/c/newvolume_test/Genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

----run------

Run it in partition --normal

split-pipe \ 
   -m all \
   --chemistry v2 \
   --nthreads 36 \
   --fq1 /home/c/c/newvolume_test/parse_data/Parse_2_R1_001.fastq.gz \
   --output_dir /home/c/newvolume_test/Analysis/S2_out/ \
   --sample 41_1_C A1 \
   --sample 48_1_C A2 \
   --sample 11_1_C A3 \
   --sample 12_1_C A4 \
   --sample 16_1_C A5\
   --sample 10_1_C A6 \
   --sample 14_3_C A7 \
   --sample 98_1_C A8 \
   --sample 10_1_C A9 \
   --sample 93_4_C A10 \
   --sample 81_4_C A11 \
   --sample 82_6_C A12 \
   --genome_dir /home/c/c/newvolume_test/Genomes/hg38/

Combined mode:

split-pipe \
    --mode comb \
    --output_dir /home/c/c/newvolume_test/Analysis/combined/ \
    --sublibraries \
        /home/c/c/newvolume_test/Analysis/S1_out/ \
        /home/c/c/newvolume_test/Analysis/S2_out/
