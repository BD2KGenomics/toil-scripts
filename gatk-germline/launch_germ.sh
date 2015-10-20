#!/usr/bin/env bash
export TMPDIR=/home/jacob/fake_mnt
python germ.py \
$TMPDIR/jstore \
--config "config.txt" \
--reference "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz" \
--phase "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/1000G_phase1.indels.b37.vcf.gz" \
--mills "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz" \
--dbsnp "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.vcf.gz" \
--hapmap "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/hapmap_3.3.b37.vcf.gz" \
--omni "ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/1000G_omni2.5.b37.vcf.gz" \
--output_dir '/home/jacob/' \
