#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -l h_vmem=2G

line_num=1
echo "Gene_1_name	Gene_1_breakpoint	pre_breakpoint_diff	pre_breakpoint_pval	post_breakpoint_diff	post_breakpoint_pval	log2_diff	Mean_abs_expr	Gene_2_name	Gene_2_breakpoint	pre_breakpoint_diff	pre_breakpoint_pval	post_breakpoint_diff	post_breakpoint_pval	log2_diff	Mean_abs_expr	sample_name" > all_fusions_results.tsv
cat $1 | awk '{print $31,$38,$32,$39,$46}' | sort -u | while read line
do 
line_num=`expr $line_num + 1`
sample=`echo $line | awk '{print $5}'| sed s/[A-Z]*_[0-9]*\.// | sed s/\.sorted_genome_alignments.bam_out//`
sample=`grep $sample TCGA_all_samples_data.tsv | awk '{print $1}'`
gene1=`echo $line | awk '{print $1}'`
gene2=`echo $line | awk '{print $3}'`
break1=`echo $line | awk '{print $2}'`
break2=`echo $line | awk '{print $4}'`

#if [ `grep -c $gene1 unc_v2_exon_hg18_probe` == 0 -o `grep -c $gene2 unc_v2_exon_hg18_probe` == 0 ]
#then
#echo $gene1 "or" $gene2 "not in TCGA exon counts file"
#continue
#fi
#
#echo "processing" $gene1 "and" $gene2 "on line" $line_num
#if [ -f $gene1"_data.tsv" ]
#then
#echo "file exists	 continuing..."
#else
#grep $gene1 unc_v2_exon_hg18_probe | awk '{print $1}' | while read line1; do grep $line1 HiSeqV2_exon ; done > $gene1"_data.tsv"
#fi
#
#if [ -f $gene2"_data.tsv" ]
#then
#echo "file exists	 continuing..."
#else
#grep $gene2 unc_v2_exon_hg18_probe | awk '{print $1}' | while read line1; do grep $line1 HiSeqV2_exon ; done > $gene2"_data.tsv"
#fi


#echo $sample
Rscript fusion_script_with_function_aesthetic.R $gene1 $break1 $gene2 $break2 $sample all_fusions_results.tsv
done
echo -e "\n" >> all_fusions_results.tsv



#pval=0.00001
#log2=0.58
#echo "Gene_1_name	Gene_1_breakpoint	pre_breakpoint_diff	pre_breakpoint_pval	post_breakpoint_diff	post_breakpoint_pval	log2_diff	Gene_2_name	Gene_2_breakpoint	pre_breakpoint_diff	pre_breakpoint_pval	post_breakpoint_diff	post_breakpoint_pval	log2_diff	sample_name" > filtered_hits.tsv

#awk -v pv=$pval '($4<pv)||($6<pv)||($12<pv)||($14<pv){print}' all_fusions_results.tsv | awk -F'\t' -v l2=$log2 'function abs(x){return ((x < 0.0) ? -x : x)} {if (\
#((abs($7)>l2)&&($7!="NaN"))||\
#((abs($15)>l2)&&($15!="NaN"))) print}'  >> filtered_hits.tsv

