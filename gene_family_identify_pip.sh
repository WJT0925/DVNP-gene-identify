while getopts "g:q:f:p:c:d:r:I:O:n:t:" OPTION
do
    case $OPTION in
      "g")
        genome=$OPTARG
        ;;
      "q")
        query=$OPTARG
        ;;
      "f")
        gff=$OPTARG
        ;;
      "p")
        pep=$OPTARG
        ;;
      "c")
        extend_cds=$OPTARG
        ;;
      "d")
        extend_pep=$OPTARG
        ;;
      "r")
        rate=$OPTARG
        ;;
      "I")
        identify=$OPTARG
        ;;
      "O")
        overlap=$OPTARG
        ;;
      "n")
        name=$OPTARG
        ;;
      "t")
        cpu=$OPTARG
        ;;
	esac
    #echo "option index is $OPTIND"
done
#默认值
#extend_cds=${c:="2000"}
#extend_pep=${d:="1000"}
#rate=${r:="0.3"}
#identify=${I:="80"}
#overlap=${O:="0.05"}

#query序列比对参考基因组
ref=$(basename $query)
mkdir "$name" 
cd "$name"
perl ~/pip/genefamily_pipline/bin/blast_cut_new.pl $genome $query -eval 1e-5 -cpu $cpu
cat $ref.cut/*m8 > all.blast.m8
perl ~/pip/genefamily_pipline/bin/solar.pl -a prot2genome2 -f m8 all.blast.m8 > all.blast.m8.solar 
perl ~/pip/genefamily_pipline/bin/solar_add_seqlen.pl $genome $query all.blast.m8.solar > all.blast.m8.solar.adden
perl ~/pip/genefamily_pipline/bin/solar_add_identity.pl --solar all.blast.m8.solar.adden --m8 all.blast.m8 > $name.fa.blast.cov
awk -v filt=$rate '$5>filt' $name.fa.blast.cov > $name.fa.blast.cov.filt
perl ~/pip/genefamily_pipline/bin/get_pos.pl $name.fa.blast.cov.filt > $name.pos
perl  ~/pip/genefamily_pipline/bin/extract_sequence.pl --pos  $name.pos --fasta $genome  --extent  $extend_cds > $name.nuc
perl ~/pip/genefamily_pipline/bin/prepare_pep.pl $name.pos $query > $name.pep
awk '{print $1"  "$3}' $name.pos  > $name.strandList
perl ~/pip/genefamily_pipline/bin/call_genewise.pl --pep $name.pep --nuc $name.nuc --list $name.strandList --key $name --out $name  --line 80
parallel -j 30 < ./$name/result/gw.sh 
find ./$name/result/ -type f -name "*.gw" -print| xargs -I {} cat {} >> $name.gw
perl ~/pip/genefamily_pipline/bin/gw_parser_change.pl --gw $name.gw --pep $name.pep --ac 0 --id 0 --type 1 >  $name.gw.alg
perl ~/pip/genefamily_pipline/bin/gw_parser_change.pl --gw $name.gw --pep $name.pep --ac 0 --id 0 --type 2 >  $name.gw.mut
perl  ~/pip/genefamily_pipline/bin/merge_overlap.pl $name.gw.alg 0.1 > $name.gw.alg.nr
perl  ~/pip/genefamily_pipline/bin/gw_to_gff.pl $name.gw $name.pep > $name.gw.gff
perl  ~/pip/genefamily_pipline/bin/Select_best_hit_genewise.pl $name.gw.mut $name.gw.alg $name.gw.gff $name.best
perl  ~/pip/genefamily_pipline/bin/frameChg.pl $name.best.genewise.gff $name.best.genewise.md.gff
perl  ~/pip/genefamily_pipline/bin/elongate.pl $genome $name.best.genewise.gff $name.best.genewise.elongate.gff -extend $extend_pep
perl ~/pip/genefamily_pipline/bin/chgName.pl $name.best.genewise.elongate.gff $name.best.genewise.final.gff
perl ~/pip/genefamily_pipline/bin/fetchCds.pl -fa $genome -gff $name.best.genewise.final.gff -out $name.best.genewise.final.cds.fa
perl ~/pip/genefamily_pipline/bin/cds2aa.pl -cds $name.best.genewise.final.cds.fa -out $name.best.genewise.final.pep.fa
perl ~/pip/genefamily_pipline/bin/addGeneToGff.pl $name.best.genewise.final.gff tmp.gff
perl ~/pip/genefamily_pipline/bin/checkCDS.pl -genome $genome -gff tmp.gff -pep $name.best.genewise.final.pep.fa -out tmp.chk
makeblastdb -in $pep -dbtype prot -title tarpep.fa -out tarpep.fa -logfile log
blastp -query $name.best.genewise.final.pep.fa -db tarpep.fa -evalue 1e-5 -outfmt 6 -out tarpep.blast.m8 -num_threads $cpu
perl ~/pip/genefamily_pipline/bin/extract_best_alin.pl tarpep.blast.m8 tarpep.blast.best.m8
perl ~/pip/genefamily_pipline/bin/gff_annotation_genewise.pl $gff $name.best.genewise.final.gff tarpep.blast.best.m8 $overlap gene_adjusted_inf.list

#Function Annotation (Note:Check seqID format)
#perl ~/pip/genefamily_pipline/bin/extra_seq.pl merge_gene.list $pep final_adijust.fiter.pep.fa
#perl ~/pip/genefamily_pipline/bin/extra_seq.pl genewise_gene.list $name.best.genewise.final.pep.fa  final_adijust.fiter.pep.fa
#/home/jingtian/software/miniconda3/bin/hmmscan -o $name.pfam.out --tblout  $name.pfam.out --noali -E 1e-5 /home/jingtian/Database/HMM/Pfam-A.hmm final_adijust.fiter.pep.fa
#/home/jingtian/software/miniconda3/bin/hmmscan --domtblout $name.pfam.out --noali -E 1e-5 /home/jingtian/Database/HMM/Pfam-A.hmm final_adijust.fiter.pep.fa
#diamond blastp -p $cpu -d ~/Database/nr/diamond/nr_2021.dmnd -q final_adijust.fiter.pep.fa -o $name_adijust.fiter.pep.NR.out -f 6 qseqid qlen sseqid length evalue score stitle -e 1e-5 -b1 -c8 -k 1
