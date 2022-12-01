#!/bin/env hash 


POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -g|--gene)
      GENE="$2"
      shift # past argument
      shift # past value
      ;;
    -b|--bedfile)
      BEDFILE="$2"
      shift # past argument
      shift # past value
      ;;
    -r|--reference)
      REFERENCE="$2"
      shift # past argument
      shift # past value
      ;;
    -t|--temp)
      TEMPDIR="$2"
      shift # past argument
      shift # past value
      ;;
    --default)
      DEFAULT=YES
      shift # past argument
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

echo "GENE NAME = ${GENE}"
echo "BED FILE = ${BEDFILE}"
echo "REFERENCE     = ${REFERENCE}"
echo "TEMPDIR         = ${TEMPDIR}"
#echo "Number files in SEARCH PATH with EXTENSION:" $(ls -1 "${SEARCHPATH}"/*."${EXTENSION}" | wc -l)

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi



data=${REFERENCE}
temp_output_dir=${TEMPDIR}
inputBED12=${BEDFILE}
gene_name=${GENE}

#todo make scriptable
#python get_bed.py |grep "ENST00000331340.8"  > ikzf1.bed
#sed 's/ \+ /\t/g'  ikzf1.bed  > ikzf1_tab.bed

mkdir -p $temp_output_dir

python create_bedfiles.py $inputBED12 $temp_output_dir 

find $temp_output_dir  -name "*.bed"  -type f -exec cat {} + >  ${gene_name}_combined.BED12

bedtools getfasta -fi $data -bed ${gene_name}_combined.BED12 -split -name > ${gene_name}_toblerone_transcriptome.fasta

sed -e '/^>.*/s/|/~/g' -e "/^>.*/s/$/|${gene_name}|3|4|5|6|7|8|9/g" ${gene_name}_toblerone_transcriptome.fasta > ${gene_name}_toblerone_transcriptome_mod.fasta
sed -e 's/::chr.*:[0-9]*-[0-9]*|/|/g' ${gene_name}_toblerone_transcriptome_mod.fasta > ${gene_name}_toblerone_transcriptome_input.fasta  



#tinyt index -i toblerone_transcriptome_cd22.tidx  toblerone_transcriptome_cd22_mod 

echo toblerone  index -i  ${gene_name}_toblerone_transcriptome_input.idx  ${gene_name}_toblerone_transcriptome_input.fasta


