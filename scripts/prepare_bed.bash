#!/bin/env hash 


POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -s|--source)
      SOURCE="$2"
      shift # past argument
      shift # past value
      ;;
    -t|--transcript)
      TRANSCRIPT="$2"
      shift # past argument
      shift # past value
      ;;
    -b|--bedfile)
      BEDFILE="$2"
      shift # past argument
      shift # past value
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

echo "SOURCE BED= ${SOURCE}"
echo "TRANSCRIPT= ${TRANSCRIPT}"
echo "BED FILE = ${BEDFILE}"
#echo "Number files in SEARCH PATH with EXTENSION:" $(ls -1 "${SEARCHPATH}"/*."${EXTENSION}" | wc -l)

if [[ -n $1 ]]; then
    echo "Last line of file specified as non-opt/last argument:"
    tail -1 "$1"
fi



#todo make scriptable
#python get_bed.py |grep "ENST00000331340.8"  > ikzf1.bed
#sed 's/ \+ /\t/g'  ikzf1.bed  > ikzf1_tab.bed


grep ${TRANSCRIPT} < ${SOURCE} > ${BEDFILE}_tmp
sed 's/ \+ /\t/g'  ${BEDFILE}_tmp > ${BEDFILE}



