#!/bin/bash

RUN_DIR=$PWD
DATADIR=$RUN_DIR

usage (){
  echo "bash $(readlink -f $0) <ID>"
  exit 1
}


# inspect run script location
FD_REPO=$(readlink -f $(dirname "$0"));
CODEDIR=$FD_REPO/src/complex_assembly

# database dir
DATABASE_DIR="/mnt/db/"
ID=$1
if [[ "$ID" == "" || ! -f "${DATADIR}/${ID}_useqs.csv" || ! -f "${DATADIR}/${ID}_chains.csv" ]];then
  echo "ID is not valid \$1 or files not found ${DATADIR}/${ID}_useqs.csv"
  usage
fi

# Path and user config
# change here if you wish to use a alternate version of any database.

uniclust30_database_path="$DATABASE_DIR/uniref30_uc30/UniRef30_2020_06/UniRef30_2020_06"

# It seems like we already have a fit conda env.
source activate alphafold
set -e
#############PARAMETERS#############
 #Where all scripts are run from, now the current directory


HHBLITS=$(which hhblits) #Path to hhblits
HHBLITSDB=$uniclust30_database_path

###Note!
#This runscript assumes that you dont have singularity in your path

#Run the assembly pipeline starting from predicted interactions
#########INPUTS and PATHS##########
USEQS=$DATADIR/$ID'_useqs.csv'
CHAINS=$DATADIR/$ID'_chains.csv'
if [[ ! -f $DATADIR/$ID'_ints.csv' ]];then
  INTERACTIONS='' #Leave empty if the interactions are not known - here they are not used. See the file $DATADIR/$ID'_ints.csv' for how to supply such a file
else
  INTERACTIONS=$DATADIR/$ID'_ints.csv'
fi

SUBSIZE=3 #What order the subcomplexes should be (2 or 3)
GET_ALL=1 #If to get all interactions (1) or not (0) - when the interactions are known
#########OUTPUTS#########
#The best assembled complex ,the assembly path and its score (mpDockQ)

#########Step1: MSA PIPELINE#########
#SINGIMG=$BASE/src/AF2/AF_environment.sif #Sing img

MSADIR=$DATADIR/hhblits
mkdir -p $MSADIR
#Write individual fasta files for all unique sequences
cmd="python ${FD_REPO}/src/preprocess/write_hhblits_fasta.py --unique_seq_df ${USEQS} --outdir ${MSADIR}/"
echo "$cmd"
eval "$cmd"

#Run HHblits
for file in $MSADIR/*.fasta
do
  SUBID=$(awk -F'/' '{print $NF}' <<< "echo $file")
  SUBID=$(echo $SUBID|cut -d '.' -f 1)
  if test -f $MSADIR/$SUBID'.a3m'; then
    echo $SUBID MSA exists
  else
    cmd="${HHBLITS} -i ${file} -d ${HHBLITSDB} -E 0.001 -cpu 32 -all -n 2 -oa3m ${MSADIR}/${SUBID}.a3m"
    echo "$cmd"
    eval "$cmd"
  fi
done

#########Step2: FOLDING PIPELINE#########
wait
#Write the Paired and Block Diagonalized MSAs to predict sub-components
cmd="python ${FD_REPO}/src/preprocess/prepare_folddock_run.py --msadir ${MSADIR}/ \
--complex_id ${ID} \
--outdir ${MSADIR}/ \
--useqs ${USEQS} \
--interactions ${INTERACTIONS} \
--intchain2seq ${CHAINS} \
--get_all ${GET_ALL} \
--subsize ${SUBSIZE} "

echo "$cmd"
eval "$cmd"

#Create structure dir
STRUCTURE_DIR=$DATADIR/AF
mkdir -p $STRUCTURE_DIR
#Get subids
head -n 1 $MSADIR/*.fasta|grep '|'|cut -d '|' -f 1|cut -d '>' -f 2 > $DATADIR/sub_ids.txt
#Get lens for AF run
head -n 1 $MSADIR/*.fasta|grep '|'|cut -d '|' -f 2|cut -d '-' -f 1-2 > $DATADIR/lens.txt

#Predict the subcomponents
##### AF2 CONFIGURATION #####
PARAM=$DATABASE_DIR/alphafold/2022-03-02/
PRESET='full_dbs' #Choose preset model configuration - no ensembling (full_dbs) and (reduced_dbs) or 8 model ensemblings (casp14).
MAX_RECYCLES=10 #max_recycles (default=3)
MODEL_NAME='model_1' #model_1_ptm

#Go through all subcomponents and predict their structure
NCOMPONENTS=$(wc -l < $DATADIR/sub_ids.txt)
for ((LN=1;LN<=NCOMPONENTS;LN++))
do
  SUBID=$(sed -n $LN'p' $DATADIR/sub_ids.txt)
  echo $SUBID
  ####Get fasta file####
  FASTAFILE=$MSADIR/$SUBID'.fasta'
  ####Get chain break#### Note! This is now set for trimer subcomponents
  CB=$(sed -n $LN'p' $DATADIR/lens.txt)
  CB1=$(echo $CB|cut -d '-' -f 1)
  CB2=$(echo $CB|cut -d '-' -f 2)
  CB2=$(( $CB1 + $CB2 ))
  CB=$CB1,$CB2
  ####Get MSAs####
  #HHblits paired
  PAIREDMSA=$MSADIR/$SUBID'_paired.a3m'
  ##HHblits block diagonalized
  BLOCKEDMSA=$MSADIR/$SUBID'_blocked.a3m'
  MSAS="$PAIREDMSA,$BLOCKEDMSA" #Comma separated list of msa paths
  #Check if prediction exists
  if test -f $STRUCTURE_DIR/$SUBID/*_1.pdb; then
    echo $SUBID prediction exists
  else

    cmd="python ${FD_REPO}/src/AF2/run_alphafold.py \
                --fasta_paths=${FASTAFILE} \
                --msas=${MSAS} \
                --chain_break_list=${CB} \
                --output_dir=${STRUCTURE_DIR} \
                --model_names=${MODEL_NAME} \
                --data_dir=${PARAM} \
                --fold_only \
                --uniref90_database_path=$HHBLITSDB \
                --mgnify_database_path=$HHBLITSDB \
                --bfd_database_path=$HHBLITSDB \
                --uniclust30_database_path=$HHBLITSDB \
                --pdb70_database_path=$HHBLITSDB \
                --template_mmcif_dir=$HHBLITSDB \
                --obsolete_pdbs_path=$HHBLITSDB \
                --preset=$PRESET \
                --max_recycles=$MAX_RECYCLES "
      echo "$cmd"
      eval "$cmd"
  fi

done

#########Step3: ASSEMBLY PIPELINE#########
wait
COMPLEXDIR=$DATADIR/assembly/complex/ #Where all the output for the complex assembly will be directed


#Make complex directory
mkdir -p $COMPLEXDIR
#Rewrite the FoldDock preds to have separate chains according to the fasta file seqlens
cmd="python $CODEDIR/rewrite_fd.py --pdbdir $STRUCTURE_DIR --pdb_id $ID"
echo "$cmd"
eval "$cmd"

#Copy all predicted unique chain interactions to reflect all possible interactions
SUB_PDBDIR=$STRUCTURE_DIR/
OUTDIR=$DATADIR/assembly/
cmd="python $CODEDIR/copy_preds.py --complex_id $ID --pdbdir $SUB_PDBDIR --outdir $OUTDIR \
--useqs $USEQS --subsize $SUBSIZE --interactions $INTERACTIONS --intchain2seq $CHAINS --get_all $GET_ALL"
echo "$cmd"
eval "$cmd"

#Rewrite AF predicted complexes to have proper numbering and chain labels
PDBDIR=$OUTDIR
cmd="python $CODEDIR/rewrite_af_pdb.py --pdbdir $PDBDIR --pdb_id $ID --outdir $OUTDIR"
echo "$cmd"
eval "$cmd"

#Write all pairs
PAIRDIR=$PDBDIR/pairs/
META=$PDBDIR/meta.csv #where to write all interactions
#It is necessary that the first unique chain is named A-..N for and the second N-... and so on
mkdir -p $PAIRDIR
#Glob for all files with each chain in order (A,B,C,D) A-->B,C,D; B--> C,D; C-->D
cmd="python $CODEDIR/write_all_pairs.py --pdbdir $PDBDIR --pairdir $PAIRDIR --meta $META \
--interactions $INTERACTIONS --get_all $GET_ALL"

echo "$cmd"
eval "$cmd"

#####Clean up intermediate files#####
mv $PDBDIR/$ID'_chains.csv' $COMPLEXDIR/$ID'_chains.csv'
rm -r $PDBDIR/$ID*
mv $COMPLEXDIR/$ID'_chains.csv' $PDBDIR/$ID'_chains.csv'

#Assemble from pairs
#Find the best non-overlapping path that connect all nodes using Monte Carlo Tree search
PLDDTDIR=$PDBDIR/plddt/
CHAIN_SEQS=$PDBDIR/$ID'_chains.csv' #Updated chain seqs
cmd="python $CODEDIR/mcts.py --network $META \
--pairdir $PAIRDIR --plddt_dir $PLDDTDIR \
--useqs $USEQS --chain_seqs $CHAIN_SEQS \
--outdir $COMPLEXDIR"
echo "$cmd"
eval "$cmd"

#########Step4: SCORING#########
#Score complex
MODEL=$COMPLEXDIR/best_complex.pdb
MODEL_PATH=$COMPLEXDIR/optimal_path.csv
DT=8
OUTNAME=$COMPLEXDIR/$ID'_score.csv'
cmd="python $CODEDIR/score_entire_complex.py --model_id $ID --model $MODEL \
--model_path $MODEL_PATH \
--useqs $USEQS --chain_seqs $CHAIN_SEQS \
--outname $OUTNAME"

echo "$cmd"
eval "$cmd"

#####Clean up intermediate files#####
rm -r $PLDDTDIR
rm -r $PAIRDIR
