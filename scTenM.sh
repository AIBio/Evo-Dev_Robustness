#!/usr/bin/env bash

### ---------------
### Whole procedure
### ---------------

# 1. Build the working directory (required);
# 2. Prepare the genome files for analysis (required); 
# 3. Cell Ranger pipeline (required);
# 4. Quantify repeat element (optional);
# 5. Velocyto analyssi (optional);

### ---------------
### Running command
### ---------------

# nohup bash scTenM.sh -c scTenM.cfg -m prep/index/gene/repeat/velocyto/all &

### --------------
### Parse augments
### --------------

Help()
{
cat >&2 <<EOF
[FUNCTION:]
        Custome pipeline to perform CellRanger pipeline for RNA, ATAC and multiome data.

[USAGE:]
        ${0##*/} -c <configuration file> -h
        -c <*.cfg file>       : configuration file, includes variables and functions;
        -m <running mode>     : prep/index/gene/repeat/velocyto/all
        -h <help document>    : help document;

[EXAMPLE COMMAND:]
        nohup bash ${0##*/} -c scTenM.cfg -m prep/index/gene/repeat/velocyto/all &
        * The steps in "all" mode: prep -> index -> gene -> repeat -> velocyto

[SOLVING CONDA ENV:]
        * For whole pipeline
        $ conda create -n sc.p37 python=3.7
        $ conda activate rna.p37
        $ conda install -c bioconda 

[CONF FILE FORMAT:]
        1. workding directory
        2. full path of raw data
        3. txt file that records sample names
        4. data type: RNA, ATAC or Multiome
        5. gene gtf file: well-processed GTF file
        6. genomic sequence fa file
        7. cell ranger index (directory)
        8. motif file
        9. species
        10. genome name
        11. genome version
        12. library csv file for multiome
        13. velocyto shell script
        14. repeat bed
        15. bam string to filter samples and then run scTE
EOF
   exit 1
}
while getopts "c:m:h" opt
do
    case "$opt" in
        c ) cfg="$OPTARG" ;;
        m ) mode="$OPTARG" ;;
        h ) Help ;;
    esac
done
if [ -z "${cfg}" ]
then
    echo -e "[$(date)] [Error] Some or all of the parameters are empty"
    Help
fi
echo "> ---------------------------- <"
source ${cfg}
declare -a vars
declare -a labels
vars=("wd" "rawdata" "samples" "dt" "gtf" "fa" "cr_indx" "motif" "species"
      "gname" "gversion" "library" "run_velocyto" "repeat" "bam_str")
labels=("[1] working directory" "[2] directory of raw fastq files"
        "[3] sample file" "[4] data type" "[5] gene gtf file"
        "[6] genomic sequence fa file" "[7] cell ranger index" 
        "[8] motif file" "[9] species" "[10] genome name" 
        "[11] genome version" "[12] library csv file for multiome" 
        "[13] velocyto shell script" "[14] repeat bed"
        "[15] bam string to run scTE")
for i in ${!vars[@]}
do
    var="${vars[$i]}"
    label="${labels[$i]}"
    echo -e "${label}: ${!var}"
done
echo "> ---------------------------- <"

### ------------
### Run pipeline
### ------------

if [ "${mode}" == "prep" ] || [ "${mode}" == "all" ]
then
    echo -e "[$(date)] [Info] 1st step: Organize the workding directory"
    echo -e "[$(date)] [Info] Begin ~~~"
    mkdir ${wd}/{logs,results,rawdata,metadata,scripts}
    echo -e "[$(date)] [Info] Finished ~~~"
fi



if [ "${mode}" == "index" ] || [ "${mode}" == "all" ]
then
    echo -e "[$(date)] [Info] 2nd step: Prepare the genome files for analysis"
    echo -e "[$(date)] [Info] Begin ~~~"
    if [ "${dt}" == "ATAC" ]
    then
        if [ ! -e ${cr_indx} ] || [ -z "$(ls -A -- ${cr_indx})" ]
        then
            RefCRatac ${fa} ${gtf} ${motif} ${species} ${gname} ${gversion} ${cr_indx} ${cr_index}/../
        else
            echo "ATAC CellRanger genome index exists"
        fi
    elif [ "${dt}" == "RNA" ]
    then
        if [ ! -e ${cr_indx} ] || [ -z "$(ls -A -- ${cr_indx})" ]; then
            RefCRrna ${fa} ${gtf} ${gname} ${gversion} ${cr_indx} ${cr_index}/../
        else
            echo "RNA CellRanger genome index exists"
        fi
    elif [ "${dt}" == "Multiome" ]
    then
        if [ ! -e ${cr_indx} ] || [ -z "$(ls -A -- ${cr_indx})" ]
        then
            RefCRarc ${fa} ${gtf} ${motif} ${species} ${gname} ${gversion} ${cr_indx} ${cr_index}/../
        else
           echo "Multiome CellRanger genome index exists"
        fi
    else
        echo "Unrecognized data type!"; exit 1
    fi
    echo -e "[$(date)] [Info] Finished ~~~"
fi



if [ "${mode}" == "gene" ] || [ "${mode}" == "all" ]
then
    echo -e "[$(date)] [Info] 3rd step: CellRanger count"
    echo -e "[$(date)] [Info] Begin ~~~"
    # - extract directory name
    ls ${rawdata} | grep "fastq.gz$" \
        | sed 's,_S.*[0-9]_L00[0-9]_.._00[0-9].fastq.gz,,g' \
        | uniq > ${wd}/metadata/sample_key.txt
    # - run CellRanger
    case "$dt" in
        ATAC)
            echo "Run CellRanger-atac count module for scATAC-seq data"
            CRatac ${rawdata} ${wd}/metadata/sample_key.txt \
                   ${wd}/results/cratac/count ${wd}/logs/cratac/count ${cr_indx}/${gname}
            ;;
        RNA)
            echo "Run CellRanger-rna count module for scRNA-seq data"
            CRrna ${rawdata} ${wd}/metadata/sample_key.txt \
                  ${wd}/results/crrna/count ${wd}/logs/crrna/count ${cr_indx}/${gname}
            ;;
        Multiome)
            echo "Run CellRanger-arc count module for scATAC+RNA data"
            CRarc ${samples} ${wd}/results/crarc/count ${wd}/logs/crarc/count \
                  ${cr_indx}/${gname} ${library}
            ;;
        *)
            echo "Unknown library type to precess"
            exit 1
            ;;
    esac
    # - move directory
    mkdir -p ${wd}/results/crrna/count/${gname}
    cat ${wd}/metadata/sample_key.txt | while read dir
    do
        mv ${dir} ${wd}/results/crrna/count/${gname}/
    done
    echo -e "[$(date)] [Info] Finished ~~~"
fi



if [ "${mode}" == "repeat" ]
then
    echo -e "[$(date)] [Info] 4th step: scTE pipeline"
    echo -e "[$(date)] [Info] Begin ~~~"
    # build index
    name=$(basename -s ".bed" ${repeat})
    idxfile=${wd}/metadata/${name}.exclusive.idx
    if [ ! -e "${idxfile}" ]; then
       scTE_build -te ${repeat} -gene ${gtf} -o custom
       mv custom.exclusive.idx ${idxfile}
    fi
    # quantify
    case "$dt" in
        RNA)
            echo "Run scTE to quantify repeat elements with scRNAseq data"
            logdir=${wd}/logs/scTE
            [ ! -d "${logdir}" ] && mkdir -p ${logdir}
            outdir=${wd}/results/scTE/${name}
            [ ! -d "${outdir}" ] && mkdir -p ${outdir}
            find ${wd}/results/crrna/count/${gname} -name 'possorted_genome_bam*' \
                | grep "bam$" \
                | grep "${bam_str}" \
                | awk -F '/' '{print $(NF-2)"\t"$NF}' > ${wd}/metadata/outs_bam_prefix.txt
            while read prefix bam
            do
                bam_dir=${wd}/results/crrna/count/${gname}/${prefix}/outs
                bam_name=$(basename -s ".bam" ${bam})
                samtools view ${bam_dir}/${bam} -h \
                    | awk '/^@/ || /CB:/' \
                    | samtools view -h -b > ${bam_dir}/${bam_name}_clean.bam
                scTE -i ${bam_dir}/${bam_name}_clean.bam \
                     -p 16 -o ${prefix} -x ${idxfile} --hdf5 False \
                     -CB CB -UMI UB --expect-cells 60000 > ${logdir}/${prefix}.log 2>&1
                ls -lh ${bam_dir}/${bam_name}_clean.bam
                mv ${prefix}.csv ${outdir}/
            done < ${wd}/metadata/outs_bam_prefix.txt
            ;;
        ATAC)
            echo "Run scTE to quantify repeat elements with scATACseq data"
            
            ;;
        Multiome)
            echo "Run scTE to quantify repeat elements with scATAC+RNA data"
            
            ;;
        *)
            echo "Unknown library type to precess"
            exit 1
            ;;
    esac
    echo -e "[$(date)] [Info] Finished ~~~"
fi



if [ "${mode}" == "velocyto" ]
then
    echo -e "[$(date)] [Info] 5th step: Velocyto pipeline"
    echo -e "[$(date)] [Info] Begin ~~~"
    cat ${wd}/metadata/sample_key.txt | while read dir
    do
        indir=${wd}/results/crrna/count/${gname}/${dir}
        bam=${wd}/results/crrna/count/${gname}/${dir}/outs/possorted_genome_bam.bam
        outdir=${wd}/results/crrna/count/${gname}/${dir}/outs/velocyto
        bash ${run_velocyto} -d ${indir} -b ${bam} -g ${gtf} -m "run10x" -o ${outdir}
    done
    echo -e "[$(date)] [Info] Finished ~~~"
fi
