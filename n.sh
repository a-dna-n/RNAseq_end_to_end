#!/bin/bash
log() { echo "$*" >&2 ; }
refdir=${HOME}/star/ref
dna=$refdir/Mus_musculus.GRCm39.dna.primary_assembly.fa
rna=$refdir/Mus_musculus.GRCm39.109.gtf
index=${HOME}/star/mouse.GRCm39_109.150
nt=150
cpus=$(grep -c ^processor /proc/cpuinfo)
maxram=48000000000
if [[ -d $index ]]; then
    echo "$index exists" >&2  && exit 1
fi
for file in $dna $rna ; do
    if [[ ! -e $file ]]; then
        if [[ -e $file.gz ]]; then
            gunzip --keep $file.gz
        else
            echo "$file not found" >&2  && exit 1
        fi
    fi
done
STAR --runMode genomeGenerate --genomeDir $index \
    --genomeFastaFiles $dna --sjdbGTFfile $rna --sjdbOverhang $nt \
    --runThreadN $cpus --limitGenomeGenerateRAM $maxram > $index.log 2> $index.err

