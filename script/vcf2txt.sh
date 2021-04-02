source /mnt/disk2_workspace/chenlingxi/workspace/Bio_Projects/10X_Pipeline/sourceme.bash
python /mnt/disk2_workspace/chenlingxi/workspace/Bio_Projects/10X_Pipeline/tenxtools/sv/sv.py all \
    --sample=$1 \
    --tool=svaba \
    --in_fn=$1.svaba.sv.vcf
