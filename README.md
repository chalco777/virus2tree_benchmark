
# Virus2tree pipeline Benchmarking

This notebook is about the benchmarking of [virus2tree](https://github.com/DanielPAagustinho/r2t_2025) on CMV and HepC datasets

## Step1 HepC and CMV

## Step2 HepC

```bash
sbatch v2t_step2_hepc_test3.slurm
```

Now let's examine the result of our 3-day run on the cluster through SLURM
Get the Number of log "metrics" file with the words 'virus2tree+   COMPLETED '
```bash
cd results_hepc_test2
srun bash -lc 'awk '\''($2 ~ /^virus2tre\+$/ || $2 ~ /^virus2tree\+$/) && $3=="COMPLETED"'\'' *.metrics | wc -l'
```
See which files don't have the word "COMPLETED"
```bash
srun \
>   bash -lc '
> missing=0
> for f in *.metrics; do
>   if ! awk '\''(($2 ~ /^virus2tre\+$/ || $2 ~ /^virus2tree\+$/) && $3=="COMPLETED") {found=1; exit} END { if (!found) exit 1 }'\'' "$f"; then
>     echo "$f"
>     missing=$((missing+1))
>   fi
> done
> echo "FILES WITHOUT COMPLETED: $missing"'
```
Redirect the SAMPLEs and JOB_IDs of those files to missing_labels.tsv to re-run them
```bash
srun \
  bash -lc '
grep -l -E "virus2tre\\+\\s+COMPLETED|virus2tree\\+\\s+COMPLETED" *.metrics | sort > /tmp/with.txt
ls *.metrics | sort > /tmp/all.txt
comm -23 /tmp/all.txt /tmp/with.txt > /tmp/missing.txt
echo -e "FILE\tLABEL\tSAMPLE\tJOB_ID"
while IFS= read -r f; do
  label=$(awk -F": " '\''/^LABEL[[:space:]]*:/ {print $2; exit}'\'' "$f")
  sample=$(awk -F": " '\''/^SAMPLE[[:space:]]*:/ {print $2; exit}'\'' "$f")
  jobid=$(awk -F": " '\''/^JOB_ID[[:space:]]*:/ {print $2; exit}'\'' "$f")
  echo -e "${f}\t${label}\t${sample}\t${jobid}"
done < /tmp/missing.txt
' > missing_labels.tsv
```
Now, to avoid problems with the r2t dir for those incompleted samples, we delete them from the r2t dir

```bash
# Assuming we are in results_hepc_test2 dir
# We count the number of IDs for 04 and 05 and then the number of missing ids present
cd r2t_ref
srun ls -d 04* | cut -d "_" -f 4-| sort -u | wc -l #779
srun ls -d 05* | cut -d "_" -f 5-|rev | cut -d "_" -f2-|rev|  sort -u | wc -l  #743
cd ..
srun bash -lc 'awk -F"\t" '\''NR>1{print $3}'\'' missing_labels.tsv | sort -u | xargs -I{} find r2t_ref -maxdepth 1 -type d \( -name "04_*{}*" -o -name "05_*{}*" \) | sort -u | wc -l' #36

#Deleting
srun bash -lc 'R2T="r2t_ref"; MISSING="missing_labels.tsv"; awk -F"\t" '\''NR>1{print $3}'\'' "$MISSING" | sort -u | while IFS= read -r s; do
  find "$R2T" -maxdepth 1 -type d \( -name "04_*${s}*" -o -name "05_*${s}*" \) -exec rm -rf {} +;
done'
# if u are unsure, u can verify what will be deleted replacing -exec rm -rf {} + with -print

#Now check the left folders
cd r2t_ref
srun ls -d 04* | cut -d "_" -f 4-| sort -u | wc -l #743 Here 36 were deleted
srun ls -d 05* | cut -d "_" -f 5-|rev | cut -d "_" -f2-|rev|  sort -u | wc -l #743
cd ..

```

After checking, better generate a more comprehensive file report with all the FAILED, COMPLETED and those that didn't finish: the NO FLAG samples.

```bash
srun \
  bash -lc '
shopt -s nullglob

OUT="final_report_hepc.tsv"
: > "$OUT"

echo -e "FILE\tLABEL\tSAMPLE\tJOB_ID\tSTATUS" | tee -a "$OUT"

for f in $(ls *.metrics 2>/dev/null | sort); do
  label=$(awk -F": " '\''/^LABEL[[:space:]]*:/ {print $2; exit}'\'' "$f")
  sample=$(awk -F": " '\''/^SAMPLE[[:space:]]*:/ {print $2; exit}'\'' "$f")
  jobid=$(awk -F": " '\''/^JOB_ID[[:space:]]*:/ {print $2; exit}'\'' "$f")

  status="NO_FLAG"
  if grep -Eq "virus2tre\\+\\s+FAILED|virus2tree\\+\\s+FAILED" "$f"; then
    status="FAILED"
  elif grep -Eq "virus2tre\\+\\s+COMPLETED|virus2tree\\+\\s+COMPLETED" "$f"; then
    status="COMPLETED"
  fi

  echo -e "${f}\t${label}\t${sample}\t${jobid}\t${status}" | tee -a "$OUT"
done

echo "[OK] Generated $OUT"
'
...
v2t_step2_945910_85.metrics     HepC_SRR5122831 SRR5122831      945998_85    COMPLETED
v2t_step2_945910_8.metrics      HepC_ERR245514  ERR245514       945919_8     FAILED
v2t_step2_945910_9.metrics      HepC_ERR245515  ERR245515       945920_9     FAILED
[OK] Generated final_report_hepc.tsv

```
We use the script  `v2t_step2_hepc_test3_recovery.slurm` for re-running not "COMPLETED" samples. In it we also re-ran de PacBio samples. We put as input to read2tree the specific minimap2 option `-ax map-pb` for CLR (no hi-fi) PacBio read sets. If you are familiar with read2tree, you know that required the following successful [pull request](https://github.com/DessimozLab/read2tree/pull/69) to polish a bit the --read_type argument. 

Now, first we also should deleted all PacBio samples from the read2tree folder `r2t_ref`. I have the 50 samples that were from PacBio in a `pacbio_samples` file. 
```bash
# Enter again the r2t folder
cd r2t_ref
grep -vE '^\s*(#|$)' .../pacbio_samples \
  | xargs -I{} rm -rf r2t_ref/04_*{}* r2t_ref/05_*{}*

```

Let's execute it with:

```bash
sbatch v2t_step2_hepc_test3_recovery.slurm
```

Luego de volver a generar el ``final_report_hepc.tsv` para que analice el resultado de la nuevas corridas, veamos las lineas con SAMPLEs duplicadas. Luego de contarlas, las imprimimos en un nuevo archivo, para verificar si la recorrida tuvo algún efecto, como cambio de NO FLAG a FAILED o COMPLETED.

```bash
cut -f3 final_report_hepc.tsv | sort | uniq -d | wc -l
## Here we redirect the rows were duplicates have been found to a new file 
cut -f3 final_report_hepc.tsv | sort | uniq -d | grep -F -f - final_report_hepc.tsv | sort -k3,3 -k1,1 > final_report_duplicated_samples_hepc.tsv
```

So, looking at final_report_duplicated_samples_hepc.tsv, we see that most samples have no changed flag as we would liked to!!, The changes are mainly from NO_FLAG to FAILED. Although a few indeed changed from NO_FLAG to COMPLETED.

I wanna now if all PacBio were completed successfully. So will mark in final_report_duplicated_samples_hepc.tsv with a new column that says PACBIO if the samples belongs to that technology. I'll use a simple awk script

```bash
awk 'NR==FNR{
  pb[$1]=1
  next
}
{ OFS="\t"
flag = ($3 in pb) ? "PACBIO" : "SHORT_READS"
print $0, flag                                  # Here we add the column that says if is PacBio or SR
}' pacbio_samples final_report_duplicated_samples_hepc.tsv > duplicates.tmp
mv duplicates.tmp final_report_duplicated_samples_hepc.tsv
```
Then we inspect the number of PACBIO ...
```bash
grep 'PACBIO' final_report_duplicated_samples_hepc.tsv | cut -f5 | sort | uniq -c
      95 COMPLETED
      5 FAILED
```
... and discover a few have failed!! We should run them again, but understand the problem with the SLURM jobs first!!

Let's check quickly again how many of the **unique** samples in final_report_hepc.tsv were completed.

```bash
cut -f3 final_report_hepc.tsv | sort | uniq -c | awk '$1>1' | rev | cut -d' ' -f1 | rev | grep -vF -f - final_report_hepc.tsv | tail -n +2 | cut -f5 | sort | uniq -c
    693 COMPLETED
```

Yep!! That is the exact number that should be there. So, now we find the samples in `final_report_duplicated_samples_hepc.tsv` that need to be re-ran. That is, the samples that haven't failed because no reads were mapped to the orthologs (these ones simply have no reads similar enough to the ref), but that were truncated in the previous SLURM runs.

```bash
#First let's see how many failed
grep -c 'FAILED' final_report_duplicated_samples_hepc.tsv
31
# And now let's see the reason they failed.
grep 'FAILED' final_report_duplicated_samples_hepc.tsv | cut -f1 | xargs -I{} grep -Hin "mapped_reads_species is empty" {} | wc -l
31
```
Therefore, all the FAILED read sets are due to "mapped_reads_species is empty". 

Next, let's check the COMPLETED IDs we have in our logs correspond to those already present in the 05_og_map files in r2t dir

```bash
grep '945910_.*COMPLETED' final_report_duplicated_samples_hepc.tsv | cut -f3 | wc -l
51 
#These are the new ones completed, we should sum them with the 693 previously completed samples
echo "51 + 693" | bc
744
# And now let's see which samples passed all and generated their completed fol 05_ogs_map folder inside the r2t dir
cd rt2_ref
ls -d 05* | cut -d "_" -f 5-|rev | cut -d "_" -f2-|rev|  sed 's/_1$//' |sort -u | wc -l
745
cd..
```
symbolo de "ojo" subnote:
grep '945910_.*FAILED' final_report_dup
licated_samples_hepc.tsv | cut -f3  | wc -l 
19
grep '945910_.*NO_FLAG' final_report_duplicated_samples_hepc.tsv | cut -f3 | wc -l
16 
#Then 
echo "51 + 693+19+16" | bc
779
And this is exactly the number of read sets for HepC!
end of subnote

Oh no!!, which is the single different sample that seems to be completed but is not being **labeled** as completed by the metric files logs?
```bash
# Get all first 693 sample IDs
cut -f3 final_report_hepc.tsv | sort | uniq -d | grep -vF -f - final_report_hepc.tsv \ #discard duplicated samples and count just the completed on the first round (the 693)
| grep 'COMPLETED' | cut -f3 > samples_completed_second_run.tsv
# Add also sample IDs "COMPLETED" from the second run (which we identify by the job id 945910)
grep '945910_' final_report_hepc.tsv | grep 'COMPLETED' | cut -f3 >> samples_completed_second_run.tsv
cd r2t
# Second, get the IDs that are actually completed
srun ls -d 05* | cut -d "_" -f 5-|rev | cut -d "_" -f2-|rev|  sed 's/_1$//' |sort -u > ../actual_completed_second_run.tsv
cd ..
# Find the diff
cat samples_completed_second_run.tsv actual_completed_second_run.tsv \ # Concatenate, so all will be duplicated except the one we like
  | sort \
  | uniq -u  #identify the single sample not duplicated
SRR1170678
```
And this is how the metrics file from SRR1170678 looks at the end 
![figure1](imgs/image.png)

So it has a curious situation. Read2tree had finished the mapping, which seems to correspond to an 05_ogs_map folder already generated, but then this job was cut by Slurm. We won't complicate and add them together with the samples that have NO_FLAG to re-run it again.

We now need to put all the samples that have NO_FLAG  in a file for their second run. The second run samples are from the JOB_ID "945910"

```bash
grep '945910_.*NO_FLAG' final_report_duplicated_samples_hepc.tsv | cut -f3 > no_flag_second_run.tsv
```
Again, we delete any rest of the r2t dirs associated to these samples
```bash
cut -f1 no_flag_second_run.tsv | sort -u \
>   | while read sample; do
>       find r2t_ref -maxdepth 1 -type d \
>         \( -name "04_*${sample}*" -o -name "05_*${sample}*" \) \
>         -exec rm -rf {} +
>     don
```

Let's run a second recovery script
```bash
nohup bash -lc '
cut -f1 ../results_hepc_test2/no_flag_second_run.tsv \
  | sort -u \
  | xargs -r -n1 -I{} sbatch --export=ALL,SAMPLE="{}" v2t_step2_hepc_test3_recovery2.slurm' > nohup_hepc_$(date +%F_%H%M).log 2>&1 &
```
