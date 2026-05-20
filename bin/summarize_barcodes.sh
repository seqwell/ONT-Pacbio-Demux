#!/bin/bash
set -euo pipefail
set -x

pool_id=${1:?Missing pool_id}
error_rate=${2:?Missing error_rate}
summary_file="${pool_id}_demux_report.csv"

declare -A barcode_counts
total_reads=0
unknown_reads=0
unclassified_me_tagged_reads=0
me_tagged_reads=0
other_barcodes_reads=0

for file in *.fastq.gz; do
    [[ -f "$file" ]] || continue
    barcode="${file%.fastq.gz}"
    reads=$(zcat "$file" 2>/dev/null | awk 'END{print int(NR/4)}' || echo 0)
    reads="${reads//[$'\t\r\n ']}"

    if [[ "$barcode" == "unknown"* ]] || [[ "$barcode" == "unclassified"* ]]; then
        if [[ "$barcode" == "unclassified.ME.tagged" ]]; then
            unclassified_me_tagged_reads=$((unclassified_me_tagged_reads + reads))
        else
            unknown_reads=$((unknown_reads + reads))
        fi
    elif [[ "$barcode" == *".ME.tagged" ]]; then
        me_tagged_reads=$((me_tagged_reads + reads))
    else
        barcode_counts["$barcode"]=$reads
    fi
done

for count in "${barcode_counts[@]}"; do
    total_reads=$((total_reads + count))
done
total_reads=$((total_reads + unknown_reads + unclassified_me_tagged_reads + me_tagged_reads))

num_barcodes=${#barcode_counts[@]}
avg_reads=0
if (( num_barcodes > 0 )); then
    sum_reads=0
    for r in "${barcode_counts[@]}"; do sum_reads=$((sum_reads + r)); done
    avg_reads=$((sum_reads / num_barcodes))
fi

threshold=$((avg_reads / 10))

declare -A final_barcode_counts
for barcode in "${!barcode_counts[@]}"; do
    count=${barcode_counts[$barcode]}
    if (( count < threshold )); then
        other_barcodes_reads=$((other_barcodes_reads + count))
    else
        final_barcode_counts["$barcode"]=$count
    fi
done

regular_barcoded_reads=0
for count in "${final_barcode_counts[@]}"; do
    regular_barcoded_reads=$((regular_barcoded_reads + count))
done
regular_barcodes_sum=$regular_barcoded_reads

unknown_pct=$(awk -v u="$unknown_reads" -v t="$total_reads" 'BEGIN{printf "%.2f", (t>0)?u*100/t:0}')
unclassified_me_pct=$(awk -v m="$unclassified_me_tagged_reads" -v t="$total_reads" 'BEGIN{printf "%.2f", (t>0)?m*100/t:0}')
regular_barcoded_pct=$(awk -v r="$regular_barcoded_reads" -v t="$total_reads" 'BEGIN{printf "%.2f", (t>0)?r*100/t:0}')
me_tagged_pct=$(awk -v m="$me_tagged_reads" -v t="$total_reads" 'BEGIN{printf "%.2f", (t>0)?m*100/t:0}')
other_barcodes_pct=$(awk -v o="$other_barcodes_reads" -v t="$total_reads" 'BEGIN{printf "%.2f", (t>0)?o*100/t:0}')

{
    echo "Total Reads:,${total_reads}"
    echo "Total Demuxed Reads:,${regular_barcoded_reads},${regular_barcoded_pct}%"
    echo "ME-tagged Reads (demuxed reads: removed):,${me_tagged_reads},${me_tagged_pct}%"
    echo "Unclassified Reads No ME-tagged:,${unknown_reads},${unknown_pct}%"
    echo "Unclassified ME-tagged Reads:,${unclassified_me_tagged_reads},${unclassified_me_pct}%"
    echo "Other Barcodes (count < 10% of avg):,${other_barcodes_reads},${other_barcodes_pct}%"
    echo ""
    echo ""
    echo "Barcode,Read_Count,Percent_of_Total"
    for barcode in $(printf '%s\n' "${!final_barcode_counts[@]}" | sort -V); do
        count=${final_barcode_counts[$barcode]}
        pct=$(awk -v c="$count" -v t="$regular_barcodes_sum" 'BEGIN{printf "%.2f", (t>0)?c*100/t:0}')
        echo "${barcode},${count},${pct}"
    done
} > "$summary_file"

echo "Successfully generated: $summary_file"
