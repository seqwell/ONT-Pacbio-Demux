process NANOFILT {
  
  input:
    tuple val(pair_id),  path(fq) 

  output:
    tuple  val(pair_id), path('*.fastq.gz') 

     """
     gunzip -c $fq | NanoFilt -q 12   | gzip > ${pair_id}.filtered.fastq.gz
     """
  
  
}