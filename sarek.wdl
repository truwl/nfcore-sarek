version 1.0

workflow sarek {
	input{
		String step = "mapping"
		File samplesheet
		String outdir
		String? tools
		Boolean? no_intervals
		Float nucleotides_per_second = 1000
		Boolean? sentieon
		String? skip_tools
		File? target_bed
		Boolean? wes
		Boolean? trim_fastq
		Int? clip_r1
		Int? clip_r2
		Int? three_prime_clip_r1
		Int? three_prime_clip_r2
		Int? trim_nextseq
		Boolean? save_trimmed
		Int? split_fastq
		Boolean? save_split_fastqs
		String? umi_read_structure
		String group_by_umi_strategy = "Adjacency"
		String aligner = "bwa-mem"
		String? use_gatk_spark
		Boolean? save_bam_mapped
		Boolean? only_paired_variant_calling
		Float ploidy = 2
		Float? ascat_purity
		Float cf_coeff = 0.05
		Boolean? cf_contamination_adjustment
		Float? cf_contamination
		Float? cf_minqual
		Float? cf_mincov
		Float? cf_window
		Boolean? joint_germline
		Boolean? generate_gvcf
		Boolean? no_strelka_bp
		File? pon
		String? pon_tbi
		Boolean? ignore_soft_clipped_bases
		Boolean? vep_dbnsfp
		String? dbnsfp
		String? dbnsfp_tbi
		Boolean? vep_loftee
		Boolean? vep_spliceai
		String? spliceai_snv
		String? spliceai_snv_tbi
		String? spliceai_indel
		String? spliceai_indel_tbi
		Boolean? vep_spliceregion
		Boolean? annotation_cache
		String? snpeff_cache
		String? vep_cache
		String genome = "GATK.GRCh38"
		File? ac_loci
		File? ac_loci_gc
		File? bwa
		String? bwamem2
		File? chr_dir
		File? dbsnp
		String? dbsnp_tbi
		File? dict
		String? dragmap
		File? fasta
		File? fasta_fai
		File? germline_resource
		String? germline_resource_tbi
		File? intervals
		File? known_indels
		String? known_indels_tbi
		File? mappability
		String? snpeff_db
		String? snpeff_genome
		String? vep_genome
		String? vep_species
		Float? vep_cache_version
		Boolean? save_reference
		String igenomes_base = "s3://ngi-igenomes/igenomes/"
		Boolean? igenomes_ignore
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url
		String? seq_center
		String seq_platform = "ILLUMINA"
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		Boolean? help
		String publish_dir_mode = "copy"
		String? email
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_title
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		Boolean validate_params = true
		Boolean? show_hidden_params
		Boolean? enable_conda

	}

	call make_uuid as mkuuid {}
	call touch_uuid as thuuid {
		input:
			outputbucket = mkuuid.uuid
	}
	call run_nfcoretask as nfcoretask {
		input:
			step = step,
			samplesheet = samplesheet,
			outdir = outdir,
			tools = tools,
			no_intervals = no_intervals,
			nucleotides_per_second = nucleotides_per_second,
			sentieon = sentieon,
			skip_tools = skip_tools,
			target_bed = target_bed,
			wes = wes,
			trim_fastq = trim_fastq,
			clip_r1 = clip_r1,
			clip_r2 = clip_r2,
			three_prime_clip_r1 = three_prime_clip_r1,
			three_prime_clip_r2 = three_prime_clip_r2,
			trim_nextseq = trim_nextseq,
			save_trimmed = save_trimmed,
			split_fastq = split_fastq,
			save_split_fastqs = save_split_fastqs,
			umi_read_structure = umi_read_structure,
			group_by_umi_strategy = group_by_umi_strategy,
			aligner = aligner,
			use_gatk_spark = use_gatk_spark,
			save_bam_mapped = save_bam_mapped,
			only_paired_variant_calling = only_paired_variant_calling,
			ploidy = ploidy,
			ascat_purity = ascat_purity,
			cf_coeff = cf_coeff,
			cf_contamination_adjustment = cf_contamination_adjustment,
			cf_contamination = cf_contamination,
			cf_minqual = cf_minqual,
			cf_mincov = cf_mincov,
			cf_window = cf_window,
			joint_germline = joint_germline,
			generate_gvcf = generate_gvcf,
			no_strelka_bp = no_strelka_bp,
			pon = pon,
			pon_tbi = pon_tbi,
			ignore_soft_clipped_bases = ignore_soft_clipped_bases,
			vep_dbnsfp = vep_dbnsfp,
			dbnsfp = dbnsfp,
			dbnsfp_tbi = dbnsfp_tbi,
			vep_loftee = vep_loftee,
			vep_spliceai = vep_spliceai,
			spliceai_snv = spliceai_snv,
			spliceai_snv_tbi = spliceai_snv_tbi,
			spliceai_indel = spliceai_indel,
			spliceai_indel_tbi = spliceai_indel_tbi,
			vep_spliceregion = vep_spliceregion,
			annotation_cache = annotation_cache,
			snpeff_cache = snpeff_cache,
			vep_cache = vep_cache,
			genome = genome,
			ac_loci = ac_loci,
			ac_loci_gc = ac_loci_gc,
			bwa = bwa,
			bwamem2 = bwamem2,
			chr_dir = chr_dir,
			dbsnp = dbsnp,
			dbsnp_tbi = dbsnp_tbi,
			dict = dict,
			dragmap = dragmap,
			fasta = fasta,
			fasta_fai = fasta_fai,
			germline_resource = germline_resource,
			germline_resource_tbi = germline_resource_tbi,
			intervals = intervals,
			known_indels = known_indels,
			known_indels_tbi = known_indels_tbi,
			mappability = mappability,
			snpeff_db = snpeff_db,
			snpeff_genome = snpeff_genome,
			vep_genome = vep_genome,
			vep_species = vep_species,
			vep_cache_version = vep_cache_version,
			save_reference = save_reference,
			igenomes_base = igenomes_base,
			igenomes_ignore = igenomes_ignore,
			custom_config_version = custom_config_version,
			custom_config_base = custom_config_base,
			config_profile_name = config_profile_name,
			config_profile_description = config_profile_description,
			config_profile_contact = config_profile_contact,
			config_profile_url = config_profile_url,
			seq_center = seq_center,
			seq_platform = seq_platform,
			max_cpus = max_cpus,
			max_memory = max_memory,
			max_time = max_time,
			help = help,
			publish_dir_mode = publish_dir_mode,
			email = email,
			email_on_fail = email_on_fail,
			plaintext_email = plaintext_email,
			max_multiqc_email_size = max_multiqc_email_size,
			monochrome_logs = monochrome_logs,
			multiqc_title = multiqc_title,
			multiqc_config = multiqc_config,
			tracedir = tracedir,
			validate_params = validate_params,
			show_hidden_params = show_hidden_params,
			enable_conda = enable_conda,
			outputbucket = thuuid.touchedbucket
            }
		output {
			Array[File] results = nfcoretask.results
		}
	}
task make_uuid {
	meta {
		volatile: true
}

command <<<
        python <<CODE
        import uuid
        print("gs://truwl-internal-inputs/nf-sarek/{}".format(str(uuid.uuid4())))
        CODE
>>>

  output {
    String uuid = read_string(stdout())
  }
  
  runtime {
    docker: "python:3.8.12-buster"
  }
}

task touch_uuid {
    input {
        String outputbucket
    }

    command <<<
        echo "sentinel" > sentinelfile
        gsutil cp sentinelfile ~{outputbucket}/sentinelfile
    >>>

    output {
        String touchedbucket = outputbucket
    }

    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task fetch_results {
    input {
        String outputbucket
        File execution_trace
    }
    command <<<
        cat ~{execution_trace}
        echo ~{outputbucket}
        mkdir -p ./resultsdir
        gsutil cp -R ~{outputbucket} ./resultsdir
    >>>
    output {
        Array[File] results = glob("resultsdir/*")
    }
    runtime {
        docker: "google/cloud-sdk:latest"
    }
}

task run_nfcoretask {
    input {
        String outputbucket
		String step = "mapping"
		File samplesheet
		String outdir
		String? tools
		Boolean? no_intervals
		Float nucleotides_per_second = 1000
		Boolean? sentieon
		String? skip_tools
		File? target_bed
		Boolean? wes
		Boolean? trim_fastq
		Int? clip_r1
		Int? clip_r2
		Int? three_prime_clip_r1
		Int? three_prime_clip_r2
		Int? trim_nextseq
		Boolean? save_trimmed
		Int? split_fastq
		Boolean? save_split_fastqs
		String? umi_read_structure
		String group_by_umi_strategy = "Adjacency"
		String aligner = "bwa-mem"
		String? use_gatk_spark
		Boolean? save_bam_mapped
		Boolean? only_paired_variant_calling
		Float ploidy = 2
		Float? ascat_purity
		Float cf_coeff = 0.05
		Boolean? cf_contamination_adjustment
		Float? cf_contamination
		Float? cf_minqual
		Float? cf_mincov
		Float? cf_window
		Boolean? joint_germline
		Boolean? generate_gvcf
		Boolean? no_strelka_bp
		File? pon
		String? pon_tbi
		Boolean? ignore_soft_clipped_bases
		Boolean? vep_dbnsfp
		String? dbnsfp
		String? dbnsfp_tbi
		Boolean? vep_loftee
		Boolean? vep_spliceai
		String? spliceai_snv
		String? spliceai_snv_tbi
		String? spliceai_indel
		String? spliceai_indel_tbi
		Boolean? vep_spliceregion
		Boolean? annotation_cache
		String? snpeff_cache
		String? vep_cache
		String genome = "GATK.GRCh38"
		File? ac_loci
		File? ac_loci_gc
		File? bwa
		String? bwamem2
		File? chr_dir
		File? dbsnp
		String? dbsnp_tbi
		File? dict
		String? dragmap
		File? fasta
		File? fasta_fai
		File? germline_resource
		String? germline_resource_tbi
		File? intervals
		File? known_indels
		String? known_indels_tbi
		File? mappability
		String? snpeff_db
		String? snpeff_genome
		String? vep_genome
		String? vep_species
		Float? vep_cache_version
		Boolean? save_reference
		String igenomes_base = "s3://ngi-igenomes/igenomes/"
		Boolean? igenomes_ignore
		String custom_config_version = "master"
		String custom_config_base = "https://raw.githubusercontent.com/nf-core/configs/master"
		String? config_profile_name
		String? config_profile_description
		String? config_profile_contact
		String? config_profile_url
		String? seq_center
		String seq_platform = "ILLUMINA"
		Int max_cpus = 16
		String max_memory = "128.GB"
		String max_time = "240.h"
		Boolean? help
		String publish_dir_mode = "copy"
		String? email
		String? email_on_fail
		Boolean? plaintext_email
		String max_multiqc_email_size = "25.MB"
		Boolean? monochrome_logs
		String? multiqc_title
		String? multiqc_config
		String tracedir = "./results/pipeline_info"
		Boolean validate_params = true
		Boolean? show_hidden_params
		Boolean? enable_conda

	}
	command <<<
		export NXF_VER=21.10.5
		export NXF_MODE=google
		echo ~{outputbucket}
		/nextflow -c /truwl.nf.config run /sarek-5d0a0b97551a9700b9532ed4d32c4e9c543a998e  -profile truwl   --input ~{samplesheet} 	~{"--step '" + step + "'"} \ 
	~{"--samplesheet '" + samplesheet + "'"} \ 
	~{"--outdir '" + outdir + "'"} \ 
	~{"--tools '" + tools + "'"} \ 
	~{true="--no_intervals  " false="" no_intervals} \ 
	~{"--nucleotides_per_second " + nucleotides_per_second} \ 
	~{true="--sentieon  " false="" sentieon} \ 
	~{"--skip_tools '" + skip_tools + "'"} \ 
	~{"--target_bed '" + target_bed + "'"} \ 
	~{true="--wes  " false="" wes} \ 
	~{true="--trim_fastq  " false="" trim_fastq} \ 
	~{"--clip_r1 " + clip_r1} \ 
	~{"--clip_r2 " + clip_r2} \ 
	~{"--three_prime_clip_r1 " + three_prime_clip_r1} \ 
	~{"--three_prime_clip_r2 " + three_prime_clip_r2} \ 
	~{"--trim_nextseq " + trim_nextseq} \ 
	~{true="--save_trimmed  " false="" save_trimmed} \ 
	~{"--split_fastq " + split_fastq} \ 
	~{true="--save_split_fastqs  " false="" save_split_fastqs} \ 
	~{"--umi_read_structure '" + umi_read_structure + "'"} \ 
	~{"--group_by_umi_strategy '" + group_by_umi_strategy + "'"} \ 
	~{"--aligner '" + aligner + "'"} \ 
	~{"--use_gatk_spark '" + use_gatk_spark + "'"} \ 
	~{true="--save_bam_mapped  " false="" save_bam_mapped} \ 
	~{true="--only_paired_variant_calling  " false="" only_paired_variant_calling} \ 
	~{"--ploidy " + ploidy} \ 
	~{"--ascat_purity " + ascat_purity} \ 
	~{"--cf_coeff " + cf_coeff} \ 
	~{true="--cf_contamination_adjustment  " false="" cf_contamination_adjustment} \ 
	~{"--cf_contamination " + cf_contamination} \ 
	~{"--cf_minqual " + cf_minqual} \ 
	~{"--cf_mincov " + cf_mincov} \ 
	~{"--cf_window " + cf_window} \ 
	~{true="--joint_germline  " false="" joint_germline} \ 
	~{true="--generate_gvcf  " false="" generate_gvcf} \ 
	~{true="--no_strelka_bp  " false="" no_strelka_bp} \ 
	~{"--pon '" + pon + "'"} \ 
	~{"--pon_tbi '" + pon_tbi + "'"} \ 
	~{true="--ignore_soft_clipped_bases  " false="" ignore_soft_clipped_bases} \ 
	~{true="--vep_dbnsfp  " false="" vep_dbnsfp} \ 
	~{"--dbnsfp '" + dbnsfp + "'"} \ 
	~{"--dbnsfp_tbi '" + dbnsfp_tbi + "'"} \ 
	~{true="--vep_loftee  " false="" vep_loftee} \ 
	~{true="--vep_spliceai  " false="" vep_spliceai} \ 
	~{"--spliceai_snv '" + spliceai_snv + "'"} \ 
	~{"--spliceai_snv_tbi '" + spliceai_snv_tbi + "'"} \ 
	~{"--spliceai_indel '" + spliceai_indel + "'"} \ 
	~{"--spliceai_indel_tbi '" + spliceai_indel_tbi + "'"} \ 
	~{true="--vep_spliceregion  " false="" vep_spliceregion} \ 
	~{true="--annotation_cache  " false="" annotation_cache} \ 
	~{"--snpeff_cache '" + snpeff_cache + "'"} \ 
	~{"--vep_cache '" + vep_cache + "'"} \ 
	~{"--genome '" + genome + "'"} \ 
	~{"--ac_loci '" + ac_loci + "'"} \ 
	~{"--ac_loci_gc '" + ac_loci_gc + "'"} \ 
	~{"--bwa '" + bwa + "'"} \ 
	~{"--bwamem2 '" + bwamem2 + "'"} \ 
	~{"--chr_dir '" + chr_dir + "'"} \ 
	~{"--dbsnp '" + dbsnp + "'"} \ 
	~{"--dbsnp_tbi '" + dbsnp_tbi + "'"} \ 
	~{"--dict '" + dict + "'"} \ 
	~{"--dragmap '" + dragmap + "'"} \ 
	~{"--fasta '" + fasta + "'"} \ 
	~{"--fasta_fai '" + fasta_fai + "'"} \ 
	~{"--germline_resource '" + germline_resource + "'"} \ 
	~{"--germline_resource_tbi '" + germline_resource_tbi + "'"} \ 
	~{"--intervals '" + intervals + "'"} \ 
	~{"--known_indels '" + known_indels + "'"} \ 
	~{"--known_indels_tbi '" + known_indels_tbi + "'"} \ 
	~{"--mappability '" + mappability + "'"} \ 
	~{"--snpeff_db '" + snpeff_db + "'"} \ 
	~{"--snpeff_genome '" + snpeff_genome + "'"} \ 
	~{"--vep_genome '" + vep_genome + "'"} \ 
	~{"--vep_species '" + vep_species + "'"} \ 
	~{"--vep_cache_version " + vep_cache_version} \ 
	~{true="--save_reference  " false="" save_reference} \ 
	~{"--igenomes_base '" + igenomes_base + "'"} \ 
	~{true="--igenomes_ignore  " false="" igenomes_ignore} \ 
	~{"--custom_config_version '" + custom_config_version + "'"} \ 
	~{"--custom_config_base '" + custom_config_base + "'"} \ 
	~{"--config_profile_name '" + config_profile_name + "'"} \ 
	~{"--config_profile_description '" + config_profile_description + "'"} \ 
	~{"--config_profile_contact '" + config_profile_contact + "'"} \ 
	~{"--config_profile_url '" + config_profile_url + "'"} \ 
	~{"--seq_center '" + seq_center + "'"} \ 
	~{"--seq_platform '" + seq_platform + "'"} \ 
	~{"--max_cpus " + max_cpus} \ 
	~{"--max_memory '" + max_memory + "'"} \ 
	~{"--max_time '" + max_time + "'"} \ 
	~{true="--help  " false="" help} \ 
	~{"--publish_dir_mode '" + publish_dir_mode + "'"} \ 
	~{"--email '" + email + "'"} \ 
	~{"--email_on_fail '" + email_on_fail + "'"} \ 
	~{true="--plaintext_email  " false="" plaintext_email} \ 
	~{true="--monochrome_logs  " false="" monochrome_logs} \ 
	~{"--multiqc_title '" + multiqc_title + "'"} \ 
	~{"--multiqc_config '" + multiqc_config + "'"} \ 
	~{"--tracedir '" + tracedir + "'"} \ 
	~{true="--validate_params  " false="" validate_params} \ 
	~{true="--show_hidden_params  " false="" show_hidden_params} \ 
	~{true="--enable_conda  " false="" enable_conda} \ 
	-w ~{outputbucket}
	>>>
        
    output {
        File execution_trace = "pipeline_execution_trace.txt"
        Array[File] results = glob("results/*/*html")
    }
    runtime {
        docker: "truwl/nfcore-sarek:5d0a0b97551a9700b9532ed4d32c4e9c543a998e_0.1.0"
        memory: "2 GB"
        cpu: 1
    }
}
    