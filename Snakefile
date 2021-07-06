TMPDIR = "/localscratch/msauria"
SYSTEM="linux.x86_64"
MAXMEM = "1024"
MAXTHREADS = 100
THRESHOLD = 100000
KMERS = ["75"]
HOR_KMERS = ["100"]
BSAT_KMERS = ["21"]
ANNOTATION_FILE="data/t2t_cenAnnotation.v2.021921.FinalColors.bed"

ARRAYS = set()
CHROMS = set()
REGIONS = []
ALL_REGIONS = []
HOR_REGIONS = []
BSAT_REGIONS = []
for line in open(ANNOTATION_FILE):
    line = line.split()
    CHROMS.add(line[0])
    ALL_REGIONS.append(line[3].split("(")[0])
    if line[3].count("hor") > 0:
        HOR_REGIONS.append(line[3].split("(")[0])
    if line[3].count("bsat") > 0:
        BSAT_REGIONS.append(line[3].split("(")[0])
    if (not line[3].count('arm') > 0
        and not line[3].startswith('gap')
        and not line[3].startswith('rDNA')
        and int(line[2]) - int(line[1]) >= THRESHOLD):
        REGIONS.append(line[3].split("(")[0])
        ARRAYS.add(line[3].split('_')[0])
ARRAYS = list(ARRAYS)
ARRAYS.sort()
NON_CT_REGIONS = list(REGIONS)
CT_REGIONS = []
for i in range(len(NON_CT_REGIONS))[::-1]:
    if NON_CT_REGIONS[i].count('ct') > 0:
        CT_REGIONS.append(NON_CT_REGIONS.pop(i))
CHROMS = list(CHROMS)
CHROMS.sort()
PAIR_DICT = {}
for i in range(len(NON_CT_REGIONS)):
    PAIR_DICT[NON_CT_REGIONS[i]] = []
    for j in range(i + 1, len(NON_CT_REGIONS)):
        PAIR_DICT[NON_CT_REGIONS[i]].append(NON_CT_REGIONS[j])

CT_PAIR_DICT = {}
for i in range(len(CT_REGIONS)):
    CT_PAIR_DICT[CT_REGIONS[i]] = []
    for j in range(i + 1, len(CT_REGIONS)):
        CT_PAIR_DICT[CT_REGIONS[i]].append(CT_REGIONS[j])

HOR_PAIR_DICT = {}
for i in range(len(HOR_REGIONS)):
    HOR_PAIR_DICT[HOR_REGIONS[i]] = []
    for j in range(i + 1, len(HOR_REGIONS)):
        HOR_PAIR_DICT[HOR_REGIONS[i]].append(HOR_REGIONS[j])

BSAT_PAIR_DICT = {}
for i in range(len(BSAT_REGIONS)):
    BSAT_PAIR_DICT[BSAT_REGIONS[i]] = []
    for j in range(i + 1, len(BSAT_REGIONS)):
        BSAT_PAIR_DICT[BSAT_REGIONS[i]].append(BSAT_REGIONS[j])
ALL_KMERS = list(set(KMERS + HOR_KMERS + BSAT_KMERS))

ALL_REGION_SET = "|".join(ALL_REGIONS)
KMER_SET = "|".join(KMERS)
HOR_KMER_SET = "|".join(HOR_KMERS)
BSAT_KMER_SET = "|".join(BSAT_KMERS)
ALL_KMER_SET = "|".join(ALL_KMERS)
ARRAY_SET = "|".join(ARRAYS)
NON_CT_REGION_SET = "|".join(NON_CT_REGIONS)
CT_REGION_SET = "|".join(CT_REGIONS)
CHROM_SET = "|".join(CHROMS)
REGION_SET = "|".join(REGIONS)
PAIR_SET = {}
for key, value in PAIR_DICT.items():
    PAIR_SET[key] = "|".join(value)
CT_PAIR_SET = {}
for key, value in CT_PAIR_DICT.items():
    CT_PAIR_SET[key] = "|".join(value)
HOR_PAIR_SET = {}
for key, value in HOR_PAIR_DICT.items():
    HOR_PAIR_SET[key] = "|".join(value)
HOR_REGION_SET = "|".join(HOR_REGIONS)
BSAT_PAIR_SET = {}
for key, value in BSAT_PAIR_DICT.items():
    BSAT_PAIR_SET[key] = "|".join(value)
BSAT_REGION_SET = "|".join(BSAT_REGIONS)


rule all:
    input:
        expand("plots/circos_{array}_{kmer}.svg", kmer=KMERS,
               array=['all', 'mismatch', 'hor', 'hsat',
                      'censat', 'mon', 'bsat', 'ct']),
        "results/chm13v1_100mer_coverage.bed"
        #expand("results/chm13v1_HOR_{kmer}mers.txt", kmer=HOR_KMERS),
        #expand("results/chm13v1_BSAT_{kmer}mers.txt", kmer=BSAT_KMERS)


######## Get software #########

rule build_KMC:
    output:
        "bin/kmc",
        "bin/kmc_genome_counts",
        "bin/kmc_tools",
        "bin/kmc_info",
        "bin/kmc_dump"
    log:
        "logs/build_kmc.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        git clone https://github.com/msauria/KMC.git
        git checkout kmer_mapping
        cd KMC && make
        cp KMC/bin/kmc bin/
        chmod a+rx bin/kmc
        cp KMC/bin/kmc_genome_counts bin/
        chmod a+rx bin/kmc_genome_counts
        cp KMC/bin/kmc_tools bin/
        chmod a+rx bin/kmc_tools
        cp KMC/bin/kmc_info bin/
        chmod a+rx bin/kmc_info
        cp KMC/bin/kmc_dump bin/
        chmod a+rx bin/kmc_dump
        """

rule download_software:
    output:
        "bin/{sw}"
    wildcard_constraints:
        sw="wigToBigWig"
    params:
        sw="{sw}",
        system=SYSTEM
    log:
        "logs/Encode/download/{sw}.log"
    shell:
        """
        wget http://hgdownload.soe.ucsc.edu/admin/exe/{params.system}/{params.sw} -O {output}
        chmod a+rx {output}
        """



######## Get genomic data #########

rule get_chm13v1_fasta:
    output:
        fa="fasta/chm13v1.fasta",
        sizes="fasta/chm13v1.chrom.sizes"
    log:
        "logs/fasta/chm13v1_download.log"
    shell:
        """
        wget -O {output.fa}.gz http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/genome/t2t-chm13-v1.0.fa.gz
        gzip -d -c {output.fa}.gz > {output.fa}
        cat {output.fa} | awk 'BEGIN{{ OFS=""; ORS=""; tab="\t"; nl="\n"; tl=0 }}\
                               {{ if ( NR == 1 ){{ split($1,A,">"); print A[2],tab; }} \
                                  else {{ if ( $1 ~ /^>/ ){{ split($1,A,">"); print tl,nl,A[2],tab; tl = 0;}} \
                                          else {{ tl = tl + length; }} }} }}\
                               END{{ print tl,nl }}' | sort -r -k2,2n > {output.sizes}
        """

rule index_fasta:
    input:
        "fasta/{genome}.fasta"
    output:
        "fasta/{genome}.fasta.fai"
    params:
        ""
    log:
        "logs/fasta/{genome}_index.log"
    wrapper:
        "0.72.0/bio/samtools/faidx"

rule get_repeatmask:
    output:
        "data/t2t_chm13v1_repeatmask.bb"
    log:
        "logs/repeatmask.log"
    shell:
        """
        wget -O {output} http://t2t.gi.ucsc.edu/chm13/hub/t2t-chm13-v1.0/rmskV2/rmskV2.bigBed
        """

rule phony:
    output:
        "data/t2t_cenAnnotation.v2.021921.FinalColors.bed"
    shell:
        """
        touch {output}
        """


######## Fasta manipulations #########

rule preprocess_annotations:
    output:
        "data/t2t_cenSatAnnotation.bed"
    params:
        anno=ANNOTATION_FILE
    log:
        "logs/preprocess/cenSatAnnotation.log"
    shell:
        """
        grep -v "arm" {params.anno} | \
        awk 'BEGIN{{ OFS="\\t"; pchr=""; ps=-1; pe=-1; pn=""; pa=""; pcol="" }} \
             {{ split($4,a,"_"); \
                if(pa == "ct") {{ \
                    if(pe > $2 && pchr == $1) pe = $2; \
                    print pchr,ps,pe,pn,pa,pcol; }} \
                if (a[1] == "ct"){{ \
                    if($2 < pe) $2 = pe }} \
                else {{ print $1,$2,$3,$4,a[1],$9; }} \
                pchr=$1; ps=$2; pe=$3; pn=$4; pa=a[1]; pcol=$9; }} \
             END{{ if(pa == "ct") print pchr,ps,pe,pn,pcol; }}' > {output}
        """

rule extract_array_regions:
    input:
        fa="fasta/chm13v1.fasta",
        fai="fasta/chm13v1.fasta.fai",
        anno="data/t2t_cenSatAnnotation.bed",
        rptmask="data/t2t_chm13v1_repeatmask.bb"
    output:
        "fasta/regions/chm13v1_{region}.fasta"
    wildcard_constraints:
        region=ALL_REGION_SET
    params:
        region="{region}"
    conda:
        "envs/deepsam.yaml"
    shell:
        """
        region=$(eval grep -w "{params.region}" {input.anno} | awk '{{ printf "%s:%s-%s", $1, $2, $3 }}')
        bin/extract_masked_region.py {input.fa} {input.rptmask} $region > {output}
        """

rule extract_arrays:
    input:
        fa="fasta/chm13v1.fasta",
        fai="fasta/chm13v1.fasta.fai",
        anno="data/t2t_cenSatAnnotation.bed",
    output:
        "fasta/regions/chm13v1_{array}.fasta"
    wildcard_constraints:
        array=ARRAY_SET
    params:
        array="{array}"
    conda:
        "envs/deepsam.yaml"
    shell:
        """
        rm -f tmp.{params.array}.regions
        regions=($(eval grep -w {params.array} {input.anno} | awk '{{ printf "%s:%s-%s ", $1, $2, $3 }}'))
        for R in ${{regions[*]}}; do
            echo $R >> tmp.{params.array}.regions
        done
        samtools faidx {input.fa} -r tmp.{params.array}.regions > {output}
        rm -f tmp.{params.array}.regions
        """

rule extract_fasta_noncen_regions:
    input:
        fa="fasta/chm13v1.fasta",
        fai="fasta/chm13v1.fasta.fai",
        sizes="fasta/chm13v1.chrom.sizes",
        anno="data/t2t_cenSatAnnotation.bed"
    output:
        "fasta/regions/chm13v1_noncen.fasta"
    params:
        chroms=CHROMS
    conda:
        "envs/deepsam.yaml"
    shell:
        """
        rm -f region.tmp {output}.tmp
        for CHROM in {params.chroms}; do
            grep -w $CHROM {input.anno} > region.tmp
            SIZE=$(eval cat {input.sizes} | grep -w "$CHROM" | awk '{{ print $2 }}')
            cat region.tmp | \
                awk -v size="$SIZE" '{{ if ( NR == 1 ) {{ if ( $2 > 0 ) {{ printf "%s:0-%s\\n", $1, $2 }} end=$3; chrom=$1 }} \
                                        else {{ printf "%s:%s-%s\\n", $1, end, $2; end=$3 }} }} \
                                     END{{ if ( end < size ) {{ printf "%s:%s-%s\\n", chrom, end, size }} }}' >> {output}.tmp
        done
        samtools faidx {input.fa} -r {output}.tmp > {output}
        rm -f region.tmp {output}.tmp
        """

rule extract_fasta_centro_regions:
    input:
        fa="fasta/chm13v1.fasta",
        fai="fasta/chm13v1.fasta.fai",
        anno="data/t2t_cenSatAnnotation.bed"
    output:
        "fasta/regions/chm13v1_centro.fasta"
    conda:
        "envs/deepsam.yaml"
    shell:
        """
        rm -f tmp.centro.regions
        regions=($(eval cat {input.anno} | awk '{{ printf "%s:%s-%s ", $1, $2, $3 }}'))
        for R in ${{regions[*]}}; do
            echo $R >> tmp.centro.regions
        done
        samtools faidx {input.fa} -r tmp.centro.regions > {output}
        rm -f tmp.centro.regions
        """

######## Create KMC DBs #########

rule create_full_kmc_db:
    input:
        fa="fasta/chm13v1.fasta",
        kmc="bin/kmc"
    output:
        "KMC_db/chm13v1_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{kmer}.kmc_suf"
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="KMC_db/chm13v1_{kmer}"
    threads:
        MAXTHREADS
    wildcard_constraints:
        kmer=ALL_KMER_SET
    log:
        "logs/kmc/chm13v1_{kmer}.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc} -k{params.kmer} -m{params.maxmem} -fm -ci2 -t{threads} {input.fa} {params.prefix} {params.tmpdir}
        """

rule create_region_kmc_db:
    input:
        fa="fasta/regions/chm13v1_{region}.fasta",
        kmc="bin/kmc"
    output:
        "KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{region}_{kmer}.kmc_suf"
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="KMC_db/chm13v1_{region}_{kmer}"
    threads:
        MAXTHREADS
    wildcard_constraints:
        kmer=KMER_SET,
        region=REGION_SET
    log:
        "logs/kmc/chm13v1_{region}_{kmer}.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc} -k{params.kmer} -m{params.maxmem} -fm -ci1 -t{threads} {input.fa} {params.prefix} {params.tmpdir}
        """

rule create_array_kmc_db:
    input:
        fa="fasta/regions/chm13v1_{array}.fasta",
        kmc="bin/kmc"
    output:
        "KMC_db/chm13v1_{array}_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{array}_{kmer}.kmc_suf"
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="KMC_db/chm13v1_{array}_{kmer}"
    threads:
        MAXTHREADS
    wildcard_constraints:
        kmer=KMER_SET,
        array=ARRAY_SET
    log:
        "logs/kmc/chm13v1_{array}_{kmer}.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc} -k{params.kmer} -m{params.maxmem} -fm -ci1 -t{threads} {input.fa} {params.prefix} {params.tmpdir}
        """

rule create_array_pair_kmc_db:
    input:
        fa1="fasta/regions/chm13v1_{array1}.fasta",
        fa2="fasta/regions/chm13v1_{array2}.fasta",
        kmc="bin/kmc",
        array_db="KMC_db/chm13v1_{array1}_{kmer}.kmc_pre"
    output:
        "KMC_db/chm13v1_{array1}-{array2}_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{array1}-{array2}_{kmer}.kmc_suf",
        "KMC_db/chm13v1_{array2}-{array1}_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{array2}-{array1}_{kmer}.kmc_suf",
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="KMC_db/chm13v1_{array1}_{kmer}",
        prefix1="KMC_db/chm13v1_{array1}-{array2}_{kmer}",
        prefix2="KMC_db/chm13v1_{array2}-{array1}_{kmer}",
    wildcard_constraints:
        kmer=KMER_SET,
        array1=ARRAY_SET,
        array2=ARRAY_SET
    threads:
        MAXTHREADS
    log:
        "logs/kmc/chm13v1_{array1}_{array2}_{kmer}.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        rm -f {output}
        if [[ {params.prefix1} == {params.prefix2} ]]; then
            ln -s $PWD/{params.prefix}.kmc_pre {params.prefix1}.kmc_pre
            ln -s $PWD/{params.prefix}.kmc_suf {params.prefix1}.kmc_suf
        else
            cat {input.fa1} {input.fa2} > {params.prefix1}.fa
            {input.kmc} -k{params.kmer} -m{params.maxmem} -fm -ci1 -t{threads} {params.prefix1}.fa {params.prefix1} {params.tmpdir}
            ln -s $PWD/{params.prefix1}.kmc_pre {params.prefix2}.kmc_pre
            ln -s $PWD/{params.prefix1}.kmc_suf {params.prefix2}.kmc_suf
            rm {params.prefix1}.fa
        fi
        """

rule create_array_pair_specific_kmc_db:
    input:
        kmc_tools="bin/kmc_tools",
        full_db="KMC_db/chm13v1_{kmer}.kmc_pre",
        array_db="KMC_db/chm13v1_{array1}-{array2}_{kmer}.kmc_pre"
    output:
        "KMC_db/chm13v1_{array1}-{array2}_unique_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{array1}-{array2}_unique_{kmer}.kmc_suf",
        "KMC_db/chm13v1_{array2}-{array1}_unique_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{array2}-{array1}_unique_{kmer}.kmc_suf",
    params:
        full_prefix="KMC_db/chm13v1_{kmer}",
        prefix="KMC_db/chm13v1_{array1}-{array2}_{kmer}",
        prefix1="KMC_db/chm13v1_{array1}-{array2}_unique_{kmer}",
        prefix2="KMC_db/chm13v1_{array2}-{array1}_unique_{kmer}"
    wildcard_constraints:
        kmer=KMER_SET,
        array1=ARRAY_SET,
        array2=ARRAY_SET
    log:
        "logs/kmc/chm13v1_{array1}_{array2}_unique_{kmer}.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        rm -f {output}
        {input.kmc_tools} simple {params.prefix} {params.full_prefix} counters_compare {params.prefix1}
        if [[ {params.prefix1} != {params.prefix2} ]]; then
            ln -s $PWD/{params.prefix1}.kmc_pre {params.prefix2}.kmc_pre
            ln -s $PWD/{params.prefix1}.kmc_suf {params.prefix2}.kmc_suf
        fi
        """

rule create_noncen_kmc_db:
    input:
        fa="fasta/regions/chm13v1_{region}.fasta",
        kmc="bin/kmc"
    output:
        "KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{region}_{kmer}.kmc_suf"
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="KMC_db/chm13v1_{region}_{kmer}"
    threads:
        MAXTHREADS
    wildcard_constraints:
        kmer=KMER_SET,
        region="noncen|centro"
    log:
        "logs/kmc/chm13v1_{region}_{kmer}.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc} -k{params.kmer} -m{params.maxmem} -fm -ci1 -t{threads} {input.fa} {params.prefix} {params.tmpdir}
        """

rule create_noncen_specific_kmc_db:
    input:
        fa="fasta/regions/chm13v1_{region}.fasta",
        kmc_tools="bin/kmc_tools",
        full_db="KMC_db/chm13v1_{kmer}.kmc_pre",
        region_db="KMC_db/chm13v1_{region}_{kmer}.kmc_pre"
    output:
        "KMC_db/chm13v1_{region}_unique_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{region}_unique_{kmer}.kmc_suf"
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="KMC_db/chm13v1_{region}_unique_{kmer}",
        full_prefix="KMC_db/chm13v1_{kmer}",
        region_prefix="KMC_db/chm13v1_{region}_{kmer}"
    wildcard_constraints:
        kmer=KMER_SET,
        region="noncen|centro"
    log:
        "logs/kmc/chm13v1_{region}_unique_{kmer}.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.region_prefix} {params.full_prefix} counters_compare {params.prefix}
        """

rule create_region_unique_kmc_db:
    input:
        fa="fasta/regions/chm13v1_{region}.fasta",
        kmc_tools="bin/kmc_tools",
        full_db="KMC_db/chm13v1_{kmer}.kmc_pre",
        region_db="KMC_db/chm13v1_{region}_{kmer}.kmc_pre"
    output:
        "KMC_db/chm13v1_{region}_unique_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{region}_unique_{kmer}.kmc_suf"
    params:
        prefix="KMC_db/chm13v1_{region}_unique_{kmer}",
        region_prefix="KMC_db/chm13v1_{region}_{kmer}",
        full_prefix="KMC_db/chm13v1_{kmer}"
    wildcard_constraints:
        kmer=KMER_SET,
        region=REGION_SET
    log:
        "logs/kmc/chm13v1_{region}_{kmer}.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.region_prefix} {params.full_prefix} counters_compare {params.prefix}
        """

rule count_db_kmers:
    input:
        db="KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
        kmc_info="bin/kmc_info"
    output:
        "KMC_db/chm13v1_{region}_{kmer}.count"
    params:
        prefix="KMC_db/chm13v1_{region}_{kmer}"
    wildcard_constraints:
        kmer=KMER_SET,
        region=REGION_SET
    log:
        "logs/kmc/chm13v1_{region}_{kmer}_count.log"
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_info} {params.prefix} | tail -n 1 | awk '{{ print $3 }}' > {output}
        """


######## Create Kmer BigWigs #########

rule create_kmer_wiggle:
    input:
        fa="fasta/chm13v1.fasta",
        db="KMC_db/chm13v1_{kmer}.kmc_pre",
        kmc="bin/kmc_genome_counts"
    output:
        temp("Kmer_BigWigs/chm13v1_{kmer}mer.wig")
    params:
        prefix="KMC_db/chm13v1_{kmer}",
        chroms=CHROMS,
        tmpdir=TMPDIR
    wildcard_constraints:
        kmer=ALL_KMER_SET
    threads:
        len(CHROMS)
    conda:
        "envs/kmc.yaml"
    shell:
        """
        tmpdir=$(mktemp -d -p {params.tmpdir});
        for C in {params.chroms}; do \
            {input.kmc} -ch$C {params.prefix} {input.fa} > $tmpdir/$C.wig & \
        done; wait
        wigs=""
        for C in {params.chroms}; do wigs="$wigs $tmpdir/$C.wig"; done;
        cat $wigs > {output}
        rm -r $tmpdir
        """

rule create_kmer_bigwig:
    input:
        wig="Kmer_BigWigs/chm13v1_{kmer}mer.wig",
        wigtobigwig="bin/wigToBigWig",
        sizes="fasta/chm13v1.chrom.sizes"
    output:
        "Kmer_BigWigs/chm13v1_{kmer}mer.bw"
    wildcard_constraints:
        kmer=ALL_KMER_SET
    shell:
        """
        {input.wigtobigwig} -clip {input.wig} {input.sizes} {output}
        """


######## Intersect KMC DBs #########

rule intersect_db_pairs:
    input:
        r1="KMC_db/chm13v1_{region1}_{kmer}.kmc_pre",
        r2="KMC_db/chm13v1_{region2}_{kmer}.kmc_pre",
        kmc_tools="bin/kmc_tools",
        kmc_info="bin/kmc_info"
    output:
        "intersections_{kmer}/chm13v1_{region1}/{region2}.txt"
    params:
        r1_prefix="KMC_db/chm13v1_{region1}_{kmer}",
        r2_prefix="KMC_db/chm13v1_{region2}_{kmer}",
        prefix="{region1}-{region2}_{kmer}",
        tmpdir=TMPDIR
    wildcard_constraints:
        kmer=KMER_SET,
        region1=NON_CT_REGION_SET,
        region2=NON_CT_REGION_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.r1_prefix} {params.r2_prefix} intersect {params.tmpdir}/{params.prefix}
        {input.kmc_info} {params.tmpdir}/{params.prefix} | tail -n 1 | awk '{{ print $3 }}' > {output}
        rm {params.tmpdir}/{params.prefix}*
        """

rule intersect_ct_db_pairs:
    input:
        r1="KMC_db/chm13v1_{region1}_{kmer}.kmc_pre",
        r2="KMC_db/chm13v1_{region2}_{kmer}.kmc_pre",
        kmc_tools="bin/kmc_tools",
        kmc_info="bin/kmc_info"
    output:
        "intersections_{kmer}/chm13v1_{region1}/{region2}.txt"
    params:
        r1_prefix="KMC_db/chm13v1_{region1}_{kmer}",
        r2_prefix="KMC_db/chm13v1_{region2}_{kmer}",
        prefix="{region1}-{region2}_{kmer}",
        tmpdir=TMPDIR
    wildcard_constraints:
        kmer=KMER_SET,
        region1=CT_REGION_SET,
        region2=CT_REGION_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.r1_prefix} {params.r2_prefix} intersect {params.tmpdir}/{params.prefix}
        {input.kmc_info} {params.tmpdir}/{params.prefix} | tail -n 1 | awk '{{ print $3 }}' > {output}
        rm {params.tmpdir}/{params.prefix}*
        """

def array_pair_input(wc):
    return "KMC_db/chm13v1_{}-{}_unique_{}.kmc_pre".format(wc.region.split('_')[0], wc.array, wc.kmer)

def array_pair_param(wc):
    return "KMC_db/chm13v1_{}-{}_unique_{}".format(wc.region.split('_')[0], wc.array, wc.kmer)

rule intersect_region_unique_arrays:
    input:
        r1="KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
        r2=array_pair_input,
        kmc_tools="bin/kmc_tools",
        kmc_info="bin/kmc_info"
    output:
        "intersections_{kmer}/chm13v1_{region}/{array}_unique.txt"
    params:
        r1_prefix="KMC_db/chm13v1_{region}_{kmer}",
        r2_prefix=array_pair_param,
        prefix="{region}-{array}_unique_{kmer}",
        tmpdir=TMPDIR
    wildcard_constraints:
        kmer=KMER_SET,
        region=REGION_SET,
        array=ARRAY_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.r1_prefix} {params.r2_prefix} \
            intersect {params.tmpdir}/{params.prefix}
        {input.kmc_info} {params.tmpdir}/{params.prefix} | tail -n 1 | awk '{{ print $3 }}' > {output}
        rm {params.tmpdir}/{params.prefix}*
        """

rule intersect_region_centro_pairs:
    input:
        r1="KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
        r2="KMC_db/chm13v1_centro_unique_{kmer}.kmc_pre",
        kmc_tools="bin/kmc_tools",
        kmc_info="bin/kmc_info"
    output:
        "intersections_{kmer}/chm13v1_{region}/centro_unique.txt"
    params:
        r1_prefix="KMC_db/chm13v1_{region}_{kmer}",
        r2_prefix="KMC_db/chm13v1_centro_unique_{kmer}",
        prefix="{region}-centro_{kmer}",
        tmpdir=TMPDIR
    wildcard_constraints:
        kmer=KMER_SET,
        region=REGION_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.r1_prefix} {params.r2_prefix} intersect {params.tmpdir}/{params.prefix}
        {input.kmc_info} {params.tmpdir}/{params.prefix} | tail -n 1 | awk '{{ print $3 }}' > {output}
        rm {params.tmpdir}/{params.prefix}*
        """

rule find_unique_kmers:
    input:
        r1="KMC_db/chm13v1_{region}_unique_{kmer}.kmc_pre",
        kmc_info="bin/kmc_info"
    output:
        "intersections_{kmer}/chm13v1_{region}/unique.txt"
    params:
        prefix="KMC_db/chm13v1_{region}_unique_{kmer}",
        tmpdir=TMPDIR
    wildcard_constraints:
        kmer=KMER_SET,
        region=REGION_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_info} {params.prefix} | tail -n 1 | awk '{{ print $3 }}' > {output}
        """

rule subcompile_intersections:
    input:
        "KMC_db/chm13v1_{region}_{kmer}.count",
        lambda wildcards: expand("intersections_{{kmer}}/chm13v1_{{region}}/{region2}.txt",
                                 region2=PAIR_DICT[wildcards.region]),
        expand("intersections_{{kmer}}/chm13v1_{{region}}/{array}_unique.txt", array=ARRAYS),
        "intersections_{kmer}/chm13v1_{region}/centro_unique.txt",
        "intersections_{kmer}/chm13v1_{region}/unique.txt"
    output:
        "intersections_{kmer}/chm13v1_{region}/compiled.txt"
    params:
        region="{region}",
        regions=NON_CT_REGIONS,
        arrays=ARRAYS,
        kmer="{kmer}"
    wildcard_constraints:
        kmer=KMER_SET,
        region=NON_CT_REGION_SET
    shell:
        """
        REGIONS="{params.regions}"
        for A in {params.arrays} centro; do
            REGIONS="$REGIONS ${{A}}_unique"
        done
        echo -n "{params.region}" > {output}
        for R in ${{REGIONS[*]}}; do
            if [[ $R == {params.region} ]]; then
                awk '{{ printf ",%s", $1 }}' KMC_db/chm13v1_{params.region}_{params.kmer}.count >> {output}
            else
                if [ -e intersections_{params.kmer}/chm13v1_{params.region}/${{R}}.txt ]; then
                    awk '{{ printf ",%s", $1 }}' intersections_{params.kmer}/chm13v1_{params.region}/${{R}}.txt >> {output}
                else
                    echo -n "," >> {output}
                fi
            fi
        done
        awk '{{ printf ",%s", $1 }}' intersections_{params.kmer}/chm13v1_{params.region}/unique.txt >> {output}
        awk '{{ printf ",%s", $1 }}' KMC_db/chm13v1_{params.region}_{params.kmer}.count >> {output}
        """

rule subcompile_ct_intersections:
    input:
        "KMC_db/chm13v1_{region}_{kmer}.count",
        lambda wildcards: expand("intersections_{{kmer}}/chm13v1_{{region}}/{region2}.txt",
                                 region2=CT_PAIR_DICT[wildcards.region]),
        expand("intersections_{{kmer}}/chm13v1_{{region}}/{array}_unique.txt", array=ARRAYS),
        "intersections_{kmer}/chm13v1_{region}/centro_unique.txt",
        "intersections_{kmer}/chm13v1_{region}/unique.txt"
    output:
        "intersections_{kmer}/chm13v1_{region}/compiled.txt"
    params:
        region="{region}",
        regions=CT_REGIONS,
        arrays=ARRAYS,
        kmer="{kmer}"
    wildcard_constraints:
        kmer=KMER_SET,
        region=CT_REGION_SET
    shell:
        """
        REGIONS="{params.regions}"
        for A in {params.arrays} centro; do
            REGIONS="$REGIONS ${{A}}_unique"
        done
        echo -n "{params.region}" > {output}
        for R in ${{REGIONS[*]}}; do
            if [[ $R == {params.region} ]]; then
                awk '{{ printf ",%s", $1 }}' KMC_db/chm13v1_{params.region}_{params.kmer}.count >> {output}
            else
                if [ -e intersections_{params.kmer}/chm13v1_{params.region}/${{R}}.txt ]; then
                    awk '{{ printf ",%s", $1 }}' intersections_{params.kmer}/chm13v1_{params.region}/${{R}}.txt >> {output}
                else
                    echo -n "," >> {output}
                fi
            fi
        done
        awk '{{ printf ",%s", $1 }}' intersections_{params.kmer}/chm13v1_{params.region}/unique.txt >> {output}
        awk '{{ printf ",%s", $1 }}' KMC_db/chm13v1_{params.region}_{params.kmer}.count >> {output}
        """

rule compile_intersections:
    input:
        expand("intersections_{{kmer}}/chm13v1_{region}/compiled.txt", region=NON_CT_REGIONS),
    output:
        "results/chm13v1_intersections_{kmer}.txt"
    params:
        regions=",".join(NON_CT_REGIONS),
        arrays=",".join(ARRAYS),
        kmer="{kmer}"
    wildcard_constraints:
        kmer=KMER_SET
    shell:
        """
        bin/compile_intersections.py {params.kmer} {params.regions} {params.arrays} {output}
        """

rule compile_ct_intersections:
    input:
        expand("intersections_{{kmer}}/chm13v1_{region}/compiled.txt", region=CT_REGIONS),
    output:
        "results/chm13v1_intersections_ct_{kmer}.txt"
    params:
        regions=",".join(CT_REGIONS),
        arrays=",".join(ARRAYS),
        kmer="{kmer}"
    wildcard_constraints:
        kmer=KMER_SET
    shell:
        """
        bin/compile_intersections.py {params.kmer} {params.regions} {params.arrays} {output}
        """

rule get_kmer_coverage:
    input:
        "Kmer_BigWigs/chm13v1_{kmer}mer.bw"
    output:
        "results/chm13v1_{kmer}mer_coverage.bed"
    params:
        kmer="{kmer}"
    wildcard_constraints:
        kmer=ALL_KMER_SET
    shell:
        """
        bin/get_unique_kmer_coverage.py {input} {output} {params.kmer}
        """


######## Plot results #########

def plot_input(wc):
    inputs = {"anno": "data/t2t_cenSatAnnotation.bed"}
    if wc.array == 'ct':
        inputs["results"] = "results/chm13v1_intersections_ct_{}.txt".format(wc.kmer)
    else:
        inputs["results"] = "results/chm13v1_intersections_{}.txt".format(wc.kmer)
    return inputs

rule plot_circos_data:
    input:
        unpack(plot_input)
    output:
        "plots/circos_{array}_{kmer}.svg"
    params:
        kmer="{kmer}",
        array="{array}"
    wildcard_constraints:
        kmer=KMER_SET,
        array="all|hor|hsat|bsat|censat|mon|mismatch|ct"
    log:
        "logs/circos/format_data_{array}_{kmer}.txt"
    conda:
        "envs/circos.yaml"
    shell:
        """
        mkdir -p circos_data
        bin/format_circos_data.py {input.anno} {params.kmer} {params.array} {input.results}
        cat circos_etc/circos.conf | \
            sed "s/KMER/{params.kmer}/g" | \
            sed "s/ARRAY/{params.array}/g" > \
            circos_etc/circos_{params.array}_{params.kmer}.conf
        circos -conf circos_etc/circos_{params.array}_{params.kmer}.conf
        """


########## Get 100-kmers for HOR combinations ##############

# 1. Identify HOR-specific Kmers KMC_db/chm13v1_HOR-HOR_unique_100
# 2. Identify all Kmers shared for a given HOR pair
# 3. Identify non-HOR-specific Kmers
# 4. Identify pair-specific Kmers
# 5. Identify HOR-specific non-pair-specific Kmers

rule create_full_HOR_kmc_db:
    input:
        "KMC_db/chm13v1_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{kmer}.kmc_suf",
    output:
        "HOR_KMC_db/chm13v1_{kmer}.kmc_pre",
        "HOR_KMC_db/chm13v1_{kmer}.kmc_suf"
    params:
        outdir="HOR_KMC_db/"
    wildcard_constraints:
        kmer=HOR_KMER_SET
    shell:
        """
        ln -s ${{PWD}}/{input[0]} {params.outdir}
        ln -s ${{PWD}}/{input[1]} {params.outdir}
        """

rule extract_HOR_array_regions:
    input:
        fa="fasta/chm13v1.fasta",
        fai="fasta/chm13v1.fasta.fai",
        anno="data/t2t_cenSatAnnotation.bed"
    output:
        "HOR_regions/chm13v1_{region}.fasta"
    wildcard_constraints:
        region=HOR_REGION_SET
    params:
        region="{region}"
    conda:
        "envs/deepsam.yaml"
    shell:
        """
        region=$(eval grep -w "{params.region}" {input.anno} | awk '{{ printf "%s:%s-%s", $1, $2, $3 }}')
        samtools faidx {input.fa} $region > {output}
        """

rule extract_hor_array:
    input:
        fa="fasta/chm13v1.fasta",
        fai="fasta/chm13v1.fasta.fai",
        anno="data/t2t_cenSatAnnotation.bed",
    output:
        "HOR_regions/chm13v1_HOR_all.fasta"
    conda:
        "envs/deepsam.yaml"
    shell:
        """
        rm -f tmp.HOR.regions
        regions=($(eval grep hor {input.anno} | awk '{{ printf "%s:%s-%s ", $1, $2, $3 }}'))
        for R in ${{regions[*]}}; do
            echo $R >> tmp.HOR.regions
        done
        samtools faidx {input.fa} -r tmp.HOR.regions > {output}
        rm -f tmp.HOR.regions
        """

rule create_HOR_region_kmc_db:
    input:
        fa="HOR_regions/chm13v1_{region}.fasta",
        kmc="bin/kmc"
    output:
        "HOR_KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
        "HOR_KMC_db/chm13v1_{region}_{kmer}.kmc_suf"
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="HOR_KMC_db/chm13v1_{region}_{kmer}"
    threads:
        MAXTHREADS
    wildcard_constraints:
        kmer=HOR_KMER_SET,
        region=HOR_REGION_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc} -k{params.kmer} -m{params.maxmem} -fm -ci1 -t{threads} {input.fa} {params.prefix} {params.tmpdir}
        """

rule create_HOR_kmc_db:
    input:
        fa="HOR_regions/chm13v1_HOR_all.fasta",
        kmc="bin/kmc"
    output:
        "HOR_KMC_db/chm13v1_HOR_all_{kmer}.kmc_pre",
        "HOR_KMC_db/chm13v1_HOR_all_{kmer}.kmc_suf"
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="HOR_KMC_db/chm13v1_HOR_all_{kmer}"
    threads:
        MAXTHREADS
    wildcard_constraints:
        kmer=HOR_KMER_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc} -k{params.kmer} -m{params.maxmem} -fm -ci1 -t{threads} {input.fa} {params.prefix} {params.tmpdir}
        """

rule create_HOR_unique_kmc_db:
    input:
        kmc_tools="bin/kmc_tools",
        full_db="HOR_KMC_db/chm13v1_{kmer}.kmc_pre",
        HOR_db="HOR_KMC_db/chm13v1_HOR_all_{kmer}.kmc_pre"
    output:
        "HOR_KMC_db/chm13v1_HOR_all_unique_{kmer}.kmc_pre",
        "HOR_KMC_db/chm13v1_HOR_all_unique_{kmer}.kmc_suf"
    params:
        prefix="HOR_KMC_db/chm13v1_HOR_all_unique_{kmer}",
        HOR_prefix="HOR_KMC_db/chm13v1_HOR_all_{kmer}",
        full_prefix="HOR_KMC_db/chm13v1_{kmer}"
    wildcard_constraints:
        kmer=HOR_KMER_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.HOR_prefix} {params.full_prefix} counters_compare {params.prefix}
        """

rule create_HOR_nonunique_kmc_db:
    input:
        kmc_tools="bin/kmc_tools",
        unique_db="HOR_KMC_db/chm13v1_HOR_all_unique_{kmer}.kmc_pre",
        HOR_db="HOR_KMC_db/chm13v1_HOR_all_{kmer}.kmc_pre"
    output:
        "HOR_KMC_db/chm13v1_HOR_all_nonunique_{kmer}.kmc_pre",
        "HOR_KMC_db/chm13v1_HOR_all_nonunique_{kmer}.kmc_suf"
    params:
        prefix="HOR_KMC_db/chm13v1_HOR_all_nonunique_{kmer}",
        HOR_prefix="HOR_KMC_db/chm13v1_HOR_all_{kmer}",
        unique_prefix="HOR_KMC_db/chm13v1_HOR_all_unique_{kmer}"
    wildcard_constraints:
        kmer=HOR_KMER_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.HOR_prefix} {params.unique_prefix} kmers_subtract {params.prefix}
        """

rule intersect_HOR_unique_kmc_db:
    input:
        kmc_tools="bin/kmc_tools",
        kmc_dump="bin/kmc_dump",
        unique_db="HOR_KMC_db/chm13v1_HOR_all_unique_{kmer}.kmc_pre",
        region_db="HOR_KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
    output:
        "HOR_{kmer}/chm13v1_{region}/HORspecific.txt"
    params:
        region_prefix="HOR_KMC_db/chm13v1_{region}_{kmer}",
        unique_prefix="HOR_KMC_db/chm13v1_HOR_all_unique_{kmer}",
        prefix="chm13v1_{region}_unique_{kmer}",
        tmpdir=TMPDIR
    wildcard_constraints:
        kmer=HOR_KMER_SET,
        region=HOR_REGION_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.region_prefix} {params.unique_prefix} \
            intersect {params.tmpdir}/{params.prefix}
        {input.kmc_dump} {params.tmpdir}/{params.prefix} {output}
        """

rule intersect_HOR_nonunique_kmc_db:
    input:
        kmc_tools="bin/kmc_tools",
        kmc_dump="bin/kmc_dump",
        nonunique_db="HOR_KMC_db/chm13v1_HOR_all_nonunique_{kmer}.kmc_pre",
        region_db="HOR_KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
    output:
        "HOR_{kmer}/chm13v1_{region}/general.txt"
    params:
        region_prefix="HOR_KMC_db/chm13v1_{region}_{kmer}",
        nonunique_prefix="HOR_KMC_db/chm13v1_HOR_all_nonunique_{kmer}",
        prefix="chm13v1_{region}_nonunique_{kmer}",
        tmpdir=TMPDIR
    wildcard_constraints:
        kmer=HOR_KMER_SET,
        region=HOR_REGION_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.region_prefix} {params.nonunique_prefix} \
            intersect {params.tmpdir}/{params.prefix}
        {input.kmc_dump} {params.tmpdir}/{params.prefix} {output}
        """

rule compile_HOR_intersections:
    input:
        expand("HOR_{{kmer}}/chm13v1_{region}/{group}.txt",
               region=HOR_REGIONS, group=["general", "HORspecific"]),
    output:
        "results/chm13v1_HOR_{kmer}mers.txt"
    params:
        regions=HOR_REGIONS,
        groups=["general","HORspecific"],
        prefix="HOR_{kmer}/chm13v1_"
    wildcard_constraints:
        kmer=HOR_KMER_SET
    shell:
        """
        REGIONS=({params.regions})
        R="${{REGIONS[0]}}"
        for I in $(seq 1 1 $(expr ${{#REGIONS[*]}} - 1)); do
            R="${{R}},${{REGIONS[${{I}}]}}"
        done
        G1=({params.groups})
        G="${{G1[0]}}"
        for I in $(seq 1 1 $(expr ${{#G1[*]}} - 1)); do
            G="${{G}},${{G1[${{I}}]}}"
        done
        bin/compile_array_kmers.py {params.prefix} $R $G > {output}
        """


########## Get 21-kmers for bsat combinations ##############

# 1. Identify BSAT-specific Kmers KMC_db/chm13v1_BSAT-BSAT_unique_100
# 2. Identify all Kmers shared for a given BSAT pair
# 3. Identify non-BSAT-specific Kmers
# 4. Identify pair-specific Kmers
# 5. Identify BSAT-specific non-pair-specific Kmers

rule create_full_BSAT_kmc_db:
    input:
        "KMC_db/chm13v1_{kmer}.kmc_pre",
        "KMC_db/chm13v1_{kmer}.kmc_suf",
    output:
        "BSAT_KMC_db/chm13v1_{kmer}.kmc_pre",
        "BSAT_KMC_db/chm13v1_{kmer}.kmc_suf"
    params:
        outdir="BSAT_KMC_db/"
    wildcard_constraints:
        kmer=BSAT_KMER_SET
    shell:
        """
        ln -s ${{PWD}}/{input[0]} {params.outdir}
        ln -s ${{PWD}}/{input[1]} {params.outdir}
        """

rule extract_BSAT_array_regions:
    input:
        fa="fasta/chm13v1.fasta",
        fai="fasta/chm13v1.fasta.fai",
        anno="data/t2t_cenSatAnnotation.bed"
    output:
        "BSAT_regions/chm13v1_{region}.fasta"
    wildcard_constraints:
        region=BSAT_REGION_SET
    params:
        region="{region}"
    conda:
        "envs/deepsam.yaml"
    shell:
        """
        region=$(eval grep -w "{params.region}" {input.anno} | awk '{{ printf "%s:%s-%s", $1, $2, $3 }}')
        samtools faidx {input.fa} $region > {output}
        """

rule extract_bsat_array:
    input:
        fa="fasta/chm13v1.fasta",
        fai="fasta/chm13v1.fasta.fai",
        anno="data/t2t_cenSatAnnotation.bed",
    output:
        "BSAT_regions/chm13v1_BSAT_all.fasta"
    conda:
        "envs/deepsam.yaml"
    shell:
        """
        rm -f tmp.BSAT.regions
        regions=($(eval grep bsat {input.anno} | awk '{{ printf "%s:%s-%s ", $1, $2, $3 }}'))
        for R in ${{regions[*]}}; do
            echo $R >> tmp.BSAT.regions
        done
        samtools faidx {input.fa} -r tmp.BSAT.regions > {output}
        rm -f tmp.BSAT.regions
        """

rule create_BSAT_region_kmc_db:
    input:
        fa="BSAT_regions/chm13v1_{region}.fasta",
        kmc="bin/kmc"
    output:
        "BSAT_KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
        "BSAT_KMC_db/chm13v1_{region}_{kmer}.kmc_suf"
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="BSAT_KMC_db/chm13v1_{region}_{kmer}"
    threads:
        MAXTHREADS
    wildcard_constraints:
        kmer=BSAT_KMER_SET,
        region=BSAT_REGION_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc} -k{params.kmer} -m{params.maxmem} -fm -ci1 -t{threads} {input.fa} {params.prefix} {params.tmpdir}
        """

rule create_BSAT_kmc_db:
    input:
        fa="BSAT_regions/chm13v1_BSAT_all.fasta",
        kmc="bin/kmc"
    output:
        "BSAT_KMC_db/chm13v1_BSAT_all_{kmer}.kmc_pre",
        "BSAT_KMC_db/chm13v1_BSAT_all_{kmer}.kmc_suf"
    params:
        tmpdir=TMPDIR,
        kmer="{kmer}",
        maxmem=MAXMEM,
        prefix="BSAT_KMC_db/chm13v1_BSAT_all_{kmer}"
    threads:
        MAXTHREADS
    wildcard_constraints:
        kmer=BSAT_KMER_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc} -k{params.kmer} -m{params.maxmem} -fm -ci1 -t{threads} {input.fa} {params.prefix} {params.tmpdir}
        """

rule create_BSAT_unique_kmc_db:
    input:
        kmc_tools="bin/kmc_tools",
        full_db="BSAT_KMC_db/chm13v1_{kmer}.kmc_pre",
        BSAT_db="BSAT_KMC_db/chm13v1_BSAT_all_{kmer}.kmc_pre"
    output:
        "BSAT_KMC_db/chm13v1_BSAT_all_unique_{kmer}.kmc_pre",
        "BSAT_KMC_db/chm13v1_BSAT_all_unique_{kmer}.kmc_suf"
    params:
        prefix="BSAT_KMC_db/chm13v1_BSAT_all_unique_{kmer}",
        BSAT_prefix="BSAT_KMC_db/chm13v1_BSAT_all_{kmer}",
        full_prefix="BSAT_KMC_db/chm13v1_{kmer}"
    wildcard_constraints:
        kmer=BSAT_KMER_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.BSAT_prefix} {params.full_prefix} counters_compare {params.prefix}
        """

rule create_BSAT_nonunique_kmc_db:
    input:
        kmc_tools="bin/kmc_tools",
        unique_db="BSAT_KMC_db/chm13v1_BSAT_all_unique_{kmer}.kmc_pre",
        BSAT_db="BSAT_KMC_db/chm13v1_BSAT_all_{kmer}.kmc_pre"
    output:
        "BSAT_KMC_db/chm13v1_BSAT_all_nonunique_{kmer}.kmc_pre",
        "BSAT_KMC_db/chm13v1_BSAT_all_nonunique_{kmer}.kmc_suf"
    params:
        prefix="BSAT_KMC_db/chm13v1_BSAT_all_nonunique_{kmer}",
        BSAT_prefix="BSAT_KMC_db/chm13v1_BSAT_all_{kmer}",
        unique_prefix="BSAT_KMC_db/chm13v1_BSAT_all_unique_{kmer}"
    wildcard_constraints:
        kmer=BSAT_KMER_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.BSAT_prefix} {params.unique_prefix} kmers_subtract {params.prefix}
        """

rule intersect_BSAT_unique_kmc_db:
    input:
        kmc_tools="bin/kmc_tools",
        kmc_dump="bin/kmc_dump",
        unique_db="BSAT_KMC_db/chm13v1_BSAT_all_unique_{kmer}.kmc_pre",
        region_db="BSAT_KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
    output:
        "BSAT_{kmer}/chm13v1_{region}/BSATspecific.txt"
    params:
        region_prefix="BSAT_KMC_db/chm13v1_{region}_{kmer}",
        unique_prefix="BSAT_KMC_db/chm13v1_BSAT_all_unique_{kmer}",
        prefix="chm13v1_{region}_unique_{kmer}",
        tmpdir=TMPDIR
    wildcard_constraints:
        kmer=BSAT_KMER_SET,
        region=BSAT_REGION_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.region_prefix} {params.unique_prefix} \
            intersect {params.tmpdir}/{params.prefix}
        {input.kmc_dump} {params.tmpdir}/{params.prefix} {output}
        """

rule intersect_BSAT_nonunique_kmc_db:
    input:
        kmc_tools="bin/kmc_tools",
        kmc_dump="bin/kmc_dump",
        nonunique_db="BSAT_KMC_db/chm13v1_BSAT_all_nonunique_{kmer}.kmc_pre",
        region_db="BSAT_KMC_db/chm13v1_{region}_{kmer}.kmc_pre",
    output:
        "BSAT_{kmer}/chm13v1_{region}/general.txt"
    params:
        region_prefix="BSAT_KMC_db/chm13v1_{region}_{kmer}",
        nonunique_prefix="BSAT_KMC_db/chm13v1_BSAT_all_nonunique_{kmer}",
        prefix="chm13v1_{region}_nonunique_{kmer}",
        tmpdir=TMPDIR
    wildcard_constraints:
        kmer=BSAT_KMER_SET,
        region=BSAT_REGION_SET
    conda:
        "envs/kmc.yaml"
    shell:
        """
        {input.kmc_tools} simple {params.region_prefix} {params.nonunique_prefix} \
            intersect {params.tmpdir}/{params.prefix}
        {input.kmc_dump} {params.tmpdir}/{params.prefix} {output}
        """

rule compile_BSAT_intersections:
    input:
        expand("BSAT_{{kmer}}/chm13v1_{region}/{group}.txt",
               region=BSAT_REGIONS, group=["general", "BSATspecific"]),
    output:
        "results/chm13v1_BSAT_{kmer}mers.txt"
    params:
        regions=BSAT_REGIONS,
        groups=["general","BSATspecific"],
        prefix="BSAT_{kmer}/chm13v1_"
    wildcard_constraints:
        kmer=BSAT_KMER_SET
    shell:
        """
        REGIONS=({params.regions})
        R="${{REGIONS[0]}}"
        for I in $(seq 1 1 $(expr ${{#REGIONS[*]}} - 1)); do
            R="${{R}},${{REGIONS[${{I}}]}}"
        done
        G1=({params.groups})
        G="${{G1[0]}}"
        for I in $(seq 1 1 $(expr ${{#G1[*]}} - 1)); do
            G="${{G}},${{G1[${{I}}]}}"
        done
        bin/compile_array_kmers.py {params.prefix} $R $G > {output}
        """

