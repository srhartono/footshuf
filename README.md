#1. SYNOPSIS

---

##1.1 Shuffle

`perl 1_random_peak.pl -i <ENERGY_PIPELINE_PEAK> -l <LABEL> -r <SHUFFLE_NUMBER> -o <PEAK1_SHUFFLE_OUTPUT>`

- -i: <ENERGY_PIPELINE_PEAK> a peak file from your energy pipeline
- -l: <LABEL> user-defined [alphanumeric]
- -r: <SHUFFLE_NUMBER> number of desired shuffles [integer, default 1000]
- -o: <PEAK1_SHUFFLE_OUTPUT> name of output directory that'll contain shuffled peaks

---

##1.2 Intersect and statistics

`perl 2_shuffle_wrapper.pl -x -i 3_intersect.pl -a <PEAK1_SHUFFLE_OUTPUT> -b <FOLDER_OF_FOOTLOOP_PEAKS> -o <OUTPUTDIR>`

- -i: <3_intersect.pl> literally the location of 3_intersect.pl if it differs than current folder
- -a: <PEAK1_SHUFFLE_OUTPUT> output directory defined in shuffle pipeline above (basically -o of 1_random_peak.pl)
- -b: <FOLDER_OF_FOOTLOOP_PEAKS> folder containing footloop peaks
- -o: <OUTPUTDIR> name of output directory tat'll contain statistics and intersect result
- -x: dry run 


If you do dry run (-x), see <OUTPUTDIR>/<LABEL>_run.sh in for scripts of what's going to be run.

---

#1.3.1 <ENERGY_PIPELINE_PEAK> format expectations.

**Please avoid using non-alphanumeric characters in filename or directories.**

Use [a-z0-9_] to be safe, as while special characters is generally fine, it can and does cause unexpected bugs in the past.

- Good names: `ea_neg7p_bbprob.bed`
- Danger zone: `ea -0.07% 3.6kbp bbprob.bed`

The script expects multiple gene blocks, where each gene block contains:

1. **One row** of browser position line with amplicon coordinate (space separated)
   - `browser position <chr.amplicon>:<start.amplicon>-<end.amplicon>`
2. **One row** of hashtag followed by gene name
   - `#<GENE>`
3. **Multiple rows** of tab-separated bed6\*.
   - `<chr.amplicon.1>   <start.peak.1>   <end.peak.1>   <name1.peak.1>   <value.peak.1>   <strand.peak.1>`
   - `<chr.amplicon.2>   <start.peak.2>   <end.peak.2>   <name2.peak.2>   <value.peak.2>   <strand.peak.2>`


*BED6 means 6 column bed format (see ![BED format explanation from UCSC](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).


Example:

```
browser position <chr.amplicon>:<start.amplicon>-<end.amplicon>
#<GENE.1>
<chr.amplicon.1>   <start.peak.1>   <end.peak.1>   <name.peak.1>   <value.peak.1>   <strand.peak.1>
<chr.amplicon.2>   <start.peak.2>   <end.peak.2>   <name.peak.2>   <value.peak.2>   <strand.peak.2>
browser position <chr.amplicon>:<start.amplicon>-<end.amplicon>
#<GENE.2>
<chr.amplicon.1>   <start.peak.1>   <end.peak.1>   <name.peak.1>   <value.peak.1>   <strand.peak.1>
<chr.amplicon.2>   <start.peak.2>   <end.peak.2>   <name.peak.2>   <value.peak.2>   <strand.peak.2>
```

- <value.peak> (5th column) is currently not used, but it's necessary for intersect purpose (as various softwares used in this script e.g. `bedtools`, expects 6 columns in bed6).
- Currently, all peaks are weighted equally, but if later you'd like to weight peak differently, <value.peak> is where you should store the weight of each peak.
- In case this isn't already clear: variable names (e.g. chr.amplicon) shouldn't be enclosed with "<" and ">", the above is explicitly shown like that only for example purposes.


**Important notes:**

- In general, variable names such as <GENE1> and <chr.amplicon> are generally case-sensitive and has to match those in data table at section 2.1. The script will still work, but the script will always produce output with names matching those in the table. The names in data table correspond to UCSC, which we use a lot, therefore it's best to adhere to UCSC naming system to prevent unexpected bugs.
- <peakname.gene> has to be unique within a gene block.

#1.3.2 Format examples

## A. Good example

```
browser position chr5:72794135-72797124
#BTF3
chr5  72794135 72794226 POS0  0  +
chr5  72794344 72794947 POS1  0  +
chr5  72795105 72795247 POS2  0  +
chr5  72795440 72795544 POS3  0  +
chr5  72796444 72796498 POS4  0  +
browser position chr3:128900662-128903317
#CNBP
chr3  128900788   128900850   NEG0  0  -
chr3  128900883   128900967   NEG1  0  -
chr3  128901151   128901214   NEG2  0  -
chr3  128903250   128903316   NEG3  0  -
browser position chr19:3982336-3985724
#EEF2
chr19 3982340  3982408  NEG0  0  -
chr19 3982445  3983188  NEG1  0  -
chr19 3984929  3985408  NEG2  0  -
chr19 3985565  3985654  NEG3  0  -
```

## B. Bad examples:

### B.1. <peakname.gene> is not unique as POS0 appears twice

```
browser position chr5:72794135-72797124
#BTF3
chr5  72794135 72794226 **POS0 ** 0  +
chr5  72794344 72794947 **POS0**  0  +
chr5  72795700 72795781 POS1  0  +
```

### B.2. <gene> block appears twice

```
browser position chr5:72794135-72797124
**#BTF3**
chr5  72794135 72794226 POS0  0  +
chr5  72795700 72795781 POS1  0  +
browser position chr3:128900662-128903317
#CNBP
chr3  128900788   128900850   NEG0  0  -
chr3  128900883   128900967   NEG1  0  -
chr3  128901151   128901214   NEG2  0  -
chr3  128903250   128903316   NEG3  0  -
browser position chr19:3982336-3985724
**#BTF3**
chr5 3982340  3982408  NEG0  0  -
chr5 3982445  3983188  NEG1  0  -
chr5 3984929  3985408  NEG2  0  -
chr5 3985565  3985654  NEG3  0  -
```

### B.3. <gene> isn't in amplicon name data table (see section 2.1)

```
browser position chr5:72794135-72797124
**#BTF3-P**
chr5  72794135 72794226 POS0  0  +
chr5  72795700 72795781 POS1  0  +
```

### B.4. <gene> and <chr.amplicon> is in amplicon name data table, but wrong case

```
browser position **ChR5**:72794135-72797124
**#btf3**
chr5  72794135 72794226 POS0  0  +
chr5  72795700 72795781 POS1  0  +
```

---

#2.1 APRIL 2018 amplicon data

## 2.1.1 amplicon names (<chr.amplicon>)

### A. IN VIVO

```
BTF3
CALM3
CEP95
CLOCK
CNBP
CPZ
EEF2
FAM219B
FUS
GADD45A
KAT5
MORF4L2
MRFAP1L1
MTERFD3
MYLIP
PIN4
POGK
PPDPF
PRPF38B
RPL13A
RPL3
RPL4
RPS24
SETD5
SNRPN70
TADA1
TOMM20
TOMM40
TRIM33
WDR3
5SrDNA3
18SrDNA1
28SdrDNA6
28SprDNA1
```

### B. IN VITRO

```
pFC53_AIRN_REVERSE
PFC8_SNRPN_REVERSE
pFC66_AIRN_FORWARD
```

## 2.1.2 amplicon coordinates

### A. IN VIVO

```
chr5	72794135	72797124	BTF3	0	+
chr19	47111177	47113789	CALM3	0	+
chr17	62531459	62535650	CEP95	0	+
chr4	56411571	56414807	CLOCK	0	-
chr3	128900661	128903317	CNBP	0	-
chr4	8613737	8616343	CPZ	0	-
chr19	3982335	3985724	EEF2	0	-
chr15	75190832	75194428	FAM219B	0	-
chr16	31200410	31202996	FUS	0	+
chr1	68149934	68153573	GADD45A	0	+
chr11	65482845	65485619	KAT5	0	+
chrX	102929428	102932836	MORF4L2	0	-
chr4	6708655	6711808	MRFAP1L1	0	-
chr12	107376879	107381257	MTERFD3	0	-
chr6	16147958	16151697	MYLIP	0	+
chrX	71401153	71404487	PIN4	0	+
chr1	166818107	166821100	POGK	0	+
chr20	62153396	62156401	PPDPF	0	+
chr1	109242027	109245458	PRPF38B	0	+
chr19	49993082	49995655	RPL13A	0	+
chr22	39712418	39716108	RPL3	0	-
chr15	66794648	66797585	RPL4	0	-
chr10	79799064	79801720	RPS24	0	+
chr3	9519703	9522929	SETD5	0	+
chr19	49609997	49612943	SNRPN70	0	+
chr1	166824952	166827904	TADA1	0	-
chr1	235272022	235275548	TOMM20	0	-
chr19	45405433	45408838	TOMM40	0	+
chr1	114931599	114934994	TRIM33	0	-
chr1	118471462	118475880	WDR3	0	+
U13369	305	2711	5SrDNA3	+
U13369	3736	5508	18SrDNA1	+
U13369	7984	10780	28SdrDNA6	+
U13369	11613	12888	28SprDNA1	+
```

### B. IN VITRO

```
pFC53FIXED	80	1829	pFC53_AIRN_REVERSE	0	-
pFC19FIXED	80	1512	PFC8_SNRPN_REVERSE	0	-
pFC66FIXED	212	2180	pFC66_AIRN_FORWARD	0	+
```
