###
# Stage test data
from contextlib import contextmanager
import numpy
from numpy import random
from pathlib import Path
import itertools
from tempfile import TemporaryDirectory


STAR_2_7_9a_log = """STAR version=dev_EoI_2.7.9a_2021-09-30
STAR compilation time,server,dir=2021-09-30T17:07:55-07:00 :/tmp/tmp.01I63TZvRF/STAR/source
STAR git: On branch dev_ExonOverIntron ; commit 12beb0c52367d568fc993c1990795676e7f9d9cb ; diff files: 
##### Command Line:
STAR --genomeDir genome --readFilesIn SRX5908538_R2.fastq.gz SRX5908538_R1.fastq.gz --readFilesCommand zcat --runThreadN 16 --genomeLoad NoSharedMemory --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD CB CR CY UB UR UY gx gn --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --sjdbScore 1 --clipAdapterType CellRanger4 --outFilterScoreMin 30 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIdedup 1MM_CR --soloUMIfiltering MultiGeneUMI_CR --soloType CB_UMI_Simple --soloCellFilter EmptyDrops_CR --soloUMIlen 10 --soloCBlen 16 --soloBarcodeReadLength 0 --soloCBwhitelist 10xv2_allowlist.txt --soloStrand Forward --soloFeatures GeneFull_Ex50pAS SJ --soloMultiMappers Unique EM --limitBAMsortRAM 51539607552 --outTmpDir _STARtmp --outFileNamePrefix ./
##### Initial USER parameters from Command Line:
outFileNamePrefix                 ./
outTmpDir                         _STARtmp
###### All USER parameters from Command Line:
genomeDir                     genome     ~RE-DEFINED
readFilesIn                   SRX5908538_R2.fastq.gz   SRX5908538_R1.fastq.gz        ~RE-DEFINED
readFilesCommand              zcat        ~RE-DEFINED
runThreadN                    16     ~RE-DEFINED
genomeLoad                    NoSharedMemory     ~RE-DEFINED
outFilterMultimapNmax         20     ~RE-DEFINED
alignSJoverhangMin            8     ~RE-DEFINED
alignSJDBoverhangMin          1     ~RE-DEFINED
outFilterMismatchNmax         999     ~RE-DEFINED
outFilterMismatchNoverReadLmax0.04     ~RE-DEFINED
alignIntronMin                20     ~RE-DEFINED
alignIntronMax                1000000     ~RE-DEFINED
alignMatesGapMax              1000000     ~RE-DEFINED
outSAMheaderCommentFile       COfile.txt     ~RE-DEFINED
outSAMheaderHD                @HD   VN:1.4   SO:coordinate        ~RE-DEFINED
outSAMunmapped                Within        ~RE-DEFINED
outFilterType                 BySJout     ~RE-DEFINED
outSAMattributes              NH   HI   AS   NM   MD   CB   CR   CY   UB   UR   UY   gx   gn        ~RE-DEFINED
outSAMstrandField             intronMotif     ~RE-DEFINED
outSAMtype                    BAM   SortedByCoordinate        ~RE-DEFINED
sjdbScore                     1     ~RE-DEFINED
clipAdapterType               CellRanger4        ~RE-DEFINED
outFilterScoreMin             30     ~RE-DEFINED
soloCBmatchWLtype             1MM_multi_Nbase_pseudocounts     ~RE-DEFINED
soloUMIdedup                  1MM_CR        ~RE-DEFINED
soloUMIfiltering              MultiGeneUMI_CR        ~RE-DEFINED
soloType                      CB_UMI_Simple     ~RE-DEFINED
soloCellFilter                EmptyDrops_CR        ~RE-DEFINED
soloUMIlen                    10     ~RE-DEFINED
soloCBlen                     16     ~RE-DEFINED
soloBarcodeReadLength         0     ~RE-DEFINED
soloCBwhitelist               10xv2_allowlist.txt        ~RE-DEFINED
soloStrand                    Forward     ~RE-DEFINED
soloFeatures                  GeneFull_Ex50pAS   SJ        ~RE-DEFINED
soloMultiMappers              Unique   EM        ~RE-DEFINED
limitBAMsortRAM               51539607552     ~RE-DEFINED
outTmpDir                     _STARtmp     ~RE-DEFINED
outFileNamePrefix             ./     ~RE-DEFINED
##### Finished reading parameters from all sources

##### Final user re-defined parameters-----------------:
runThreadN                        16
genomeDir                         genome
genomeLoad                        NoSharedMemory
readFilesIn                       SRX5908538_R2.fastq.gz   SRX5908538_R1.fastq.gz   
readFilesCommand                  zcat   
limitBAMsortRAM                   51539607552
outFileNamePrefix                 ./
outTmpDir                         _STARtmp
outSAMtype                        BAM   SortedByCoordinate   
outSAMstrandField                 intronMotif
outSAMattributes                  NH   HI   AS   NM   MD   CB   CR   CY   UB   UR   UY   gx   gn   
outSAMunmapped                    Within   
outSAMheaderHD                    @HD   VN:1.4   SO:coordinate   
outSAMheaderCommentFile           COfile.txt
outFilterType                     BySJout
outFilterMultimapNmax             20
outFilterScoreMin                 30
outFilterMismatchNmax             999
outFilterMismatchNoverReadLmax    0.04
clipAdapterType                   CellRanger4   
alignIntronMin                    20
alignIntronMax                    1000000
alignMatesGapMax                  1000000
alignSJoverhangMin                8
alignSJDBoverhangMin              1
sjdbScore                         1
soloType                          CB_UMI_Simple
soloCBlen                         16
soloUMIlen                        10
soloBarcodeReadLength             0
soloCBwhitelist                   10xv2_allowlist.txt   
soloStrand                        Forward
soloFeatures                      GeneFull_Ex50pAS   SJ   
soloUMIdedup                      1MM_CR   
soloCBmatchWLtype                 1MM_multi_Nbase_pseudocounts
soloCellFilter                    EmptyDrops_CR   
soloUMIfiltering                  MultiGeneUMI_CR   
soloMultiMappers                  Unique   EM   

-------------------------------
##### Final effective command line:
STAR   --runThreadN 16   --genomeDir genome   --genomeLoad NoSharedMemory   --readFilesIn SRX5908538_R2.fastq.gz   SRX5908538_R1.fastq.gz      --readFilesCommand zcat      --limitBAMsortRAM 51539607552   --outFileNamePrefix ./   --outTmpDir _STARtmp   --outSAMtype BAM   SortedByCoordinate      --outSAMstrandField intronMotif   --outSAMattributes NH   HI   AS   NM   MD   CB   CR   CY   UB   UR   UY   gx   gn      --outSAMunmapped Within      --outSAMheaderHD @HD   VN:1.4   SO:coordinate      --outSAMheaderCommentFile COfile.txt   --outFilterType BySJout   --outFilterMultimapNmax 20   --outFilterScoreMin 30   --outFilterMismatchNmax 999   --outFilterMismatchNoverReadLmax 0.04   --clipAdapterType CellRanger4      --alignIntronMin 20   --alignIntronMax 1000000   --alignMatesGapMax 1000000   --alignSJoverhangMin 8   --alignSJDBoverhangMin 1   --sjdbScore 1   --soloType CB_UMI_Simple   --soloCBlen 16   --soloUMIlen 10   --soloBarcodeReadLength 0   --soloCBwhitelist 10xv2_allowlist.txt      --soloStrand Forward   --soloFeatures GeneFull_Ex50pAS   SJ      --soloUMIdedup 1MM_CR      --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts   --soloCellFilter EmptyDrops_CR      --soloUMIfiltering MultiGeneUMI_CR      --soloMultiMappers Unique   EM   
----------------------------------------

Number of fastq files for each mate = 1

   Input read files for mate 1 :
-rw-r--r-- 1 diane diane 9242999357 Oct 11 16:23 SRX5908538_R2.fastq.gz

   readsCommandsFile:
exec > "_STARtmp/tmp.fifo.read1"
echo FILE 0
zcat      "SRX5908538_R2.fastq.gz"


   Input read files for mate 2 :
-rw-r--r-- 1 diane diane 2888618147 Oct 11 16:20 SRX5908538_R1.fastq.gz

   readsCommandsFile:
exec > "_STARtmp/tmp.fifo.read2"
echo FILE 0
zcat      "SRX5908538_R1.fastq.gz"

WARNING --outSAMstrandField=intronMotif, therefore STAR will output XS attribute
ParametersSolo: using hardcoded filtering parameters for --soloCellFilterType EmptyDrops_CR
ParametersSolo: --soloCellFilterType EmptyDrops_CR filtering parameters: 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000
Number of CBs in the whitelist = 737280
Oct 11 16:23:23 ... Finished reading, sorting and deduplicating CB whitelist sequences.
Finished loading and checking parameters
Reading genome generation parameters:
### /woldlab/loxcyc/home/diane/proj/STAR/bin/Linux_x86_64_static/STAR   --runMode genomeGenerate      --runThreadN 12   --genomeDir /woldlab/loxcyc/home/diane/proj/genome/GRCh38-V29-male-2.7.8a   --genomeFastaFiles GRCh38_no_alt_analysis_set_GCA_000001405.15.fa      --sjdbGTFfile gencode.vV29-tRNAs-ERCC.gff   --sjdbOverhang 100
### GstrandBit=32
versionGenome                 2.7.4a     ~RE-DEFINED
genomeType                    Full     ~RE-DEFINED
genomeFastaFiles              GRCh38_no_alt_analysis_set_GCA_000001405.15.fa        ~RE-DEFINED
genomeSAindexNbases           14     ~RE-DEFINED
genomeChrBinNbits             18     ~RE-DEFINED
genomeSAsparseD               1     ~RE-DEFINED
genomeTransformType           None     ~RE-DEFINED
genomeTransformVCF            -     ~RE-DEFINED
sjdbOverhang                  100     ~RE-DEFINED
sjdbFileChrStartEnd           -        ~RE-DEFINED
sjdbGTFfile                   gencode.vV29-tRNAs-ERCC.gff     ~RE-DEFINED
sjdbGTFchrPrefix              -     ~RE-DEFINED
sjdbGTFfeatureExon            exon     ~RE-DEFINED
sjdbGTFtagExonParentTranscripttranscript_id     ~RE-DEFINED
sjdbGTFtagExonParentGene      gene_id     ~RE-DEFINED
sjdbInsertSave                Basic     ~RE-DEFINED
genomeFileSizes               3236730736   24809459115        ~RE-DEFINED
Genome version is compatible with current STAR
Number of real (reference) chromosomes= 292
1	chr1	248956422	0
2	chr2	242193529	249036800
3	chr3	198295559	491257856
4	chr4	190214555	689700864
5	chr5	181538259	880017408
6	chr6	170805979	1061683200
7	chr7	159345973	1232601088
8	chr8	145138636	1391984640
9	chr9	138394717	1537212416
10	chr10	133797422	1675624448
11	chr11	135086622	1809580032
12	chr12	133275309	1944846336
13	chr13	114364328	2078277632
14	chr14	107043718	2192834560
15	chr15	101991189	2300051456
16	chr16	90338345	2402287616
17	chr17	83257441	2492727296
18	chr18	80373285	2576089088
19	chr19	58617616	2656567296
20	chr20	64444167	2715287552
21	chr21	46709983	2779774976
22	chr22	50818468	2826698752
23	chrX	156040895	2877554688
24	chrY	57227415	3033792512
25	chrM	16569	3091202048
26	chr1_KI270706v1_random	175055	3091464192
27	chr1_KI270707v1_random	32032	3091726336
28	chr1_KI270708v1_random	127682	3091988480
29	chr1_KI270709v1_random	66860	3092250624
30	chr1_KI270710v1_random	40176	3092512768
31	chr1_KI270711v1_random	42210	3092774912
32	chr1_KI270712v1_random	176043	3093037056
33	chr1_KI270713v1_random	40745	3093299200
34	chr1_KI270714v1_random	41717	3093561344
35	chr2_KI270715v1_random	161471	3093823488
36	chr2_KI270716v1_random	153799	3094085632
37	chr3_GL000221v1_random	155397	3094347776
38	chr4_GL000008v2_random	209709	3094609920
39	chr5_GL000208v1_random	92689	3094872064
40	chr9_KI270717v1_random	40062	3095134208
41	chr9_KI270718v1_random	38054	3095396352
42	chr9_KI270719v1_random	176845	3095658496
43	chr9_KI270720v1_random	39050	3095920640
44	chr11_KI270721v1_random	100316	3096182784
45	chr14_GL000009v2_random	201709	3096444928
46	chr14_GL000225v1_random	211173	3096707072
47	chr14_KI270722v1_random	194050	3096969216
48	chr14_GL000194v1_random	191469	3097231360
49	chr14_KI270723v1_random	38115	3097493504
50	chr14_KI270724v1_random	39555	3097755648
51	chr14_KI270725v1_random	172810	3098017792
52	chr14_KI270726v1_random	43739	3098279936
53	chr15_KI270727v1_random	448248	3098542080
54	chr16_KI270728v1_random	1872759	3099066368
55	chr17_GL000205v2_random	185591	3101163520
56	chr17_KI270729v1_random	280839	3101425664
57	chr17_KI270730v1_random	112551	3101949952
58	chr22_KI270731v1_random	150754	3102212096
59	chr22_KI270732v1_random	41543	3102474240
60	chr22_KI270733v1_random	179772	3102736384
61	chr22_KI270734v1_random	165050	3102998528
62	chr22_KI270735v1_random	42811	3103260672
63	chr22_KI270736v1_random	181920	3103522816
64	chr22_KI270737v1_random	103838	3103784960
65	chr22_KI270738v1_random	99375	3104047104
66	chr22_KI270739v1_random	73985	3104309248
67	chrY_KI270740v1_random	37240	3104571392
68	chrUn_KI270302v1	2274	3104833536
69	chrUn_KI270304v1	2165	3105095680
70	chrUn_KI270303v1	1942	3105357824
71	chrUn_KI270305v1	1472	3105619968
72	chrUn_KI270322v1	21476	3105882112
73	chrUn_KI270320v1	4416	3106144256
74	chrUn_KI270310v1	1201	3106406400
75	chrUn_KI270316v1	1444	3106668544
76	chrUn_KI270315v1	2276	3106930688
77	chrUn_KI270312v1	998	3107192832
78	chrUn_KI270311v1	12399	3107454976
79	chrUn_KI270317v1	37690	3107717120
80	chrUn_KI270412v1	1179	3107979264
81	chrUn_KI270411v1	2646	3108241408
82	chrUn_KI270414v1	2489	3108503552
83	chrUn_KI270419v1	1029	3108765696
84	chrUn_KI270418v1	2145	3109027840
85	chrUn_KI270420v1	2321	3109289984
86	chrUn_KI270424v1	2140	3109552128
87	chrUn_KI270417v1	2043	3109814272
88	chrUn_KI270422v1	1445	3110076416
89	chrUn_KI270423v1	981	3110338560
90	chrUn_KI270425v1	1884	3110600704
91	chrUn_KI270429v1	1361	3110862848
92	chrUn_KI270442v1	392061	3111124992
93	chrUn_KI270466v1	1233	3111649280
94	chrUn_KI270465v1	1774	3111911424
95	chrUn_KI270467v1	3920	3112173568
96	chrUn_KI270435v1	92983	3112435712
97	chrUn_KI270438v1	112505	3112697856
98	chrUn_KI270468v1	4055	3112960000
99	chrUn_KI270510v1	2415	3113222144
100	chrUn_KI270509v1	2318	3113484288
101	chrUn_KI270518v1	2186	3113746432
102	chrUn_KI270508v1	1951	3114008576
103	chrUn_KI270516v1	1300	3114270720
104	chrUn_KI270512v1	22689	3114532864
105	chrUn_KI270519v1	138126	3114795008
106	chrUn_KI270522v1	5674	3115057152
107	chrUn_KI270511v1	8127	3115319296
108	chrUn_KI270515v1	6361	3115581440
109	chrUn_KI270507v1	5353	3115843584
110	chrUn_KI270517v1	3253	3116105728
111	chrUn_KI270529v1	1899	3116367872
112	chrUn_KI270528v1	2983	3116630016
113	chrUn_KI270530v1	2168	3116892160
114	chrUn_KI270539v1	993	3117154304
115	chrUn_KI270538v1	91309	3117416448
116	chrUn_KI270544v1	1202	3117678592
117	chrUn_KI270548v1	1599	3117940736
118	chrUn_KI270583v1	1400	3118202880
119	chrUn_KI270587v1	2969	3118465024
120	chrUn_KI270580v1	1553	3118727168
121	chrUn_KI270581v1	7046	3118989312
122	chrUn_KI270579v1	31033	3119251456
123	chrUn_KI270589v1	44474	3119513600
124	chrUn_KI270590v1	4685	3119775744
125	chrUn_KI270584v1	4513	3120037888
126	chrUn_KI270582v1	6504	3120300032
127	chrUn_KI270588v1	6158	3120562176
128	chrUn_KI270593v1	3041	3120824320
129	chrUn_KI270591v1	5796	3121086464
130	chrUn_KI270330v1	1652	3121348608
131	chrUn_KI270329v1	1040	3121610752
132	chrUn_KI270334v1	1368	3121872896
133	chrUn_KI270333v1	2699	3122135040
134	chrUn_KI270335v1	1048	3122397184
135	chrUn_KI270338v1	1428	3122659328
136	chrUn_KI270340v1	1428	3122921472
137	chrUn_KI270336v1	1026	3123183616
138	chrUn_KI270337v1	1121	3123445760
139	chrUn_KI270363v1	1803	3123707904
140	chrUn_KI270364v1	2855	3123970048
141	chrUn_KI270362v1	3530	3124232192
142	chrUn_KI270366v1	8320	3124494336
143	chrUn_KI270378v1	1048	3124756480
144	chrUn_KI270379v1	1045	3125018624
145	chrUn_KI270389v1	1298	3125280768
146	chrUn_KI270390v1	2387	3125542912
147	chrUn_KI270387v1	1537	3125805056
148	chrUn_KI270395v1	1143	3126067200
149	chrUn_KI270396v1	1880	3126329344
150	chrUn_KI270388v1	1216	3126591488
151	chrUn_KI270394v1	970	3126853632
152	chrUn_KI270386v1	1788	3127115776
153	chrUn_KI270391v1	1484	3127377920
154	chrUn_KI270383v1	1750	3127640064
155	chrUn_KI270393v1	1308	3127902208
156	chrUn_KI270384v1	1658	3128164352
157	chrUn_KI270392v1	971	3128426496
158	chrUn_KI270381v1	1930	3128688640
159	chrUn_KI270385v1	990	3128950784
160	chrUn_KI270382v1	4215	3129212928
161	chrUn_KI270376v1	1136	3129475072
162	chrUn_KI270374v1	2656	3129737216
163	chrUn_KI270372v1	1650	3129999360
164	chrUn_KI270373v1	1451	3130261504
165	chrUn_KI270375v1	2378	3130523648
166	chrUn_KI270371v1	2805	3130785792
167	chrUn_KI270448v1	7992	3131047936
168	chrUn_KI270521v1	7642	3131310080
169	chrUn_GL000195v1	182896	3131572224
170	chrUn_GL000219v1	179198	3131834368
171	chrUn_GL000220v1	161802	3132096512
172	chrUn_GL000224v1	179693	3132358656
173	chrUn_KI270741v1	157432	3132620800
174	chrUn_GL000226v1	15008	3132882944
175	chrUn_GL000213v1	164239	3133145088
176	chrUn_KI270743v1	210658	3133407232
177	chrUn_KI270744v1	168472	3133669376
178	chrUn_KI270745v1	41891	3133931520
179	chrUn_KI270746v1	66486	3134193664
180	chrUn_KI270747v1	198735	3134455808
181	chrUn_KI270748v1	93321	3134717952
182	chrUn_KI270749v1	158759	3134980096
183	chrUn_KI270750v1	148850	3135242240
184	chrUn_KI270751v1	150742	3135504384
185	chrUn_KI270752v1	27745	3135766528
186	chrUn_KI270753v1	62944	3136028672
187	chrUn_KI270754v1	40191	3136290816
188	chrUn_KI270755v1	36723	3136552960
189	chrUn_KI270756v1	79590	3136815104
190	chrUn_KI270757v1	71251	3137077248
191	chrUn_GL000214v1	137718	3137339392
192	chrUn_KI270742v1	186739	3137601536
193	chrUn_GL000216v2	176608	3137863680
194	chrUn_GL000218v1	161147	3138125824
195	chrEBV	171823	3138387968
196	ERCC-00002	1061	3138650112
197	ERCC-00003	1023	3138912256
198	ERCC-00004	523	3139174400
199	ERCC-00007	1135	3139436544
200	ERCC-00009	984	3139698688
201	ERCC-00012	994	3139960832
202	ERCC-00013	808	3140222976
203	ERCC-00014	1957	3140485120
204	ERCC-00016	844	3140747264
205	ERCC-00017	1136	3141009408
206	ERCC-00018	1032	3141271552
207	ERCC-00019	644	3141533696
208	ERCC-00022	751	3141795840
209	ERCC-00023	273	3142057984
210	ERCC-00024	536	3142320128
211	ERCC-00025	1994	3142582272
212	ERCC-00028	1130	3142844416
213	ERCC-00031	1138	3143106560
214	ERCC-00033	2022	3143368704
215	ERCC-00034	1019	3143630848
216	ERCC-00035	1130	3143892992
217	ERCC-00039	740	3144155136
218	ERCC-00040	744	3144417280
219	ERCC-00041	1122	3144679424
220	ERCC-00042	1023	3144941568
221	ERCC-00043	1023	3145203712
222	ERCC-00044	1156	3145465856
223	ERCC-00046	522	3145728000
224	ERCC-00048	992	3145990144
225	ERCC-00051	274	3146252288
226	ERCC-00053	1023	3146514432
227	ERCC-00054	274	3146776576
228	ERCC-00057	1021	3147038720
229	ERCC-00058	1136	3147300864
230	ERCC-00059	525	3147563008
231	ERCC-00060	523	3147825152
232	ERCC-00061	1136	3148087296
233	ERCC-00062	1023	3148349440
234	ERCC-00067	644	3148611584
235	ERCC-00069	1137	3148873728
236	ERCC-00071	642	3149135872
237	ERCC-00073	603	3149398016
238	ERCC-00074	522	3149660160
239	ERCC-00075	1023	3149922304
240	ERCC-00076	642	3150184448
241	ERCC-00077	273	3150446592
242	ERCC-00078	993	3150708736
243	ERCC-00079	644	3150970880
244	ERCC-00081	534	3151233024
245	ERCC-00083	1022	3151495168
246	ERCC-00084	994	3151757312
247	ERCC-00085	844	3152019456
248	ERCC-00086	1020	3152281600
249	ERCC-00092	1124	3152543744
250	ERCC-00095	521	3152805888
251	ERCC-00096	1107	3153068032
252	ERCC-00097	523	3153330176
253	ERCC-00098	1143	3153592320
254	ERCC-00099	1350	3153854464
255	ERCC-00104	2022	3154116608
256	ERCC-00108	1022	3154378752
257	ERCC-00109	536	3154640896
258	ERCC-00111	994	3154903040
259	ERCC-00112	1136	3155165184
260	ERCC-00113	840	3155427328
261	ERCC-00116	1991	3155689472
262	ERCC-00117	1136	3155951616
263	ERCC-00120	536	3156213760
264	ERCC-00123	1022	3156475904
265	ERCC-00126	1118	3156738048
266	ERCC-00128	1133	3157000192
267	ERCC-00130	1059	3157262336
268	ERCC-00131	771	3157524480
269	ERCC-00134	274	3157786624
270	ERCC-00136	1033	3158048768
271	ERCC-00137	537	3158310912
272	ERCC-00138	1024	3158573056
273	ERCC-00142	493	3158835200
274	ERCC-00143	784	3159097344
275	ERCC-00144	538	3159359488
276	ERCC-00145	1042	3159621632
277	ERCC-00147	1023	3159883776
278	ERCC-00148	494	3160145920
279	ERCC-00150	743	3160408064
280	ERCC-00154	537	3160670208
281	ERCC-00156	494	3160932352
282	ERCC-00157	1019	3161194496
283	ERCC-00158	1027	3161456640
284	ERCC-00160	743	3161718784
285	ERCC-00162	523	3161980928
286	ERCC-00163	543	3162243072
287	ERCC-00164	1022	3162505216
288	ERCC-00165	872	3162767360
289	ERCC-00168	1024	3163029504
290	ERCC-00170	1023	3163291648
291	ERCC-00171	505	3163553792
292	phiX174	5386	3163815936
--sjdbOverhang = 100 taken from the generated genome
Started loading the genome: Mon Oct 11 16:23:23 2021

Genome: size given as a parameter = 3236730736
SA: size given as a parameter = 24809459115
SAindex: size given as a parameter = 1
Read from SAindex: pGe.gSAindexNbases=14  nSAi=357913940
nGenome=3236730736;  nSAbyte=24809459115
GstrandBit=32   SA number of indices=6014414330
Shared memory is not used for genomes. Allocated a private copy of the genome.
Genome file size: 3236730736 bytes; state: good=1 eof=0 fail=0 bad=0
Loading Genome ... done! state: good=1 eof=0 fail=0 bad=0; loaded 3236730736 bytes
SA file size: 24809459115 bytes; state: good=1 eof=0 fail=0 bad=0
Loading SA ... done! state: good=1 eof=0 fail=0 bad=0; loaded 24809459115 bytes
Loading SAindex ... done: 1565873619 bytes
Finished loading the genome: Mon Oct 11 16:25:46 2021

Processing splice junctions database sjdbN=361456,   pGe.sjdbOverhang=100 
To accommodate alignIntronMax=1000000 redefined winBinNbits=18
To accommodate alignIntronMax=1000000 and alignMatesGapMax=1000000, redefined winFlankNbins=4 and winAnchorDistNbins=8
Loaded transcript database, nTr=207507
Loaded exon database, nEx=1263774
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread0 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread0 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread1 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread1 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread2 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread2 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread3 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread3 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread4 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread4 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread5 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread5 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread6 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread6 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread7 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread7 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread8 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread8 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread9 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread9 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread10 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread10 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread11 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread11 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread12 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread12 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread13 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread13 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread14 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread14 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate1.thread15 ... ok
Opening the file: _STARtmp//FilterBySJoutFiles.mate2.thread15 ... ok
Created thread # 1
Created thread # 2
Created thread # 3
Starting to map file # 0
mate 1:   SRX5908538_R2.fastq.gz
mate 2:   SRX5908538_R1.fastq.gz
Created thread # 4
Created thread # 5
Created thread # 6
Created thread # 7
Created thread # 8
Created thread # 9
Created thread # 10
Created thread # 11
Created thread # 12
Created thread # 13
Created thread # 14
Created thread # 15
BAM sorting: 144740 mapped reads
BAM sorting bins genomic start loci:
1	0	32796727
2	0	92835200
3	0	153534777
4	0	225790044
5	1	46385031
6	1	118833911
7	1	231256214
8	2	39412040
9	2	136538517
10	3	83461277
11	4	16766073
12	4	126498376
13	4	181237608
14	5	33276055
15	5	73517728
16	6	20002986
17	6	98385848
18	6	141619428
19	7	98762485
20	8	173322
21	8	124860196
22	9	28899272
23	9	124760061
24	10	61964584
25	10	61964779
26	10	65499642
27	10	119674638
28	11	54284613
29	11	112405239
30	12	28221073
31	13	55067205
32	14	43793787
33	14	69455267
34	15	1962201
35	15	22010682
36	15	89561597
37	16	16441788
38	16	46721607
39	17	2581018
40	18	7561644
41	18	18592171
42	18	48616114
43	18	58393455
44	19	59034582
45	21	35339127
46	22	19567233
47	22	119925180
48	24	7235
Thread #2 end of input stream, nextChar=-1
Completed: thread #13
Completed: thread #9
Completed: thread #5
Completed: thread #8
Completed: thread #15
Completed: thread #4
Completed: thread #14
Completed: thread #0
Completed: thread #2
Completed: thread #1
Joined thread # 1
Joined thread # 2
Completed: thread #6
Completed: thread #12
Completed: thread #11
Completed: thread #3
Joined thread # 3
Joined thread # 4
Joined thread # 5
Joined thread # 6
Completed: thread #10
Completed: thread #7
Joined thread # 7
Joined thread # 8
Joined thread # 9
Joined thread # 10
Joined thread # 11
Joined thread # 12
Joined thread # 13
Joined thread # 14
Joined thread # 15
Completed stage 1 mapping of outFilterBySJout mapping
Detected 81352 novel junctions that passed filtering, will proceed to filter reads that contained unannotated junctions
Created thread # 1
Created thread # 2
Created thread # 3
Created thread # 4
Created thread # 5
Created thread # 6
Created thread # 7
Created thread # 8
Created thread # 9
Created thread # 10
Created thread # 11
Created thread # 12
Created thread # 13
Created thread # 14
Created thread # 15
Completed: thread #6
Completed: thread #11
Completed: thread #1
Completed: thread #3
Completed: thread #2
Completed: thread #7
Completed: thread #5
Completed: thread #13
Completed: thread #10
Completed: thread #14
Completed: thread #12
Completed: thread #8
Completed: thread #0
Joined thread # 1
Joined thread # 2
Joined thread # 3
Completed: thread #4
Joined thread # 4
Joined thread # 5
Joined thread # 6
Joined thread # 7
Joined thread # 8
Completed: thread #15
Completed: thread #9
Joined thread # 9
Joined thread # 10
Joined thread # 11
Joined thread # 12
Joined thread # 13
Joined thread # 14
Joined thread # 15
Oct 11 16:39:42 ..... started Solo counting
Oct 11 16:39:42 ... Starting Solo post-map for SJ
Oct 11 16:39:42 ... Finished allocating arrays for Solo 0.192355 GB
Oct 11 16:39:57 ... Finished reading reads from Solo files nCB=215105, nReadPerCBmax=40248, yesWLmatch=0
Oct 11 16:40:19 ... Finished collapsing UMIs
Oct 11 16:40:19 ... Solo: writing raw matrix
Solo output directory directory exists and will be overwritten: ./Solo.out/SJ//raw/
Oct 11 16:40:21 ... Solo: cell filtering
Oct 11 16:40:21 ... Starting Solo post-map for GeneFull_Ex50pAS
Oct 11 16:40:21 ... Allocated and initialized readInfo array, nReadsInput = 132228862
Oct 11 16:40:21 ... Finished allocating arrays for Solo 1.40314 GB
Oct 11 16:41:26 ... Finished reading reads from Solo files nCB=440936, nReadPerCBmax=195203, yesWLmatch=0
Oct 11 16:43:49 ... Finished collapsing UMIs
Oct 11 16:43:49 ... Solo: writing raw matrix
Solo output directory directory exists and will be overwritten: ./Solo.out/GeneFull_Ex50pAS//raw/
Oct 11 16:44:00 ... Solo: cell filtering
cellFiltering: simple: nUMImax=55568; nUMImin=5557; nCellsSimple=2715
Oct 11 16:44:00 ... starting emptyDrops_CR filtering
Oct 11 16:44:01 ... finished ambient cells counting
Oct 11 16:44:01 ... finished SGT
Oct 11 16:44:01 ... finished ambient profile
Oct 11 16:44:01 ... candidate cells: minUMI=500; number of candidate cells=445
Oct 11 16:44:01 ... finished observed logProb
Oct 11 16:44:09 ... finished simulations
Oct 11 16:44:09 ... finished emptyDrops_CR filtering: number of additional non-ambient cells=244
Solo output directory directory exists and will be overwritten: ./Solo.out/GeneFull_Ex50pAS/filtered/
Oct 11 16:44:12 ..... finished Solo counting
Oct 11 16:44:12 ..... started sorting BAM
Max memory needed for sorting = 1107643494
ALL DONE!

"""


def generate_barcodes():
    letters = "AGTC"
    for barcode in itertools.product(letters, repeat=4):
        yield "".join(barcode)


def filter_barcodes(barcodes):
    for barcode in barcodes:
        if barcode[0] in ("G", "C"):
            yield barcode


def write_barcodes(filename, barcodes):
    with open(filename, "wt") as outstream:
        for barcode in barcodes:
            outstream.write("{}\n".format(barcode))


def generate_features(feature_length):
    for feature_id in range(feature_length):
        yield (
            "GENE{:0>5}".format(feature_id),
            "symbol{:0>5}".format(feature_id),
            "Gene Expression",
        )


def write_features(filename, features):
    with open(filename, "wt") as outstream:
        for gene_id, symbol, expression in features:
            outstream.write("{}\t{}\t{}\n".format(gene_id, symbol, expression))


def generate_count_matrix(barcodes, feature_length):
    counts = []
    for barcode in barcodes:
        seed = [ord(x) for x in barcode]
        random.seed(seed)
        row = numpy.floor(random.gamma(0.5, 10, size=feature_length))
        counts.append(row)
    counts = numpy.asarray(counts).T
    return counts


def write_count_matrix(filename, features, barcodes, counts):
    with open(filename, "wt") as outstream:
        # Write header
        outstream.write("%%MatrixMarket matrix coordinate integer general\n")
        # Write comment
        outstream.write("% This is completely synthetic data\n")
        outstream.write(
            "{} {} {}\n".format(
                len(list(features)),
                len(list(barcodes)),
                numpy.count_nonzero(counts),
            )
        )

        for column in range(counts.shape[1]):
            for row in range(counts.shape[0]):
                if counts[row][column] > 0:
                    outstream.write(
                        "{} {} {}\n".format(
                            row + 1, column + 1, int(counts[row][column])
                        )
                    )


def make_sample_data(scratch_dir):
    experiment_dir = scratch_dir / "experiment"
    solo_dir = experiment_dir / "Solo.out"
    gene_dir = solo_dir / "GeneFull_Ex50pAS"
    filtered_dir = gene_dir / "filtered"
    raw_dir = gene_dir / "raw"

    dirs = [experiment_dir, solo_dir, gene_dir, filtered_dir, raw_dir]
    for d in dirs:
        if not d.exists():
            d.mkdir()

    log_out = experiment_dir / "Log.out"
    with open(log_out, "wt") as outstream:
        outstream.write(STAR_2_7_9a_log)

    raw_barcodes = raw_dir / "barcodes.tsv"
    write_barcodes(raw_barcodes, generate_barcodes())

    raw_features = raw_dir / "features.tsv"
    feature_length = 100
    write_features(raw_features, generate_features(feature_length))

    raw_matrix = raw_dir / "matrix.mtx"
    raw_counts = generate_count_matrix(generate_barcodes(), feature_length)
    write_count_matrix(
        raw_matrix, generate_features(feature_length), generate_barcodes(), raw_counts
    )

    filtered_barcodes = filtered_dir / "barcodes.tsv"
    write_barcodes(filtered_barcodes, filter_barcodes(generate_barcodes()))

    filtered_features = filtered_dir / "features.tsv"
    write_features(filtered_features, generate_features(feature_length))

    filtered_matrix = filtered_dir / "matrix.mtx"
    filtered_counts = generate_count_matrix(
        filter_barcodes(generate_barcodes()), feature_length
    )
    write_count_matrix(
        filtered_matrix,
        generate_features(feature_length),
        filter_barcodes(generate_barcodes()),
        filtered_counts,
    )

    return experiment_dir


@contextmanager
def scratch_dir():
    try:
        scratch_dir = TemporaryDirectory(suffix="_mex")
        yield Path(scratch_dir.name)
    finally:
        scratch_dir.cleanup()


if __name__ == "__main__":
    make_sample_data(Path("/tmp/sample"))
