---
layout: post
read_time: true
show_date: true
title:  cellranger-arc 使用说明及bug fix
date:   2024-01-21 13:40:20 -0600
header-img: img/20240121/千姬.jpg
tags: [计算生物学，生信软件，cellranger，fastp]
author: 孙睿
#github:  amaynez/Perceptron/
mathjax: yes
catalog: true
---

本篇博客主要记录自己最近在使用10X的cellranger系列软件(cellranger, cellranger-atac, cellranger-arc)过程中遇到的一个问题及解决过程，顺带整理下相关的软件使用步骤，方便后续工作查找。我们首先从bug说起。

## 问题概述
遇到的问题是，手头一批10X-multiome数据在使用cellranger-arc count 进行处理时出现了报错，初始报错如下

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="/img/20240121/image-1.png" width = "100%" alt=""/>
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">
    图1：犬蝠海马多组学报错，ATAC-count中出现了过长的R1reads 26-151
  	</div>
</center>

和合作方测序工程师沟通后，反馈结果如下
<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="/img/20240121/image-2.png" width = "80%" alt=""/>
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">
    图2：反馈
  	</div>
</center>

结合上面报错，看起来问题是在测序时使用了更长的测序深度，导致出现长度不一致的片段。

### 初次解决尝试

直接参考图1中10X的日志信息，在cellranger-arc count中设置参数 --r1-length=26, 报错。 原因是cellranger-arc count中根本没有这个参数，致信10X官方后也得到了这个回复，确实没有。淦！

之后考虑分别跑rna-seq和atac-seq，10X 官网上倒是有提怎么把多组学数据分别用RNA和ATAC的流程处理，但是这样做和直接cellranger-arc count之间的差异有多大，需要后续咨询10X官方。

开始先做的cellranger count看RNA-seq效果，此时需要加额外参数 --r1-length=26 和 --chemistry="ARC-v1"。最后结果RNA-seq是可行的，能出结果并且结果不太烂。

### 后续解决尝试

后面决定还是优先考虑用多组学流程去做，毕竟多组学比单组学要贵好多并且质量是要低于单组学的，要是只看RNA，实在是亏死……

#### 2024.01.16

想直接拿一个工具把测序数据中违规的片段过滤掉。尝试的fastp，按照图2中的反馈，设置扔掉片段长度小于50的。（10X报错是26，为啥要扔50的？）


只处理了ATAC-seq数据，用的工具是fastp，但是一个问题是，fastp貌似只支持R1,R2两条序列同时处理，但ATAC是R1,R2,R3,所以没法同时处理，最后决定对每个Read 单独处理，就是

```
fastp -in R1 -out R1_f --r1-length=50  
fastp -in R2 ......    # 这里的命令不是正确的fastp命令，只是一个示意
fastp -in R3 ......
```

用处理完后的ATAC数据扔到cellranger-arc count上操作，继续报错，结果如下

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="/img/20240121/image-4.png" width = "100%" alt=""/>
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">
    fastp单独处理后cellranger-arc 报错
  	</div>
</center>

可以看到问题是说过滤后数据的第八行不匹配，用 zcat file | head -n 20 命令查看R1R2R3后发现，有的片段,R1中过滤掉了，但是R2R3中没有，所以出现的问题。 猜测，要是一个片段R1R2R3同时过滤掉，就没有问题了？

后面又问昆动所工程师，说让我只用 fastp -in file1 -out file2 这样的默认参数试一下（我很怀疑，这样还是分开处理，怎么解决上面的问题？）

果不其然，鸟用没有，就是做了个质控，还是遇到不匹配问题……

### 2024.01.19 

进一步排除各种可能的报错，这里是想知道RNA-seq测序有没有问题，结果还是出问题了。下面是长翼蝠海马和长翼蝠垂体的报错，

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="/img/20240121/image.png" width = "100%" alt=""/>
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">
    长翼蝠海马报错
  	</div>
</center>

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="/img/20240121/image-3.png" width = "100%" alt=""/>
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">
    长翼蝠垂体报错
  	</div>
</center>

这是长翼蝠垂体RNA-seq处理使用的脚本，可以看到对RNA-seq，不设置--r1-length=26时会出现相同的报错。但是和图1的报错有冲突，图1报错的时候是在ATAC的步骤上出错的
```
#!/bin/sh
#SBATCH --output=CYCT.out
#SBATCH --error=CYCT.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sunrui171@mails.ucas.edu.cn

cellranger count --id=CYCT \
                 --sample=CYCT \
                 --chemistry='ARC-v1' \
                 --transcriptome=/temp/swap/sunrui/genome/mfu \
                 --fastqs=/work/swap/sunrui/rna_seq/clean_rna/CY/CY-CT-G \
                 --localcores=48 \
                 --localmem=64

```

这里是想看ATAC-seq数据有没有问题，看起来还是有问题, 报错如下

```
2024-01-19 17:32:33 [runtime] (chunks_complete) ID.CYHM.SC_ATAC_COUNTER_CS.SC_ATAC_COUNTER._BASIC_SC_ATAC_COUNTER._ATAC_MATRIX_COMPUTER.ALIGN_ATAC_READS
2024-01-19 17:32:33 [runtime] (run:local)       ID.CYHM.SC_ATAC_COUNTER_CS.SC_ATAC_COUNTER._BASIC_SC_ATAC_COUNTER._ATAC_MATRIX_COMPUTER.ALIGN_ATAC_READS.fork0.join
2024-01-19 17:32:33 [runtime] (failed)          ID.CYHM.SC_ATAC_COUNTER_CS.SC_ATAC_COUNTER._BASIC_SC_ATAC_COUNTER._ATAC_MATRIX_COMPUTER.ALIGN_ATAC_READS

[error] Pipestance failed. Error log at:
CYHM/SC_ATAC_COUNTER_CS/SC_ATAC_COUNTER/_BASIC_SC_ATAC_COUNTER/_ATAC_MATRIX_COMPUTER/ALIGN_ATAC_READS/fork0/join-uae59aa3032/_errors

Log message:
2.2% (< 10%) of read pairs have a valid 10x barcode. This could be a result of poor sequencing quality, a sample mixup, or running the wrong pipeline, for example, running `cellranger-atac` on Multiome ATAC + GEX data, or vice versa.

Waiting 6 seconds for UI to do final refresh.
Pipestance failed. Use --noexit option to keep UI running after failure.

2024-01-19 17:32:39 Shutting down.

```

运行脚本如下：

```
(base) sunrui@n07:/work/swap/sunrui/atac_count/CY$ cat CYHM.slurm
#!/bin/bash
#SBATCH --output=CYHM_atac.out
#SBATCH --error=CYHM_atac.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sunrui171@mails.ucas.edu.cn
#SBATCH --nodelist=n07

cellranger-atac count --id=CYHM \
                      --reference=/temp/swap/sunrui/genome/mfu \
                      --fastqs=/work/swap/sunrui/atac_seq/clean_data/CY/CYHM \
                     # --sample=mysample \
                      --localcores=48 \
                      --localmem=64 \
                      --chemistry='ARC-v1'

```

上面的ATAC报错居然没有先报错说ATAC这边Read有问题，而是说有效barcode太少。这种报错我感觉更像是chemistry这个参数没有设置合适导致的。我又看了下cellranger-atac count 的帮助，发现确实没有--chemistry这个参数，但是这个用法10X官网教程上有，无语，回头必须质问10X。

到这一步重新猜测，是不是ATAC是OK的，我只要把RNA截一下就能在cellranger-arc count上跑了？

截止到目前想到两个需要进一步尝试的策略

- [ ] 用截断后的RNA-seq加原始的ATAC-seq走流程
- [ ] RNA-seq截断， ATAC也截断， ATAC的截断得自己写脚本啊啊啊啊啊啊啊啊啊啊，这样再不成只能宣布或者只看RNA或者让合作者再找公司解决数据问题了，淦！
  
#### 尝试解决 

fastp 同时处理长翼蝠垂体RNA-seq双端数据, fastp结果如下（fastp双端运行中途没有任何输出，安心等15min再说）
```
(base) sunrui@n07:/work/swap/sunrui/rna_seq/clean_rna/CY/CY-CT-G$ fastp -i CYCT_S1_L001_R1_001.fastq.gz -I CYCT_S1_L001_R2_001.fastq.gz -o fCYCT_S1_L001_R1_001.fastq.gz -O fCYCT_S1_L001_R2_001.fastq.gz -b 26 -w 16
Read1 before filtering:
total reads: 406790200
total bases: 61361573650
Q20 bases: 46685950907(76.0834%)
Q30 bases: 37543113863(61.1834%)

Read2 before filtering:
total reads: 406790200
total bases: 61360777858
Q20 bases: 58763295092(95.7669%)
Q30 bases: 54787082253(89.2868%)

Read1 after filtering:
total reads: 404866206
total bases: 10526521033
Q20 bases: 10321364264(98.051%)
Q30 bases: 9909249652(94.136%)

Read2 after filtering:
total reads: 404866206
total bases: 10526486597
Q20 bases: 10292391829(97.7761%)
Q30 bases: 9817541273(93.2651%)

Filtering result:
reads passed filter: 809732412
reads failed due to low quality: 3824906
reads failed due to too many N: 0
reads failed due to too short: 23082
reads with adapter trimmed: 4
bases trimmed due to adapters: 194

Duplication rate: 5.28522%

Insert size peak (evaluated by paired-end reads): 271

JSON report: fastp.json
HTML report: fastp.html

fastp -i CYCT_S1_L001_R1_001.fastq.gz -I CYCT_S1_L001_R2_001.fastq.gz -o fCYCT_S1_L001_R1_001.fastq.gz -O fCYCT_S1_L001_R2_001.fastq.gz -b 26 -w 16
fastp v0.23.4, time used: 679 seconds

```

应用cellranger-arc count 处理应用fastp 截断 R1,R2的RNA-seq(ATAC未处理)，结果如下：

```
2024-01-19 21:30:45 [runtime] (split_complete)  ID.CYCT.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._GEX_MATRIX_COMPUT ER.MAKE_SHARD
2024-01-19 21:30:45 [runtime] (run:local)       ID.CYCT.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._GEX_MATRIX_COMPUT ER.MAKE_SHARD.fork0.chnk0.main
2024-01-19 21:31:12 [runtime] (failed)          ID.CYCT.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._GEX_MATRIX_COMPUT ER.MAKE_SHARD

[error] Pipestance failed. Error log at:
CYCT/SC_ATAC_GEX_COUNTER_CS/SC_ATAC_GEX_COUNTER/_GEX_MATRIX_COMPUTER/MAKE_SHARD/fork0/chnk0-ub9efaa7985/_errors

Log message:
IO error in FASTQ file '"/work/swap/sunrui/rna_seq/clean_rna/CY/fCY-CT-G/fCYCT_S1_L001_R2_001.fastq.gz"', line: 85795 00: unexpected end of file

Waiting 6 seconds for UI to do final refresh.
Pipestance failed. Use --noexit option to keep UI running after failure.

2024-01-19 21:31:18 Shutting down.
^Z
[2]+  Stopped                 tail -f CYCT.out
```

这里发现之前自己一个判断错误，最初的报错是error log 保存在sc_ATAC_GEX文件夹下，但是实际出问题的还是RNA数据，所以只要处理好RNA就够了。之前一直理解错了，尴尬。

看一下出问题的部分
```
(base) sunrui@n07:/work/swap/sunrui/rna_seq/clean_rna/CY/fCY-CT-G$ zcat fCYCT_S1_L001_R2_001.fastq.gz | sed -n '8579400,857960p'

gzip: fCYCT_S1_L001_R2_001.fastq.gz: unexpected end of file
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFF:FFFFFFFFFF
```

不知道啥问题，试一下截断到50？再试一下啥也不做？

这里补一个长翼蝠海马的RNA运行结果, 脚本命令和结果如下 

```
#!/bin/sh
#SBATCH --output=CYHM.out
#SBATCH --error=CYHM.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sunrui171@mails.ucas.edu.cn

cellranger count --id=CYHM \
                 --sample=CYHM \
                 --chemistry='ARC-v1' \
                 --r1-length=26 \
                 --transcriptome=/temp/swap/sunrui/genome/mfu \
                 --fastqs=/work/swap/sunrui/rna_seq/clean_rna/CY/CY-HM-G \
                 --localcores=48 \
                 --localmem=64

```

```
2024-01-20 01:13:03 [runtime] (ready)           ID.CYHM.SC_RNA_COUNTER_CS.SC_MULTI_CORE.MULTI_REPORTER.CHOOSE_CLOUPE
2024-01-20 01:13:03 [runtime] (run:local)       ID.CYHM.SC_RNA_COUNTER_CS.SC_MULTI_CORE.MULTI_REPORTER.CHOOSE_CLOUPE .fork0.chnk0.main
2024-01-20 01:13:04 [runtime] (chunks_complete) ID.CYHM.SC_RNA_COUNTER_CS.SC_MULTI_CORE.MULTI_REPORTER.CHOOSE_CLOUPE
2024-01-20 01:13:45 [runtime] (join_complete)   ID.CYHM.SC_RNA_COUNTER_CS.SC_MULTI_CORE.MULTI_GEM_WELL_PROCESSOR.COU NT_GEM_WELL_PROCESSOR._BASIC_SC_RNA_COUNTER.WRITE_POS_BAM

Outputs:
- Run summary HTML:                         /work/swap/sunrui/rna_count/CY/CYHM/outs/web_summary.html
- Run summary CSV:                          /work/swap/sunrui/rna_count/CY/CYHM/outs/metrics_summary.csv
- BAM:                                      /work/swap/sunrui/rna_count/CY/CYHM/outs/possorted_genome_bam.bam
- BAM BAI index:                            /work/swap/sunrui/rna_count/CY/CYHM/outs/possorted_genome_bam.bam.bai
- BAM CSI index:                            null
- Filtered feature-barcode matrices MEX:    /work/swap/sunrui/rna_count/CY/CYHM/outs/filtered_feature_bc_matrix
- Filtered feature-barcode matrices HDF5:   /work/swap/sunrui/rna_count/CY/CYHM/outs/filtered_feature_bc_matrix.h5
- Unfiltered feature-barcode matrices MEX:  /work/swap/sunrui/rna_count/CY/CYHM/outs/raw_feature_bc_matrix
- Unfiltered feature-barcode matrices HDF5: /work/swap/sunrui/rna_count/CY/CYHM/outs/raw_feature_bc_matrix.h5
- Secondary analysis output CSV:            /work/swap/sunrui/rna_count/CY/CYHM/outs/analysis
- Per-molecule read information:            /work/swap/sunrui/rna_count/CY/CYHM/outs/molecule_info.h5
- CRISPR-specific analysis:                 null
- Antibody aggregate barcodes:              null
- Loupe Browser file:                       /work/swap/sunrui/rna_count/CY/CYHM/outs/cloupe.cloupe
- Feature Reference:                        null
- Target Panel File:                        null
- Probe Set File:                           null

Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!

2024-01-20 01:14:15 Shutting down.
Saving pipestance info to "CYHM/CYHM.mri.tgz"
^Z
[7]+  Stopped                 tail -f CYHM.out
```

可以看到在只设置截断r1为26的情形下，RNA是能够正常运行的。

下面这个实验结果是长翼蝠垂体的多组学数据，用的RNA-seq是截断到51的长度，运行结果如下：

```
2024-01-20 01:12:36 [runtime] (update)          ID.CYCT.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._GEX_MATRIX_COMPUTR.MAKE_SHARD.fork0 join_running
2024-01-20 01:12:43 [runtime] (failed)          ID.CYCT.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._GEX_MATRIX_COMPUTR.MAKE_SHARD

[error] Pipestance failed. Error log at:
CYCT/SC_ATAC_GEX_COUNTER_CS/SC_ATAC_GEX_COUNTER/_GEX_MATRIX_COMPUTER/MAKE_SHARD/fork0/join-ub9efaa7983/_errors

Log message:
ERROR: We detected a mixture of different R1 lengths ([26-51]), which breaks assumptions in how UMIs are tabulated an corrected. To process these data, you will need to truncate to the shortest observed R1 length by providing 26 to th --r1-length argument if are running count/vdj, or via the r1-length parameter in each of the following tables of you multi config CSV if you are running multi: [gene-expression]

Waiting 6 seconds for UI to do final refresh.
Pipestance failed. Use --noexit option to keep UI running after failure.

2024-01-20 01:12:49 Shutting down.

```

从这里我们基本确定问题所在了，在RNA数据中，混入了一小部分长度不等于151（长26）的片段，只要把这些长度不等于151的片段过滤掉，程序应该就能正常运行了。

### 2024.01.20 

应用fastp过滤掉长度不等于151的片段，命令如下

```
(base) sunrui@n06:/work/swap/sunrui/rna_seq/clean_rna/CY/CY-HM-G$ fastp -i CYHM_S1_L001_R1_001.fastq.gz \
      -I CYHM_S1_L001_R2_001.fastq.gz \
      -o /work/swap/sunrui/rna_seq/clean_rna/CY/fCY-HM-G/fCYHM_S1_L001_R1_001.fastq.gz \
      -O /work/swap/sunrui/rna_seq/clean_rna/CY/fCY-HM-G/fCYHM_S1_L001_R2_001.fastq.gz \
      --length_required=150 \
      --thread=16
    Read1 before filtering:
total reads: 453120558
total bases: 68336384444
Q20 bases: 53504228535(78.2954%)
Q30 bases: 45313963764(66.3102%)

Read2 before filtering:
total reads: 453120558
total bases: 68342512478
Q20 bases: 66114787475(96.7404%)
Q30 bases: 62366018303(91.2551%)

Read1 after filtering:
total reads: 449806960
total bases: 67920823659
Q20 bases: 53176824991(78.2924%)
Q30 bases: 45031841437(66.3005%)

Read2 after filtering:
total reads: 449806960
total bases: 67920823659
Q20 bases: 65749637730(96.8034%)
Q30 bases: 62042156342(91.3448%)

Filtering result:
reads passed filter: 899613920
reads failed due to low quality: 2088
reads failed due to too many N: 0
reads failed due to too short: 6625108
reads with adapter trimmed: 0
bases trimmed due to adapters: 0

Duplication rate: 6.15658%

Insert size peak (evaluated by paired-end reads): 271

JSON report: fastp.json
HTML report: fastp.html

fastp -i CYHM_S1_L001_R1_001.fastq.gz -I CYHM_S1_L001_R2_001.fastq.gz -o /work/swap/sunrui/rna_seq/clean_rna/CY/fCY-HM-G/fCYHM_S1_L001_R1_001.fastq.gz -O /work/swap/sunrui/rna_seq/clean_rna/CY/fCY-HM-G/fCYHM_S1_L001_R2_001.fastq.gz --length_required=150 --thread=16
fastp v0.23.4, time used: 1959 seconds
```

<font color=red>**问题成功解决!!!**</font> 

```
- ATAC peak locations:                           /work/swap/sunrui/rna_seq/clean_rna/CY/CYHM/outs/atac_peaks.bed
- ATAC smoothed transposition site track:        /work/swap/sunrui/rna_seq/clean_rna/CY/CYHM/outs/atac_cut_sites.bigwig
- ATAC peak annotations based on proximal genes: /work/swap/sunrui/rna_seq/clean_rna/CY/CYHM/outs/atac_peak_annotation.tsv

Waiting 6 seconds for UI to do final refresh.
Pipestance completed successfully!

2024-01-20 23:32:21 Shutting down.

```