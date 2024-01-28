---
layout: post
read_time: true
show_date: true
title:  cellranger-arc bug fix(2)
date:   2024-01-28 13:40:20 -0600
header-img: img/20240128/茨木童子.jpg
tags: [计算生物学,生信软件,cellranger]
author: 孙睿
#github:  amaynez/Perceptron/
mathjax: yes
catalog: true
---

记录一下这周处理犬蝠10X-multi-omics数据的一个报错。报错如下

```
2024-01-26 03:36:31 [runtime] (run:local)       ID.QFPC.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._SC_ATAC_GEX_ANALYZER._GEX_CLUSTERING_COMPUTER.ATAC_RUN_GR                                                                                APH_CLUSTERING.fork0.join
2024-01-26 03:36:31 [runtime] (failed)          ID.QFPC.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._SC_ATAC_GEX_ANALYZER._PEAK_ANNOTATOR.ANNOTATE_PEAKS
2024-01-26 03:36:31 [runtime] (split_complete)  ID.QFPC.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._SC_ATAC_GEX_ANALYZER._PEAK_ANNOTATOR.GENERATE_TF_MATRIX
2024-01-26 03:36:31 [runtime] (run:local)       ID.QFPC.SC_ATAC_GEX_COUNTER_CS.SC_ATAC_GEX_COUNTER._SC_ATAC_GEX_ANALYZER._PEAK_ANNOTATOR.GENERATE_TF_MATRIX.f                                                                                ork0.join

[error] Pipestance failed. Error log at:
QFPC/SC_ATAC_GEX_COUNTER_CS/SC_ATAC_GEX_COUNTER/_SC_ATAC_GEX_ANALYZER/_PEAK_ANNOTATOR/ANNOTATE_PEAKS/fork0/join-ue55fb2b7d9/_errors

Log message:
Traceback (most recent call last):
  File "/home/sunrui/software/cellranger-arc-2.0.2/external/martian/adapters/python/martian_shell.py", line 659, in _main
    stage.main()
  File "/home/sunrui/software/cellranger-arc-2.0.2/external/martian/adapters/python/martian_shell.py", line 629, in main
    lambda: self._module.join(args, outs, chunk_defs, chunk_outs)
  File "/home/sunrui/software/cellranger-arc-2.0.2/external/martian/adapters/python/martian_shell.py", line 589, in _run
    cmd()
  File "/home/sunrui/software/cellranger-arc-2.0.2/external/martian/adapters/python/martian_shell.py", line 629, in <lambda>
    lambda: self._module.join(args, outs, chunk_defs, chunk_outs)
  File "/home/sunrui/software/cellranger-arc-2.0.2/mro/atac/stages/analysis/annotate_peaks/__init__.py", line 105, in join
    gene_id_name_map = build_gene_id_name_map(ref_mgr)
  File "/home/sunrui/software/cellranger-arc-2.0.2/mro/atac/stages/analysis/annotate_peaks/__init__.py", line 60, in build_gene_id_name_map
    fields.attrs.get("gene_name", fields.attrs["gene_id"])
  File "pybedtools/cbedtools.pyx", line 392, in pybedtools.cbedtools.Interval.attrs.__get__
  File "pybedtools/cbedtools.pyx", line 180, in pybedtools.cbedtools.Attributes.__init__
ValueError: need more than 1 value to unpack


Waiting 6 seconds for UI to do final refresh.
Pipestance failed. Use --noexit option to keep UI running after failure.
```

从上面看，是pybedtools运行的问题，个人感觉问题应该在gtf文件上，因为和犬蝠同一批次的菊头蝠、长翼蝠在上周处理完片段长度后都能顺利运行了，感觉fastq文件应该不太可能有问题了。

后面在cellranger的官方github上找到了一个相关issue https://github.com/10XGenomics/cellranger/issues/220， 给出的一些解决措施如下：

```
We have saw this error ValueError: need more than 1 value to unpack before, this error is usually related with GTF format.

Some formatting tips to use:

- retain only gene_id & transcript_id attributes (#can include gene_names if needed)
- replace/remove the empty gene_id entries
- duplicate transcripts_ids for multiple gene_ids have to be converted as unique
- Remove semi colons after gene_name and avoid bloated col9 in the GTF (extraneous information in this col can cause parsing issues)
More details can be found on our support site: https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/advanced/references.
```
二三两点自己检验了一下后没发现问题，第四点看不太懂，第一条倒是很明确，把第九列中的除了gene_id 和 transcript_id的信息都干掉，这个也容易实现,代码如下：

```
dic = {}

with open("csp.gtf", "r") as source_file, open("simple_csp.gtf", "w") as target_file:
    count = 0
    for line in source_file:
        line_info = line.split('\t')
        info_list = line_info[-1].split(' ')
        info_list = info_list[:4]
        line_info[-1] = ' '.join(info_list) + '\n'
        new_line = '\t'.join(line_info)

        target_file.write(new_line)
```

下面是处理前的犬蝠的gtf文件内容, gene_id 和 transcript_id 是第九列的前四个token，所以上面的代码中只保留前四个token，当然更严格的还是写正则表达式匹配比较好，防止有格式不完全一致的部分。

```
##gtf-version 3
ChrX    EVM     gene    47208   47795   .       -       .       gene_id "Csp.XG.00001"; ID "Csp.XG.00001"; Name "Csp.XG.00001";
ChrX    EVM     transcript      47208   47795   .       -       .       gene_id "Csp.XG.00001"; transcript_id "Csp.XG.00001.t1"; ID "Csp.XG.00001.t1"; Name "Csp.XG.00001"; Parent "Csp.XG.00001"; original_biotype "mrna";
ChrX    EVM     exon    47208   47795   .       -       .       gene_id "Csp.XG.00001"; transcript_id "Csp.XG.00001.t1"; ID "Csp.XG.00001.t1.exon1"; Parent "Csp.XG.00001.t1";
ChrX    EVM     CDS     47208   47795   .       -       0       gene_id "Csp.XG.00001"; transcript_id "Csp.XG.00001.t1"; ID "Csp.XG.00001.t1.cds.1"; Parent "Csp.XG.00001.t1";
ChrX    .       gene    82582   112937  .       -       .       gene_id "Csp.XG.00002"; ID "Csp.XG.00002"; Name "Csp.XG.00002";
ChrX    .       transcript      82582   112937  .       -       .       gene_id "Csp.XG.00002"; transcript_id "Csp.XG.00002.t1"; ID "Csp.XG.00002.t1"; Name "Csp.XG.00002"; Parent "Csp.XG.00002"; original_biotype "mrna";
ChrX    .       exon    82582   83926   .       -       .       gene_id "Csp.XG.00002"; transcript_id "Csp.XG.00002.t1"; ID "Csp.XG.00002.t1.exon4"; Parent "Csp.XG.00002.t1";
ChrX    .       exon    91871   92100   .       -       .       gene_id "Csp.XG.00002"; transcript_id "Csp.XG.00002.t1"; ID "Csp.XG.00002.t1.exon3"; Parent "Csp.XG.00002.t1";
ChrX    .       exon    100470  100688  .       -       .       gene_id "Csp.XG.00002"; transcript_id "Csp.XG.00002.t1"; ID "Csp.XG.00002.t1.exon2"; Parent "Csp.XG.00002.t1";

```

下面是处理完之后的的gtf文件
```
##gtf-version 3

ChrX    EVM     gene    47208   47795   .       -       .       gene_id "Csp.XG.00001"; ID "Csp.XG.00001";
ChrX    EVM     transcript      47208   47795   .       -       .       gene_id "Csp.XG.00001"; transcript_id "Csp.XG.00001.t1";
ChrX    EVM     exon    47208   47795   .       -       .       gene_id "Csp.XG.00001"; transcript_id "Csp.XG.00001.t1";
ChrX    EVM     CDS     47208   47795   .       -       0       gene_id "Csp.XG.00001"; transcript_id "Csp.XG.00001.t1";
ChrX    .       gene    82582   112937  .       -       .       gene_id "Csp.XG.00002"; ID "Csp.XG.00002";
ChrX    .       transcript      82582   112937  .       -       .       gene_id "Csp.XG.00002"; transcript_id "Csp.XG.00002.t1";
ChrX    .       exon    82582   83926   .       -       .       gene_id "Csp.XG.00002"; transcript_id "Csp.XG.00002.t1";
ChrX    .       exon    91871   92100   .       -       .       gene_id "Csp.XG.00002"; transcript_id "Csp.XG.00002.t1";
```

之后用处理后的gtf文件去做参考基因组，然后cellranger-arc count，问题成功解决