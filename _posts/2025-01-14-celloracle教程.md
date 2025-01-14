---
layout: post
read_time: true
show_date: true
title:  生信笔记-2(celloracle 推断基因调控网络)
date:   2025-01-14 13:40:20 -0600
header-img: img/20241202/布都御魂.png
tags: [计算生物学]
author: 孙睿
mathjax: yes
catalog: true
--- 

这篇Blog介绍celloracle原理以及怎样使用celloracle推断基因调控网络。

基因调控网络(GRN)是解读基因互作的重要工具，基于单细胞单组学或者多组学数据的基因调控网络推断也是前几年的一个热门研究方向。celloracle是2023年发表在Nature上的一篇论文，使用单细胞多组学数据进行GRN推断并进行扰动分析。

对于多组学GRN推断，分析上的一个不便之处是很少有一个完备的分析库，要在python包和R包之间来回跳。celloracle在易用性上做了很多的贡献，几乎完全可以依赖这个python库实现网络推断(nice job)。然而官网教程并没有提供一个完全的基于python的GRN推断教程，因此我决定记录下自己使用全python的一个GRN推断流程。

先简单介绍一些celloracle原理。

## celloracle简介


- Input: basic GRN, scRNA-seq data 

    这里basic GRN 要求是一个TF-TG的网络邻接矩阵，basic GRN是有向无权图。对于多组学数据，celloracle官方建议使用scATAC-seq数据构建basic GRN。如果没有scATAC-seq数据，也可以使用一些默认的basic GRN结果。

    basic GRN的作用是使用TF、TG结合的先验信息来对GRN提供约束。假设有N个基因，不加约束的GRN将会是一个$N\times N$大小的邻接矩阵，在推断时存在ill-condition问题。但在使用TF、TG结合的先验信息后，网络大小就回减少到先验TF-TG调控数目上，这是ill-condition问题就会缓解。

    scRNA-seq数据是用来推断GRN网络权重的。

- GRN inference:

    对basic GRN中的每个TG， 其对应的n个细胞表达量为 $x_{tg} \in R^{n}$, 根据basic GRN, 确定调控它的TF集合$\{x_{tf}^{i}\}_{i=1}^{p_{tg}}$。之后使用ridge regression来预测TG $x_{tg}$。

    $$
    \begin{equation}
    x_{tg} =  \sum_{i=1}^{p_{tg}} \beta_i x_{tf}^{i} + \lambda ||\beta||_2^2
    \end{equation}
    $$

    $\beta_i$是TG $x_{tg}$ 和TF $x_{tf}^{i}$的权重，$\lambda$是ridge regression的参数。

- perturbation simulation 

    celloracle 的 perturbation simulation 用一个类似于概率转移的策略实现。举例来说，扰动某个TF的表达量，根据GRN的权重，将扰动分配到TF下游的TG上，之后TG会有扰动，再根据GRN 进行传播。N次传播的扰动变化为

    $$
    \begin{equation}
    \delta^{N} = G^{N-1}\delta, \;\; \delta \in R^{p}
    \end{equation}
    $$

    其中$\delta$是模拟扰动的某个细胞的基因表达量。有了上述扰动向量，可以用类似于single cell velocity的策略来展示扰动后细胞状态的变化。

总的来说，是一个清晰的、可解释的建模策略。

## celloracle使用

celloracle 使用分为两步， base GRN 构建和 GRN inference。

### base GRN 构建 

base GRN的构建可以分为三个步骤， TF-peak 关联构建， TG-peak 关联构建， TF-TG 关联构建。其中TF-peak 关联通常是根据TF结合motif的策略构建，TG-peak关联通常根据TG-peak距离以及相关性构建，TF-TG关联构建则是将前面的两个关联merge到一起，某个TF调控TG必须要有一个公共peak。

#### TG-peak 

celloracle给出的教程中，TG-peak的关联构建使用的是cicero，这需要用户去切换R包。注意到10X cellranger-arc运行结果中，在outs/analysis/feature_linkage/ 目录下给出了feature_linkage.bedpe文件，这个文件的内容格式如下：

```
chr1    3113217 3114070 chr1    3671497 3671498 <_intergenic><Xkr4>     0.445   .       .       6.6603  557854  peak-gene
chr1    3119411 3120213 chr1    3671497 3671498 <_intergenic><Xkr4>     0.4388  .       .       9.3434  551685  peak-gene
chr1    3198065 3199007 chr1    3671497 3671498 <_intergenic><Xkr4>     0.4309  .       .       5.1168  472961  peak-gene
chr1    3216835 3217748 chr1    3671497 3671498 <Xkr4_distal><Xkr4>     0.3979  .       .       5.191   454206  peak-gene
chr1    3299386 3300261 chr1    3671497 3671498 <Xkr4_distal;Gm1992_distal><Xkr4>       0.4404  .       .       8.3617  371674  peak-gene
chr1    3309660 3310543 chr1    3671497 3671498 <Gm1992_distal;Xkr4_distal><Xkr4>       0.4146  .       .       6.4114  361396  peak-gene
chr1    3399529 3400411 chr1    3671497 3671498 <Xkr4_distal;Gm1992_distal><Xkr4>       0.5397  .       .       12.3441 271527  peak-gene
chr1    3407103 3408014 chr1    3671497 3671498 <Gm1992_distal;Xkr4_distal><Xkr4>       0.4074  .       .       8.0994  263939  peak-gene
chr1    3433632 3434520 chr1    3671497 3671498 <Gm1992_distal;Xkr4_distal><Xkr4>       0.4199  .       .       7.9354  237421  peak-gene
chr1    3473876 3474722 chr1    3671497 3671498 <Xkr4_distal;Gm1992_distal><Xkr4>       0.4417  .       .       5.0917  197198  peak-gene
```

可以看到里面提供了TG_peak关联的信息，因此可以直接使用这个文件构建TG-peak关联。下面是我使用的提取函数：

```python
###############################
# generate the relationship between peak(re) and target gene from 10X peak_annot_file and linkage_file
# this function is modified based on multivelo auxilary.py prepare_gene_mat function

def get_re_tg_10x(peak_annot_file, linkage_file, 
                 peak_dist = 10000, min_corr = 0.5,
                 gene_body = False, return_dict=False,
                 parallel=False, n_jobs=1):

    """Peak to gene aggregation.

    This function aggregates promoter and enhancer peaks to genes based on the
    10X linkage file.

    Parameters
    ----------
    adata_atac: :class:`~anndata.AnnData`
        ATAC anndata object which stores raw peak counts.
    peak_annot_file: `str`
        Peak annotation file from 10X CellRanger ARC.
    linkage_file: `str`
        Peak-gene linkage file from 10X CellRanger ARC. This file stores highly
        correlated peak-peak and peak-gene pair information.
    peak_dist: `int` (default: 10000)
        Maximum distance for peaks to be included for a gene.
    min_corr: `float` (default: 0.5)
        Minimum correlation for a peak to be considered as enhancer.
    gene_body: `bool` (default: `False`)
        Whether to add gene body peaks to the associated promoters.
    return_dict: `bool` (default: `False`)
        Whether to return promoter and enhancer dictionaries.

    Returns
    -------
    A new ATAC anndata object which stores gene aggreagted peak counts.
    Additionally, if `return_dict==True`:
        A dictionary which stores genes and promoter peaks.
        And a dictionary which stores genes and enhancer peaks.
    """
    promoter_dict = {}
    distal_dict = {}
    gene_body_dict = {}
    corr_dict = {}

    with open(peak_annot_file) as f:
        header = next(f)
        tmp = header.split('\t')
        if len(tmp) == 4:
            cellranger_version = 1
        elif len(tmp) == 6:
            cellranger_version = 2
        else:
            raise ValueError('Peak annotation file should contain 4 columns '
                             '(CellRanger ARC 1.0.0) or 6 columns (CellRanger '
                             'ARC 2.0.0)')

        #logg.update(f'CellRanger ARC identified as {cellranger_version}.0.0',
        #            v=1)

        if cellranger_version == 1:
            for line in f:
                tmp = line.rstrip().split('\t')
                tmp1 = tmp[0].split('_')
                peak = f'{tmp1[0]}:{tmp1[1]}-{tmp1[2]}'
                if tmp[1] != '':
                    genes = tmp[1].split(';')
                    dists = tmp[2].split(';')
                    types = tmp[3].split(';')
                    for i, gene in enumerate(genes):
                        dist = dists[i]
                        annot = types[i]
                        if annot == 'promoter':
                            if gene not in promoter_dict:
                                promoter_dict[gene] = [peak]
                            else:
                                promoter_dict[gene].append(peak)
                        elif annot == 'distal':
                            if dist == '0':
                                if gene not in gene_body_dict:
                                    gene_body_dict[gene] = [peak]
                                else:
                                    gene_body_dict[gene].append(peak)
                            else:
                                if gene not in distal_dict:
                                    distal_dict[gene] = [peak]
                                else:
                                    distal_dict[gene].append(peak)
        else:
            for line in f:
                tmp = line.rstrip().split('\t')
                peak = f'{tmp[0]}:{tmp[1]}-{tmp[2]}'
                gene = tmp[3]
                dist = tmp[4]
                annot = tmp[5]
                if annot == 'promoter':
                    if gene not in promoter_dict:
                        promoter_dict[gene] = [peak]
                    else:
                        promoter_dict[gene].append(peak)
                elif annot == 'distal':
                    if dist == '0':
                        if gene not in gene_body_dict:
                            gene_body_dict[gene] = [peak]
                        else:
                            gene_body_dict[gene].append(peak)
                    else:
                        if gene not in distal_dict:
                            distal_dict[gene] = [peak]
                        else:
                            distal_dict[gene].append(peak)
    all_info = []
    # read linkages
    with open(linkage_file) as f:
        for line in f:
            tmp = line.rstrip().split('\t')
            all_info.append(tmp)
            if tmp[12] == "peak-peak":
                peak1 = f'{tmp[0]}:{tmp[1]}-{tmp[2]}'
                peak2 = f'{tmp[3]}:{tmp[4]}-{tmp[5]}'
                tmp2 = tmp[6].split('><')[0][1:].split(';')
                tmp3 = tmp[6].split('><')[1][:-1].split(';')
                corr = float(tmp[7])
                for t2 in tmp2:
                    gene1 = t2.split('_')
                    for t3 in tmp3:
                        gene2 = t3.split('_')
                        # one of the peaks is in promoter, peaks belong to the
                        # same gene or are close in distance
                        if (((gene1[1] == "promoter") !=
                            (gene2[1] == "promoter")) and
                            ((gene1[0] == gene2[0]) or
                             (float(tmp[11]) < peak_dist))):

                            if gene1[1] == "promoter":
                                gene = gene1[0]
                            else:
                                gene = gene2[0]
                            if gene in corr_dict:
                                # peak 1 is in promoter, peak 2 is not in gene
                                # body -> peak 2 is added to gene 1
                                if (peak2 not in corr_dict[gene] and
                                    gene1[1] == "promoter" and
                                    (gene2[0] not in gene_body_dict or
                                     peak2 not in gene_body_dict[gene2[0]])):

                                    corr_dict[gene][0].append(peak2)
                                    corr_dict[gene][1].append(corr)
                                # peak 2 is in promoter, peak 1 is not in gene
                                # body -> peak 1 is added to gene 2
                                if (peak1 not in corr_dict[gene] and
                                    gene2[1] == "promoter" and
                                    (gene1[0] not in gene_body_dict or
                                     peak1 not in gene_body_dict[gene1[0]])):

                                    corr_dict[gene][0].append(peak1)
                                    corr_dict[gene][1].append(corr)
                            else:
                                # peak 1 is in promoter, peak 2 is not in gene
                                # body -> peak 2 is added to gene 1
                                if (gene1[1] == "promoter" and
                                    (gene2[0] not in
                                     gene_body_dict
                                     or peak2 not in
                                     gene_body_dict[gene2[0]])):

                                    corr_dict[gene] = [[peak2], [corr]]
                                # peak 2 is in promoter, peak 1 is not in gene
                                # body -> peak 1 is added to gene 2
                                if (gene2[1] == "promoter" and
                                    (gene1[0] not in
                                     gene_body_dict
                                     or peak1 not in
                                     gene_body_dict[gene1[0]])):

                                    corr_dict[gene] = [[peak1], [corr]]
            elif tmp[12] == "peak-gene":
                peak1 = f'{tmp[0]}:{tmp[1]}-{tmp[2]}'
                tmp2 = tmp[6].split('><')[0][1:].split(';')
                gene2 = tmp[6].split('><')[1][:-1]
                corr = float(tmp[7])
                for t2 in tmp2:
                    gene1 = t2.split('_')
                    # peak 1 belongs to gene 2 or are close in distance
                    # -> peak 1 is added to gene 2
                    if ((gene1[0] == gene2) or (float(tmp[11]) < peak_dist)):
                        gene = gene1[0]
                        if gene in corr_dict:
                            if (peak1 not in corr_dict[gene] and
                                gene1[1] != "promoter" and
                                (gene1[0] not in gene_body_dict or
                                 peak1 not in gene_body_dict[gene1[0]])):

                                corr_dict[gene][0].append(peak1)
                                corr_dict[gene][1].append(corr)
                        else:
                            if (gene1[1] != "promoter" and
                                (gene1[0] not in gene_body_dict or
                                 peak1 not in gene_body_dict[gene1[0]])):
                                corr_dict[gene] = [[peak1], [corr]]
            elif tmp[12] == "gene-peak":
                peak2 = f'{tmp[3]}:{tmp[4]}-{tmp[5]}'
                gene1 = tmp[6].split('><')[0][1:]
                tmp3 = tmp[6].split('><')[1][:-1].split(';')
                corr = float(tmp[7])
                for t3 in tmp3:
                    gene2 = t3.split('_')
                    # peak 2 belongs to gene 1 or are close in distance
                    # -> peak 2 is added to gene 1
                    if ((gene1 == gene2[0]) or (float(tmp[11]) < peak_dist)):
                        gene = gene1
                        if gene in corr_dict:
                            if (peak2 not in corr_dict[gene] and
                                gene2[1] != "promoter" and
                                (gene2[0] not in gene_body_dict or
                                 peak2 not in gene_body_dict[gene2[0]])):

                                corr_dict[gene][0].append(peak2)
                                corr_dict[gene][1].append(corr)
                        else:
                            if (gene2[1] != "promoter" and
                                (gene2[0] not in gene_body_dict or
                                 peak2 not in gene_body_dict[gene2[0]])):

                                corr_dict[gene] = [[peak2], [corr]]

    gene_dict = promoter_dict
    enhancer_dict = {}
    promoter_genes = list(promoter_dict.keys())
    #logg.update(f'Found {len(promoter_genes)} genes with promoter peaks', 1)
    for gene in promoter_genes:
        if gene_body:  # add gene-body peaks
            if gene in gene_body_dict:
                for peak in gene_body_dict[gene]:
                    if peak not in gene_dict[gene]:
                        gene_dict[gene].append(peak)
        enhancer_dict[gene] = []
        if gene in corr_dict:  # add enhancer peaks
            for j, peak in enumerate(corr_dict[gene][0]):
                corr = corr_dict[gene][1][j]
                if corr > min_corr:
                    if peak not in gene_dict[gene]:
                        gene_dict[gene].append(peak)
                        enhancer_dict[gene].append(peak)
    return corr_dict, enhancer_dict, all_info
```

其他的处理函数如下：

```python
def info_display(corr_dict, enhancer_dict, all_info):
    print(f'all peak correlation in linkage_file {len(all_info)}')
    a = 0
    for _ in corr_dict:
        a += len(corr_dict[_])
    print(f'the number of peak-peaks passed distance and correlation threshold {a}')
    
    zero_count = 0
    total_count = 0

    print(f'the number of detected genes {len(enhancer_dict)}')
    
    for _ in enhancer_dict:
        if enhancer_dict[_] == []:
            zero_count += 1
        else:
            total_count += len(enhancer_dict[_])
    print(f'total RE-TG count : {total_count}')
    print(f'the number of genes with no correlated peaks: {zero_count}')
    return None




def info_to_df(all_info):
    tmp = []
    for x in all_info:
        a1 = '_'.join(x[:3])
        a2 = '_'.join(x[3:6])
        peak_info = x[6]
        cor_value = x[7]
        significance = x[10]
        distance = x[11]
        peak_type = x[12]
        tmp.append([a1,a2,peak_info,cor_value,significance, distance, peak_type])
    df = pd.DataFrame(tmp, columns = ['peak_1', 'peak_2', 'peak_info', 'cor_value', 'significance', 'distance', 'peak_type'])
    return df



def get_tg_peaks(enhancer_dict):
    tg_peaks = []
    for _ in enhancer_dict:
        if enhancer_dict[_] == []:
            continue
        else:
            for peak in enhancer_dict[_]:
                tg_peaks.append([_,peak])
    tg_peaks = pd.DataFrame(tg_peaks)
    
    tmp = []
    for _ in tg_peaks.iloc[:,1].values:
        a,x = _.split(':')
        begin,end = x.split('-')
        tmp.append('_'.join([a,begin,end]))
    
    tg_peaks.iloc[:,1] = tmp
    tg_peaks.columns = ['gene','peak']
    return tg_peaks




def peak_tg_pipeline(save_dir, peak_annot_file, linkage_file, 
                 peak_dist , min_corr ):
    
    corr_dict, enhancer_dict, all_info = get_re_tg_10x(
                peak_annot_file, linkage_file, 
                 peak_dist = peak_dist, min_corr = min_corr)
    
    info_display(corr_dict, enhancer_dict, all_info)
    
    df = info_to_df(all_info)
    tg_peaks = get_tg_peaks(enhancer_dict)
    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        
    all_info_path = os.path.join(save_dir,'all_info.csv')
    tg_peaks_path = os.path.join(save_dir, 'tg_peaks.csv')
    region_path = os.path.join(save_dir, 'region.txt')
    df.to_csv(all_info_path)
    tg_peaks.to_csv(tg_peaks_path)
    
    peak_info = tg_peaks.loc[:,'peak'].values 
    a,b,c = [], [], []

    for _ in peak_info:
        x1,x2,x3 = _.split('_')
        a.append(x1)
        b.append(x2)
        c.append(x3) 

    with open(region_path, 'w') as f:
        for i,x1 in enumerate(a):
            f.write(a[i] +'\t' + b[i] + '\t' + c[i] + '\t' + peak_info[i])
            f.write('\n')
        f.close()
        
    print('peak_tg_pipeline is over '  +'===='*15)
    return None
```

在peak_tg_pipeline结束后，可以在save_dir中找到region.txt(记录各个peak的起始位置和结束位置)， tg_peaks.csv（记录TG-peaks），all_info.csv（feature_linkage.bedpe文件的csv格式）

#### TF-peak 

这里看celloracle官方给出的教程，对模式物种应该是没问题的，但自己的数据中有一些非模式物种，因此选择用Homer来构建TF-peak关联。 

Homer 也是一个很经典的软件了，以下是使用命令

```bash
# 执行Perl脚本 findMotifsGenome.pl
if [ ! -d "$homer_output_dir" ]; then
    # 如果目录不存在，使用 `mkdir -p` 来创建目录
    # `-p` 选项允许你创建多级目录
    mkdir -p "$homer_output_dir"
    echo "Directory created: $homer_output_dir"
else
    echo "Directory already exists: $homer_output_dir"
fi

if [ ! -d "$homer_preparsedDir" ]; then
    # 如果目录不存在，使用 `mkdir -p` 来创建目录
    # `-p` 选项允许你创建多级目录
    mkdir -p "$homer_preparsedDir"
    echo "Directory created: $homer_preparsedDir"
else
    echo "Directory already exists: $homer_preparsedDir"
fi

findMotifsGenome.pl "$homer_peaks_file" \
                    "$homer_genome" \
                    "$homer_output_dir" -p 32 -size 200 -mask \
                    -find "$homer_motif_file" \
                    -preparsedDir "$homer_preparsedDir" \
                    > "$homer_output"

# 等待Perl脚本完成
wait
echo "homer task is over"
```

这里参数如下：

|参数|描述|
|---|---|
|homer_peaks_file|前面生成的region.txt文件|
|homer_genome|基因组注释文件, .fa|
|homer_output_dir|输出目录|
|homer_preparsedDir|预处理目录，保存homer运行的中间文件|

Homer运行结束后，调用如下代码生成TF-peak关联：

```python
import numpy as np
import tqdm
import pandas as pd 
import os
import argparse

def get_peaks_tf(homer_bed_path,save_dir,motif2tf_path,K = 1000000):
    # homer_bed_file: the file path of homer output result .bed file
    # save_dir: tf_peaks.csv save path 
    # motif2tf_path: motif_TF.txt file path 
    # K: top K motif result
    
    f = open(homer_bed_path,'r')

    info_list = []
    for i,line in tqdm.tqdm(enumerate(f)):
        if i == 0:
            continue
        info = line.strip().split('\t')
        idx,tf, score = info[0], info[3], eval(info[-1])
        info_list.append([idx,tf,score])
    info = pd.DataFrame(info_list) 


    # 根据最后一列降序排序
    df_sorted = info.sort_values(by=2, ascending=False)

    # 选取Top K的样本
    K = min(K, df_sorted.shape[0])
    top_k_samples = df_sorted.head(K)


    f = open(motif2tf_path,'r') #

    tf_dic = {}
    for line in f:
        tmp = line.strip().split('\t')
        tf_dic[tmp[0]] = tmp[1]
    print(len(tf_dic))

    top_k_samples.loc[:,'tf'] = top_k_samples.iloc[:,1].map(tf_dic).values 
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    top_k_samples.to_csv(os.path.join(save_dir,'tf_peaks.csv'))
    print('get_peaks_tf is over ' + '===='*15)
    return None
```

最终返回 TF_peaks.csv, 记录 TF-peak 关联。

#### TF-TG 

根据有无公共peaks，构建TF-TG关联。下面是我的代码，其中tg_peaks_load里面用到一些gene_map.json文件，这里主要是把非模式物种的基因转成同源基因，方便后续差异分析。最后返回tf_tg.csv和tf_tg_peaks.csv。

```python
def tf_peaks_load(tf_peaks_path):
    tf_peaks = pd.read_csv(tf_peaks_path)
    tf_peaks = tf_peaks.iloc[:,1:]
    tf_peaks = tf_peaks[~tf_peaks.iloc[:,3].isnull()]
    print(f'filtered tf_peaks file shape {tf_peaks.shape}')

    tf_peaks = tf_peaks.iloc[:,[3,0]]

    upper_gene = []
    tmp = tf_peaks.loc[:,'tf'].values
    for _ in tmp:
        upper_gene.append(_.upper())
    tf_peaks.iloc[:,0] = upper_gene
    tf_peaks.columns = ['tf','peak']
    return tf_peaks 

def tg_peaks_load(tg_peaks_path, species):
    tg_peaks = pd.read_csv(tg_peaks_path)
    tg_peaks = tg_peaks.iloc[:,1:]

            
        # for changyi bat genome
    if species == 'CY':
        with open("/home/rsun@ZHANGroup.local/sly_data/genome_transform/mfu_gene_map.json",'r') as f:
            gene_map = json.load(f)
        f.close()
        
    elif species == 'JT':
        # for jutou bat genome
        with open("/home/rsun@ZHANGroup.local/sly_data/genome_transform/rsi_gene_map.json",'r') as f:
            gene_map = json.load(f)
        f.close()

    elif species == 'QF':
        # for quanfu bat genome 
        qf_map=pd.read_excel("/home/rsun@ZHANGroup.local/sly_data/genome_transform/csp_gene.xlsx",sheet_name='CSP')
        gene_map = {}
        for i in range(qf_map.shape[0]):
            gene_map[qf_map.iloc[i,1]] = qf_map.iloc[i,0].upper()
    elif species == 'M':
        with open("/home/rsun@ZHANGroup.local/sly_data/genome_transform/m_gene_map.json",'r') as f:
            gene_map = json.load(f)
        f.close()
    elif species =='T':
        with open("/home/rsun@ZHANGroup.local/sly_data/genome_transform/tda_gene_map.json",'r') as f:
            gene_map = json.load(f)
        f.close()
    else:
        print('Wrong Species')

    tg_peaks.loc[:,'gene_symbol'] = tg_peaks.iloc[:,0].map(gene_map).values
    print(f'filtered tg_peaks file shape {tg_peaks.shape}')
    return tg_peaks, gene_map

def joint_process(tg_peaks, tf_peaks,gene_map, save_dir):
    uni_tf_p = tf_peaks.loc[:,'peak'].unique()
    uni_tg_p = tg_peaks.loc[:,'peak'].unique()

    uni_p = np.intersect1d(uni_tf_p, uni_tg_p)
    print(f' unique peaks in tf_peaks {uni_tf_p.shape}')
    print(f' unique peaks in tg_peaks {uni_tg_p.shape}')
    print(f' mutual peaks {uni_p.shape}')
    
    
    tf_genes  = tf_peaks.loc[:,'tf'].unique()

    mapped_genes = list(gene_map.values())
    mapped_genes = np.unique(np.array(mapped_genes))

    mutual_genes = np.intersect1d(tf_genes, mapped_genes)
    print(f'unique genes in tf_genes {tf_genes.shape}')
    print(f'unique genes in gene_map {mapped_genes.shape}')
    print(f'mutual genes {mutual_genes.shape}')

    ### get tf_tg_peak_info
    tf_tg_peak = []

    for peak in tqdm.tqdm(uni_p):
        sel_tf_id = tf_peaks.loc[:,'peak'] == peak
        sel_tg_id = tg_peaks.loc[:,'peak'] == peak 

        sel_tf = tf_peaks[sel_tf_id].loc[:,'tf'].unique()
        sel_tg = tg_peaks[sel_tg_id].loc[:,'gene_symbol'].unique()

        for tf in sel_tf:
            for tg in sel_tg:
                tf_tg_peak.append([tf,tg,peak])
    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    tf_tg_peak = pd.DataFrame(tf_tg_peak)
    tf_tg_peak.columns = ['tf', 'tg', 'peak']
    print(f' tf_tg_peaks num {tf_tg_peak.shape}')
    tf_tg_peak.to_csv(os.path.join(save_dir, 'tf_tg_peak.csv'))
    
    
    tf_tg = tf_tg_peak.iloc[:,:2]
    tf_tg = tf_tg.drop_duplicates()
    print(f'tf_tg edge num {tf_tg.shape} tf num {tf_tg.iloc[:,0].nunique()} tg num {tf_tg.iloc[:,1].nunique()}' )

    tf_tg.to_csv(os.path.join(save_dir, 'tf_tg.csv'))
    return None 
```

完成上面三个步骤后，得到base GRN(tf-tg 关联)。


### GRN inference  

GRN inference 可以follow celloracle的教程。下面是构建oracle object时，怎样应用前面得到的base GRN的代码

```python
def load_base_grn(base_grn_path):
    tf_tg = pd.read_csv(base_grn_path)
    tf_tg = tf_tg.iloc[:,1:]
    
    tf_tg = tf_tg[~tf_tg.iloc[:,1].isnull()]
    tf_tg = tf_tg[~tf_tg.iloc[:,0].isnull()]
    
    tf_set = tf_tg.iloc[:,0].unique()
    tf_tg_dic = {}

    for tf in tf_set:
        sub_df = tf_tg[tf_tg.iloc[:,0] == tf]
        tf_tg_dic[tf] = sub_df.iloc[:,1].unique()

    tg_tf_dic = co.utility.inverse_dictionary(tf_tg_dic)
    return tf_tg,tf_tg_dic, tg_tf_dic

def inference(scdata,tg_tf_dic, save_dir):

    os.makedirs(save_dir, exist_ok=True)
    # Instantiate Oracle object
    oracle = co.Oracle()

    # Instantiate Oracle object.
    oracle.import_anndata_as_raw_count(adata=scdata,
                                    cluster_column_name="new_anno",
                                    embedding_name="X_umap")
    
    new_tg_tf_dic = {}
    for ele in tg_tf_dic:
        tmp = [key.capitalize() for key in tg_tf_dic[ele]]
        new_tg_tf_dic[ele.capitalize()] = tmp

    # Add TF information 
    oracle.addTFinfo_dictionary(new_tg_tf_dic)
        # Perform PCA
    oracle.perform_PCA()
    # Select important PCs and imputation
    n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
    n_comps = min(n_comps, 50)
    n_cell = oracle.adata.shape[0]
    k = int(0.025*n_cell)
    oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                        b_maxl=k*4, n_jobs=4)
    

    # Calculate GRN for each population in "louvain_annot" clustering unit.
    # This step may take some time.(~30 minutes)
    links = oracle.get_links(cluster_name_for_GRN_unit="new_anno", 
                            alpha=2.5,
                            verbose_level=10)

    oracle.to_hdf5(os.path.join(save_dir,"grn.celloracle.oracle"))
    links.to_hdf5(os.path.join(save_dir, "grn.celloracle.links"))
    return None
```

以上就是celloracle GRN推断的全部内容。