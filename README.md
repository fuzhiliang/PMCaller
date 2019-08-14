# PMCaller
此程序用于低深度胚系突变的calling，以及计算数据本身的一些特征的脚本 

方法：
对SNP和InDel calling采用贝叶斯模型，通过数据本身的关系对先验概率进行校正，计算不同genotype的likelihood，取后验概率最大的genotype。

优势： 
对比GATK的先验概率是一个经验值，该方法对数据本身进行了计算，得到一个更准确的先验概率，对计算后验概率，理论上也更准确。

 
