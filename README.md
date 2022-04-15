# MapRRCon - Sun18

ChIP-seq analysis on repetitive elements.

## Database of TF binders on LINE-1

MapRRCon database is available at: https://maprrcon.org

## Reference

The methodology and analyzed data is described in the following paper:

Xiaoji Sun, Xuya Wang, Zuojian Tang, Mark Grivainis, David Kahler, Chi Yun, Paolo Mita, David Fenyö, Jef D. Boeke. Transcription factor profiling reveals molecular choreography and key regulators of human retrotransposon expression. Proc Natl Acad Sci U S A. 2018 May 25. pii: 201722565. doi: 10.1073/pnas.1722565115. 
https://www.ncbi.nlm.nih.gov/pubmed/29802231

**Abstract**
Transposable elements (TEs) represent a substantial fraction of many eukaryotic genomes, and transcriptional regulation of these factors is important to determine TE activities in human cells. However, due to the repetitive nature of TEs, identifying transcription factor (TF)-binding sites from ChIP-sequencing (ChIP-seq) datasets is challenging. Current algorithms are focused on subtle differences between TE copies and thus bias the analysis to relatively old and inactive TEs. Here we describe an approach termed “MapRRCon” (mapping repeat reads to a consensus) which allows us to identify proteins binding to TE DNA sequences by mapping ChIP-seq reads to the TE consensus sequence after whole-genome alignment. Although this method does not assign binding sites to individual insertions in the genome, it provides a landscape of interacting TFs by capturing factors that bind to TEs under various conditions. We applied this method to screen TFs’ interaction with L1 in human cells/tissues using ENCODE ChIP-seq datasets and identified 178 of the 512 TFs tested as bound to L1 in at least one biological condition with most of them (138) localized to the promoter. Among these L1-binding factors, we focused on Myc and CTCF, as they play important roles in cancer progression and 3D chromatin structure formation. Furthermore, we explored the transcriptomes of The Cancer Genome Atlas breast and ovarian tumor samples in which a consistent anti-/correlation between L1 and Myc/CTCF expression was observed, suggesting that these two factors may play roles in regulating L1 transcription during the development of such tumors.

**Significance**
Retrotransposons replicate through RNA intermediates that are reverse transcribed and inserted at new genomic locations. LINE-1 (L1) elements constitute ∼17% of the human genome, making them the most successful retrotransposons in the human genome by mass. The activity of L1s was shown first in the germline or during early embryogenesis. More recent studies demonstrate a wider prevalence of L1 expression in somatic cells including neurons, aging cells, and different types of cancer. In this study, we developed the MapRRCon pipeline and performed a comprehensive computational analysis of L1 transcriptional regulators using ENCODE ChIP-seq datasets. We revealed the binding of various transcription factors, including Myc and CTCF, to the 5′ UTR promoter of the youngest human L1 family (L1HS) and their potential functional impact on L1HS expression.

