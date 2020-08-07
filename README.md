# JAMSE
## (Just Another Molecular Search Engine)

This repository contains the proof-of-concept embodiment of the theory described in the publication referenced below. No descriptions or documentation are given for its implementation at this time. Please contact the manuscript authors.   

### A Pre-computed Probabilistic Molecular Search Engine for Tandem Mass Spectrometry Proteomics

[https://www.biorxiv.org/content/10.1101/2020.02.06.937870v2](https://www.biorxiv.org/content/10.1101/2020.02.06.937870v2)

Mass spectrometry methods of peptide identification involve comparing observed tandem spectra with in-silico derived spectrum models. Presented here is a proteomics search engine that offers a new variation of the standard approach, with improved results. The proposed method employs information theory and probabilistic information retrieval on a pre-computed and indexed fragmentation database generating a peptide-to-spectrum match (PSM) score modeled on fragment ion frequency. As a result, the direct application of modern document mining, allows for treating the collection of peptides as a corpus and corresponding fragment ions as indexable words, leveraging ready-built search engines and common predefined ranking algorithms. Fast and accurate PSM matches are achieved yielding a 5-10% higher rate of peptide identities than current database mining methods. Immediate applications of this search engine are aimed at identifying peptides from large sequence databases consisting of homologous proteins with minor sequence variations, such as genetic variation expected in the human population.
